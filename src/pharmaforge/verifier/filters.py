"""Substructure filters and alert catalogs for PAINS and Brenk toxicophores."""
from __future__ import annotations

# Representative SMARTS patterns for common PAINS (Pan-Assay Interference Compounds)
# Reference: Baell & Holloway, J. Med. Chem. 2010
PAINS_SMARTS: dict[str, str] = {
    "quinone": "[#6]1(=[O,S])[#6]=[#6][#6](=[O,S])[#6]=[#6]1",
    "rhodanine": "O=C1NC(=S)SC1=*",
    "ene_rhodanine": "O=C1NC(=S)S/C1=C",
    "catechol": "c1c(O)c(O)ccc1",
    "azo": "[#6]-[#7]=[#7]-[#6]",
    "alkyl_halide_reactive": "[CX4][Cl,Br,I]",
    "michael_acceptor": "[#6]=[#6]-[#6](=[O,S])",
    "thioester": "[#6]-[#6](=[O,S])-[S]-[#6]",
    "disulfide": "[#16]-[#16]",
    "epoxide": "C1OC1",
}

# Representative Brenk structural alerts (Brenk et al., ChemMedChem 2008)
BRENK_SMARTS: dict[str, str] = {
    "nitro_aromatic": "c[N+](=O)[O-]",
    "aliphatic_ester": "[CX4][C](=O)O[CX4]",
    "thiocarbonyl": "[#6](=[S])",
    "hydrazine": "[#7]-[#7]",
    "peroxide": "[#8]-[#8]",
    "acid_halide": "[CX3](=[OX1])[F,Cl,Br,I]",
    "sulfonyl_halide": "S(=O)(=O)[F,Cl,Br,I]",
    "isocyanide": "[C-]#[N+]",
}


class SubstructureFilter:
    """Evaluates molecules against toxicophore and interference alert catalogs."""

    def __init__(self):
        self._rdkit_available = False
        self._pains_mols = {}
        self._brenk_mols = {}
        self._initialize_oracles()

    def _initialize_oracles(self) -> None:
        """Attempt to compile SMARTS patterns using RDKit."""
        try:
            from rdkit import Chem
            self._rdkit_available = True
            for name, smarts in PAINS_SMARTS.items():
                mol = Chem.MolFromSmarts(smarts)
                if mol is not None:
                    self._pains_mols[name] = mol
            for name, smarts in BRENK_SMARTS.items():
                mol = Chem.MolFromSmarts(smarts)
                if mol is not None:
                    self._brenk_mols[name] = mol
        except ImportError:
            self._rdkit_available = False

    def check_pains(self, smiles: str, mol: any = None) -> tuple[bool, list[str]]:
        """Check if molecule triggers PAINS alerts.
        
        Returns:
            Tuple of (passed: bool, matched_alerts: list[str])
        """
        matched = []
        if self._rdkit_available:
            from rdkit import Chem
            if mol is None:
                mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                for name, pattern in self._pains_mols.items():
                    if mol.HasSubstructMatch(pattern):
                        matched.append(f"PAINS:{name}")
                return len(matched) == 0, matched

        # Heuristic fallback if RDKit is not installed
        # Check simple characteristic string fragments
        if "N=N" in smiles or "[N+]#[C-]" in smiles:
            matched.append("PAINS:azo_or_isocyanide_heuristic")
        if "C1OC1" in smiles:
            matched.append("PAINS:epoxide_heuristic")
        return len(matched) == 0, matched

    def check_brenk(self, smiles: str, mol: any = None) -> tuple[bool, list[str]]:
        """Check if molecule triggers Brenk structural alerts.
        
        Returns:
            Tuple of (passed: bool, matched_alerts: list[str])
        """
        matched = []
        if self._rdkit_available:
            from rdkit import Chem
            if mol is None:
                mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                for name, pattern in self._brenk_mols.items():
                    if mol.HasSubstructMatch(pattern):
                        matched.append(f"Brenk:{name}")
                return len(matched) == 0, matched

        # Heuristic fallback for Brenk
        if "[N+](=O)[O-]" in smiles:
            matched.append("Brenk:nitro_heuristic")
        if "OO" in smiles:
            matched.append("Brenk:peroxide_heuristic")
        return len(matched) == 0, matched
