"""Substructure filters and alert catalogs using RDKit FilterCatalog and verified SMARTS patterns."""
from __future__ import annotations

from pharmaforge.core.exceptions import RDKitUnavailableError

# Supplementary medicinal chemistry structural alerts
# Documented alert patterns for reactive or interfering substructures
SUPPLEMENTARY_PAINS_SMARTS: dict[str, str] = {
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

SUPPLEMENTARY_BRENK_SMARTS: dict[str, str] = {
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
    """Evaluates molecules against official RDKit FilterCatalogs (PAINS, Brenk) and supplementary SMARTS.
    
    Strictly requires RDKit. Does NOT perform silent heuristic or string-based fallback.
    """

    def __init__(self):
        try:
            from rdkit import Chem
            from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
            self._Chem = Chem
        except ImportError as e:
            raise RDKitUnavailableError("SubstructureFilter initialization (FilterCatalog)") from e

        # 1. Initialize RDKit Official PAINS Catalog (PAINS_A, PAINS_B, PAINS_C)
        pains_params = FilterCatalogParams()
        pains_params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
        pains_params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
        pains_params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
        self._pains_catalog = FilterCatalog(pains_params)

        # 2. Initialize RDKit Official Brenk Catalog
        brenk_params = FilterCatalogParams()
        brenk_params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
        self._brenk_catalog = FilterCatalog(brenk_params)

        # 3. Initialize Supplementary Curated SMARTS Patterns
        self._supp_pains = {}
        for name, smarts in SUPPLEMENTARY_PAINS_SMARTS.items():
            pattern = Chem.MolFromSmarts(smarts)
            if pattern is not None:
                self._supp_pains[name] = pattern

        self._supp_brenk = {}
        for name, smarts in SUPPLEMENTARY_BRENK_SMARTS.items():
            pattern = Chem.MolFromSmarts(smarts)
            if pattern is not None:
                self._supp_brenk[name] = pattern

    def _ensure_mol(self, mol_or_smiles: any) -> any:
        """Helper to ensure input is a valid RDKit Mol object."""
        if isinstance(mol_or_smiles, str):
            mol = self._Chem.MolFromSmiles(mol_or_smiles)
            if mol is None:
                raise ValueError(f"Cannot parse invalid SMILES for substructure filtering: {mol_or_smiles}")
            return mol
        return mol_or_smiles

    def check_pains(self, mol_or_smiles: any) -> tuple[bool, list[str]]:
        """Check if molecule triggers PAINS alerts via RDKit FilterCatalog and supplementary SMARTS.
        
        Returns:
            Tuple of (passed: bool, matched_alerts: list[str])
        """
        mol = self._ensure_mol(mol_or_smiles)
        matched: list[str] = []

        # Official RDKit FilterCatalog matches
        for entry in self._pains_catalog.GetMatches(mol):
            matched.append(f"PAINS_Catalog:{entry.GetDescription()}")

        # Supplementary SMARTS pattern matches
        for name, pattern in self._supp_pains.items():
            if mol.HasSubstructMatch(pattern):
                alert = f"PAINS_Supp:{name}"
                if alert not in matched:
                    matched.append(alert)

        return len(matched) == 0, matched

    def check_brenk(self, mol_or_smiles: any) -> tuple[bool, list[str]]:
        """Check if molecule triggers Brenk structural alerts via RDKit FilterCatalog and supplementary SMARTS.
        
        Returns:
            Tuple of (passed: bool, matched_alerts: list[str])
        """
        mol = self._ensure_mol(mol_or_smiles)
        matched: list[str] = []

        # Official RDKit FilterCatalog matches
        for entry in self._brenk_catalog.GetMatches(mol):
            matched.append(f"Brenk_Catalog:{entry.GetDescription()}")

        # Supplementary SMARTS pattern matches
        for name, pattern in self._supp_brenk.items():
            if mol.HasSubstructMatch(pattern):
                alert = f"Brenk_Supp:{name}"
                if alert not in matched:
                    matched.append(alert)

        return len(matched) == 0, matched
