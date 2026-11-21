"""Scientific Test-Oracle (STO-IV): Independent deterministic chemical verifier using RDKit."""
from __future__ import annotations

import time
from typing import Any
from pharmaforge.core.exceptions import RDKitUnavailableError
from pharmaforge.core.types import MoleculeCandidate, VerificationReport
from pharmaforge.verifier.filters import SubstructureFilter


class ScientificTestOracle:
    """Independent deterministic verifier and property evaluation oracle backed by RDKit.
    
    Strictly performs genuine chemoinformatics calculations. Does NOT silently fall back
    to non-scientific character-counting heuristics.
    """

    def __init__(self):
        self._init_rdkit()
        self._filter = SubstructureFilter()

    def _init_rdkit(self) -> None:
        """Validate that all mandatory RDKit modules are importable."""
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors, QED, FilterCatalog, AllChem
            self._Chem = Chem
            self._Descriptors = Descriptors
            self._QED = QED
            self._AllChem = AllChem
            self._FilterCatalog = FilterCatalog
        except ImportError as e:
            raise RDKitUnavailableError("ScientificTestOracle initialization") from e

        # Probe for the genuine Ertl SAScore module in rdkit.Contrib.SA_Score
        self._sascorer = None
        try:
            from rdkit.Contrib.SA_Score import sascorer
            self._sascorer = sascorer
        except Exception:
            self._sascorer = None

    @property
    def is_rdkit_backed(self) -> bool:
        """True if RDKit is currently verified and available."""
        return True

    @property
    def is_sascore_available(self) -> bool:
        """True if the genuine Ertl SAScore module is loaded with its fragment database."""
        return self._sascorer is not None

    def verify_smiles(self, smiles: str) -> VerificationReport:
        """Perform comprehensive deterministic chemical verification on a SMILES string.
        
        Evaluates:
        1. SMILES syntactic parsing.
        2. RDKit sanitization (valence, aromaticity, kekulization).
        3. Substructure alerts (official RDKit FilterCatalog for PAINS and Brenk).
        4. Exact physicochemical descriptors (MW, LogP, TPSA, QED, and authentic SAScore).
        
        Raises:
            RDKitUnavailableError: If invoked without RDKit installed.
        """
        start_time = time.perf_counter()
        error_reasons: list[str] = []

        if not smiles or not isinstance(smiles, str) or smiles.strip() == "":
            return VerificationReport(
                is_valid=False,
                error_reasons=["Empty or invalid SMILES string format."],
                execution_time_ms=(time.perf_counter() - start_time) * 1000,
            )

        smiles = smiles.strip()

        # Step 1: Parse SMILES without sanitization to isolate syntax errors
        mol = self._Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return VerificationReport(
                is_valid=False,
                passed_sanitization=False,
                error_reasons=["RDKit MolFromSmiles failed: invalid chemical syntax."],
                execution_time_ms=(time.perf_counter() - start_time) * 1000,
            )

        # Step 2: Sanitize molecular structure (enforces valence and aromaticity)
        try:
            self._Chem.SanitizeMol(mol)
            passed_sanitization = True
        except Exception as e:
            return VerificationReport(
                is_valid=False,
                passed_sanitization=False,
                error_reasons=[f"RDKit SanitizeMol rejected molecule: {str(e)}"],
                execution_time_ms=(time.perf_counter() - start_time) * 1000,
            )

        # Step 3: Explicit valence inspection
        passed_valence = True
        for atom in mol.GetAtoms():
            try:
                explicit_v = atom.GetValence(self._Chem.ValenceType.EXPLICIT)
            except (AttributeError, TypeError):
                explicit_v = atom.GetExplicitValence()

            if explicit_v > atom.GetTotalValence():
                passed_valence = False
                error_reasons.append(
                    f"Atom {atom.GetSymbol()}(idx={atom.GetIdx()}) has invalid explicit valence "
                    f"({explicit_v} > {atom.GetTotalValence()})."
                )
                break

        # Step 4: Substructure FilterCatalog matching (PAINS & Brenk)
        passed_pains, pains_matches = self._filter.check_pains(mol)
        passed_brenk, brenk_matches = self._filter.check_brenk(mol)

        error_reasons.extend(pains_matches)
        error_reasons.extend(brenk_matches)

        # Step 5: Exact physicochemical property calculation
        mw = float(self._Descriptors.MolWt(mol))
        logp = float(self._Descriptors.MolLogP(mol))
        tpsa = float(self._Descriptors.TPSA(mol))
        hbd = float(self._Descriptors.NumHDonors(mol))
        hba = float(self._Descriptors.NumHAcceptors(mol))
        rot_bonds = float(self._Descriptors.NumRotatableBonds(mol))
        qed_val = float(self._QED.qed(mol))
        sascore = self._calculate_sascore(mol)

        properties: dict[str, float] = {
            "molecular_weight": round(mw, 4),
            "logp": round(logp, 4),
            "tpsa": round(tpsa, 4),
            "hbd": float(hbd),
            "hba": float(hba),
            "rotatable_bonds": float(rot_bonds),
            "qed": round(qed_val, 4),
        }
        if sascore is not None:
            properties["sascore"] = round(sascore, 4)

        # Chemical validity requires valid sanitization, valid valence, and no PAINS assay interference
        is_valid = passed_sanitization and passed_valence and passed_pains

        return VerificationReport(
            is_valid=is_valid,
            passed_sanitization=passed_sanitization,
            passed_valence=passed_valence,
            passed_pains=passed_pains,
            passed_brenk=passed_brenk,
            property_scores=properties,
            error_reasons=error_reasons,
            execution_time_ms=(time.perf_counter() - start_time) * 1000,
        )

    def _calculate_sascore(self, mol: any) -> float | None:
        """Calculate genuine Ertl & Schuffenhauer Synthetic Accessibility score.
        
        If the established RDKit contrib sascorer module or its fragment database
        is unavailable, explicitly returns None. Does NOT invent a heuristic formula.
        """
        if self._sascorer is not None:
            try:
                score = float(self._sascorer.calculateScore(mol))
                return float(min(10.0, max(1.0, score)))
            except Exception:
                return None
        return None

    def generate_morgan_fingerprint(
        self,
        smiles_or_mol: str | any,
        radius: int = 2,
        n_bits: int = 2048,
    ) -> list[int]:
        """Generate canonical Morgan fingerprint (ECFP4 equivalent for radius=2) as bit list."""
        if isinstance(smiles_or_mol, str):
            mol = self._Chem.MolFromSmiles(smiles_or_mol)
            if mol is None:
                raise ValueError(f"Cannot generate fingerprint for invalid SMILES: {smiles_or_mol}")
        else:
            mol = smiles_or_mol

        try:
            from rdkit.Chem import rdFingerprintGenerator
            gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
            fp = gen.GetFingerprint(mol)
        except (ImportError, AttributeError):
            fp = self._AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)

        # Convert bit vector to standard Python integer list
        return [int(b) for b in fp.ToBitString()]

    def execute_reaction(
        self,
        reaction_smarts: str,
        reactant_smiles: list[str],
    ) -> list[str]:
        """Execute and validate a chemical reaction via Reaction SMARTS.
        
        Validates:
        1. Reaction SMARTS parsing.
        2. Reactant structure and count matching.
        3. Product generation and RDKit sanitization.
        4. Canonical SMILES conversion.
        
        Returns:
            List of valid, sanitized product SMILES strings.
        """
        rxn = self._AllChem.ReactionFromSmarts(reaction_smarts)
        if rxn is None:
            raise ValueError(f"Invalid reaction SMARTS pattern: {reaction_smarts}")
        rxn.Initialize()

        expected_reactants = rxn.GetNumReactantTemplates()
        if len(reactant_smiles) != expected_reactants:
            raise ValueError(
                f"Reaction requires {expected_reactants} reactants, but {len(reactant_smiles)} were provided."
            )

        reactants = []
        for s in reactant_smiles:
            mol = self._Chem.MolFromSmiles(s)
            if mol is None:
                raise ValueError(f"Invalid reactant SMILES: {s}")
            self._Chem.SanitizeMol(mol)
            if not rxn.IsMoleculeReactant(mol):
                raise ValueError(f"Reactant molecule '{s}' does not match any reaction template in SMARTS.")
            reactants.append(mol)

        product_sets = rxn.RunReactants(tuple(reactants))
        sanitized_products: list[str] = []

        for prod_tuple in product_sets:
            for prod_mol in prod_tuple:
                try:
                    self._Chem.SanitizeMol(prod_mol)
                    prod_smiles = self._Chem.MolToSmiles(prod_mol, canonical=True)
                    if prod_smiles not in sanitized_products:
                        sanitized_products.append(prod_smiles)
                except Exception:
                    continue

        return sanitized_products

    def create_candidate(
        self,
        smiles: str,
        parent_smiles: str | None = None,
        transformation: str | None = None,
        affinity_score: float | None = None,
    ) -> MoleculeCandidate:
        """Validate SMILES and wrap into an immutable MoleculeCandidate entity."""
        report = self.verify_smiles(smiles)
        if not report.is_valid:
            raise ValueError(f"Cannot create candidate from invalid SMILES '{smiles}': {report.error_reasons}")

        props = report.property_scores
        return MoleculeCandidate(
            smiles=smiles,
            molecular_weight=props.get("molecular_weight"),
            logp=props.get("logp"),
            tpsa=props.get("tpsa"),
            qed=props.get("qed"),
            sascore=props.get("sascore"),
            affinity_score=affinity_score,
            parent_smiles=parent_smiles,
            transformation_applied=transformation,
            metadata={"verification_report": report.model_dump()},
        )
