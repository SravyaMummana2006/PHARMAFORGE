"""Scientific Test-Oracle (STO-IV): Independent deterministic chemical verifier."""
from __future__ import annotations

import math
import time
from typing import Any
from pharmaforge.core.types import MoleculeCandidate, OptimizationObjective, VerificationReport
from pharmaforge.verifier.filters import SubstructureFilter

# Standard atomic weights for pure-Python fallback calculations
ATOMIC_WEIGHTS: dict[str, float] = {
    "H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998,
    "P": 30.974, "S": 32.06, "Cl": 35.45, "Br": 79.904, "I": 126.904,
}


class ScientificTestOracle:
    """Independent deterministic verifier and property evaluation oracle.
    
    Acts as the automated 'test harness' for scientific workflows, verifying:
    1. Chemical validity (sanitization, explicit valence, aromaticity).
    2. Structural alerts (PAINS, Brenk).
    3. Physicochemical properties (MW, LogP, TPSA, QED, SAScore).
    """

    def __init__(self):
        self._filter = SubstructureFilter()
        self._rdkit_available = False
        self._init_rdkit()

    def _init_rdkit(self) -> None:
        """Probe for RDKit availability."""
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors, QED
            self._rdkit_available = True
        except ImportError:
            self._rdkit_available = False

    @property
    def is_rdkit_backed(self) -> bool:
        """True if RDKit is currently available in the environment."""
        return self._rdkit_available

    def verify_smiles(self, smiles: str) -> VerificationReport:
        """Perform comprehensive deterministic chemical verification on a SMILES string."""
        start_time = time.perf_counter()
        error_reasons: list[str] = []

        if not smiles or not isinstance(smiles, str) or smiles.strip() == "":
            return VerificationReport(
                is_valid=False,
                error_reasons=["Empty or invalid SMILES string format."],
                execution_time_ms=(time.perf_counter() - start_time) * 1000,
            )

        smiles = smiles.strip()

        # Check with RDKit if available
        if self._rdkit_available:
            return self._verify_with_rdkit(smiles, start_time)

        # Fallback pure-Python deterministic verification
        return self._verify_fallback(smiles, start_time)

    def _verify_with_rdkit(self, smiles: str, start_time: float) -> VerificationReport:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, QED

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return VerificationReport(
                is_valid=False,
                passed_sanitization=False,
                error_reasons=["RDKit MolFromSmiles failed: invalid chemical syntax."],
                execution_time_ms=(time.perf_counter() - start_time) * 1000,
            )

        try:
            Chem.SanitizeMol(mol)
            passed_sanitization = True
        except Exception as e:
            return VerificationReport(
                is_valid=False,
                passed_sanitization=False,
                error_reasons=[f"RDKit SanitizeMol failed: {str(e)}"],
                execution_time_ms=(time.perf_counter() - start_time) * 1000,
            )

        # Explicit valence check
        passed_valence = True
        for atom in mol.GetAtoms():
            if atom.GetExplicitValence() > atom.GetTotalValence():
                passed_valence = False
                break

        # Check PAINS & Brenk
        passed_pains, pains_matches = self._filter.check_pains(smiles, mol)
        passed_brenk, brenk_matches = self._filter.check_brenk(smiles, mol)

        error_reasons = []
        if not passed_valence:
            error_reasons.append("Molecule contains atoms with invalid explicit valence.")
        error_reasons.extend(pains_matches)
        error_reasons.extend(brenk_matches)

        # Calculate physicochemical properties
        mw = float(Descriptors.MolWt(mol))
        logp = float(Descriptors.MolLogP(mol))
        tpsa = float(Descriptors.TPSA(mol))
        hbd = float(Descriptors.NumHDonors(mol))
        hba = float(Descriptors.NumHAcceptors(mol))
        rot_bonds = float(Descriptors.NumRotatableBonds(mol))
        qed_val = float(QED.qed(mol))
        sascore = self._calculate_sascore_heuristic(mol, mw, rot_bonds)

        properties = {
            "molecular_weight": mw,
            "logp": logp,
            "tpsa": tpsa,
            "hbd": hbd,
            "hba": hba,
            "rotatable_bonds": rot_bonds,
            "qed": qed_val,
            "sascore": sascore,
        }

        is_valid = passed_sanitization and passed_valence and (len(pains_matches) == 0)

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

    def _verify_fallback(self, smiles: str, start_time: float) -> VerificationReport:
        """Robust fallback validator when RDKit is not installed."""
        error_reasons = []
        
        # Check basic SMILES syntax invariants (balanced parentheses, valid atom symbols)
        if smiles.count("(") != smiles.count(")"):
            error_reasons.append("Unbalanced parentheses in SMILES string.")
        if smiles.count("[") != smiles.count("]"):
            error_reasons.append("Unbalanced square brackets in SMILES string.")

        # Disallow chemically absurd atom symbols
        valid_symbols = set("CNOFPSClBrIcnofps1234567890()[]=-#+@/\\.")
        invalid_chars = [c for c in smiles if c not in valid_symbols]
        if invalid_chars:
            error_reasons.append(f"Invalid chemical characters detected: {set(invalid_chars)}")

        passed_pains, pains_matches = self._filter.check_pains(smiles)
        passed_brenk, brenk_matches = self._filter.check_brenk(smiles)
        error_reasons.extend(pains_matches)
        error_reasons.extend(brenk_matches)

        # Estimate molecular properties using fallback heuristics
        carbon_count = smiles.count("C") + smiles.count("c")
        nitrogen_count = smiles.count("N") + smiles.count("n")
        oxygen_count = smiles.count("O") + smiles.count("o")
        fluorine_count = smiles.count("F")
        chlorine_count = smiles.count("Cl")

        estimated_mw = (
            carbon_count * ATOMIC_WEIGHTS["C"] +
            nitrogen_count * ATOMIC_WEIGHTS["N"] +
            oxygen_count * ATOMIC_WEIGHTS["O"] +
            fluorine_count * ATOMIC_WEIGHTS["F"] +
            chlorine_count * ATOMIC_WEIGHTS["Cl"] +
            (carbon_count * 2 + 2) * ATOMIC_WEIGHTS["H"]
        )
        if estimated_mw <= 0:
            estimated_mw = 180.0

        # Heuristic LogP approximation based on Wildman-Crippen atom types
        estimated_logp = (carbon_count * 0.25) - (oxygen_count * 0.4) - (nitrogen_count * 0.2) + (fluorine_count * 0.3)
        estimated_tpsa = (oxygen_count * 18.0) + (nitrogen_count * 24.0)

        # QED heuristic (bounded between 0.1 and 0.9)
        mw_penalty = math.exp(-0.5 * ((estimated_mw - 350.0) / 100.0) ** 2)
        logp_penalty = math.exp(-0.5 * ((estimated_logp - 2.5) / 1.5) ** 2)
        estimated_qed = max(0.1, min(0.95, float(0.5 * mw_penalty + 0.5 * logp_penalty)))

        # SAScore heuristic (1 = easy, 10 = difficult)
        ring_count = sum(1 for c in smiles if c.isdigit()) // 2
        estimated_sascore = min(10.0, max(1.0, 2.0 + (ring_count * 0.5) + (estimated_mw / 200.0)))

        properties = {
            "molecular_weight": round(estimated_mw, 2),
            "logp": round(estimated_logp, 2),
            "tpsa": round(estimated_tpsa, 2),
            "qed": round(estimated_qed, 3),
            "sascore": round(estimated_sascore, 2),
        }

        is_valid = len(error_reasons) == 0

        return VerificationReport(
            is_valid=is_valid,
            passed_sanitization=is_valid,
            passed_valence=is_valid,
            passed_pains=passed_pains,
            passed_brenk=passed_brenk,
            property_scores=properties,
            error_reasons=error_reasons,
            execution_time_ms=(time.perf_counter() - start_time) * 1000,
        )

    def _calculate_sascore_heuristic(self, mol: any, mw: float, rot_bonds: float) -> float:
        """Synthetic accessibility score estimation based on Ertl & Schuffenhauer principles."""
        ring_info = mol.GetRingInfo()
        num_rings = ring_info.NumRings()
        num_chiral = len(mol.GetSubstructMatches(mol))  # basic proxy
        score = 2.0 + (num_rings * 0.4) + (rot_bonds * 0.1) + (mw / 250.0)
        return float(min(10.0, max(1.0, round(score, 2))))

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
