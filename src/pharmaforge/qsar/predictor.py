"""Integrated EGFR QSAR Predictor Interface for PHARMAFORGE.

Provides:
- In silico surrogate potency estimation (predicted pIC50).
- Heuristic ensemble predictive dispersion (sigma_trees).
- Structural applicability domain diagnostic (nearest-neighbor Tanimoto T_max).
- Explicit epistemic caveats (never claims experimental validation or clinical efficacy).
"""

from __future__ import annotations
import os
import joblib
from pathlib import Path
from typing import Any
import numpy as np
from rdkit import Chem

from pharmaforge.core.types import MoleculeCandidate
from pharmaforge.qsar.features import MolecularFeaturizer
from pharmaforge.qsar.evaluator import compute_max_tanimoto

DEFAULT_MODEL_PATH = Path("artifacts/qsar/egfr_qsar_best_model.joblib")
DEFAULT_FEATURES_PATH = Path("data/qsar_precomputed_features.npz")


class EGFRQSARPredictor:
    """Surrogate affinity prediction engine integrated into the PHARMAFORGE multi-agent system."""

    def __init__(
        self,
        model_path: str | Path | None = None,
        reference_features_path: str | Path | None = None,
    ):
        self.model_path = Path(model_path) if model_path else DEFAULT_MODEL_PATH
        self.reference_features_path = Path(reference_features_path) if reference_features_path else DEFAULT_FEATURES_PATH

        self.featurizer = MolecularFeaturizer(radius=2, n_bits=2048)
        self.model = None
        self._ref_morgan: np.ndarray | None = None

        if self.model_path.exists():
            self._load_model()
        if self.reference_features_path.exists():
            self._load_reference_data()

    def _load_model(self) -> None:
        """Load serialized model artifact."""
        self.model = joblib.load(self.model_path)

    def _load_reference_data(self) -> None:
        """Load reference training fingerprints to evaluate structural applicability domain."""
        try:
            data = np.load(self.reference_features_path)
            if "X_m_core" in data:
                self._ref_morgan = data["X_m_core"]
        except Exception:
            self._ref_morgan = None

    @property
    def is_ready(self) -> bool:
        """True if model is loaded and available for inference."""
        return self.model is not None

    def predict(self, smiles: str) -> dict[str, Any]:
        """Generate predicted EGFR pIC50 with uncertainty and applicability domain telemetry.
        
        Note: The output is an in silico statistical surrogate estimate.
        It is NOT an experimental measurement or clinical efficacy claim.
        """
        if not self.is_ready:
            raise RuntimeError(
                f"QSAR model artifact not loaded from {self.model_path}. Model must be trained first."
            )

        if not smiles or not isinstance(smiles, str):
            return {"success": False, "error": "Empty or non-string SMILES provided."}

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"success": False, "error": f"Invalid SMILES string: {smiles}"}

        canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
        fp, desc = self.featurizer.featurize_smiles(canonical_smiles)
        if fp is None or desc is None:
            return {"success": False, "error": "Featurization failed for molecule."}

        fp_batch = fp.reshape(1, -1)
        desc_batch = desc.reshape(1, -1)

        # Run model prediction
        if hasattr(self.model, "predict_with_dispersion"):
            # Check model signature
            if hasattr(self.model, "_prepare_features"):  # M2 Hybrid
                preds, disp = self.model.predict_with_dispersion(fp_batch, desc_batch)
            elif hasattr(self.model, "model") and hasattr(self.model, "scaler"):
                preds, disp = self.model.predict_with_dispersion(fp_batch, desc_batch)
            else:
                try:
                    preds, disp = self.model.predict_with_dispersion(fp_batch, desc_batch)
                except TypeError:
                    preds, disp = self.model.predict_with_dispersion(fp_batch)
        else:
            preds = self.model.predict(fp_batch)
            disp = np.zeros_like(preds)

        pred_pchembl = float(preds[0])
        pred_dispersion = float(disp[0])

        # Compute Applicability Domain (T_max)
        t_max = 0.0
        domain_status = "unassessed"
        if self._ref_morgan is not None:
            t_max_arr = compute_max_tanimoto(fp_batch, self._ref_morgan)
            t_max = float(t_max_arr[0])
            if t_max >= 0.60:
                domain_status = "high_similarity"
            elif t_max >= 0.45:
                domain_status = "moderate_similarity"
            elif t_max >= 0.35:
                domain_status = "low_similarity"
            else:
                domain_status = "out_of_domain_novel"

        # Convert pIC50 to predicted IC50 in nM: IC50 = 10^(9.0 - pIC50)
        pred_ic50_nm = float(10.0 ** (9.0 - pred_pchembl))

        return {
            "success": True,
            "smiles": canonical_smiles,
            "predicted_pchembl": round(pred_pchembl, 3),
            "predicted_ic50_nm": round(pred_ic50_nm, 2),
            "heuristic_dispersion_sigma": round(pred_dispersion, 3),
            "applicability_domain": {
                "max_tanimoto_train": round(t_max, 3),
                "domain_status": domain_status,
                "is_within_domain": t_max >= 0.35,
            },
            "epistemic_metadata": {
                "is_experimentally_validated": False,
                "target_target_chembl_id": "CHEMBL203",
                "target_name": "Epidermal growth factor receptor (EGFR)",
                "caveat": (
                    "In silico QSAR statistical surrogate estimate. Not an experimental "
                    "measurement or clinical efficacy validation."
                ),
            },
        }

    def predict_candidate(self, candidate: MoleculeCandidate) -> MoleculeCandidate:
        """Enrich a MoleculeCandidate with predicted affinity score and QSAR metadata."""
        pred = self.predict(candidate.smiles)
        if pred.get("success"):
            candidate.affinity_score = pred["predicted_pchembl"]
            if candidate.metadata is None:
                candidate.metadata = {}
            candidate.metadata["qsar_prediction"] = pred
        return candidate
