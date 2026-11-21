"""Molecular feature extraction for QSAR modeling.

Computes:
1. Morgan circular fingerprints (ECFP4 equivalent, radius=2, 2048 bits).
2. Interpretable physicochemical descriptors (10 physical properties).
3. Hybrid representation with standard scaling for continuous features.
"""

from __future__ import annotations
import numpy as np
from typing import Sequence
from rdkit import Chem
from rdkit.Chem import Descriptors, QED, rdFingerprintGenerator

DESCRIPTOR_NAMES = [
    "MolWt",
    "MolLogP",
    "TPSA",
    "NumHDonors",
    "NumHAcceptors",
    "NumRotatableBonds",
    "NumAromaticRings",
    "HeavyAtomCount",
    "FractionCSP3",
    "QED",
]

class MolecularFeaturizer:
    """Deterministic featurizer extracting Morgan fingerprints and physicochemical descriptors."""

    def __init__(self, radius: int = 2, n_bits: int = 2048):
        self.radius = radius
        self.n_bits = n_bits
        self._fp_gen = rdFingerprintGenerator.GetMorganGenerator(radius=self.radius, fpSize=self.n_bits)

    def compute_morgan_fp(self, mol: Chem.Mol) -> np.ndarray:
        """Compute Morgan circular fingerprint (ECFP4 equivalent)."""
        fp = self._fp_gen.GetFingerprint(mol)
        arr = np.zeros((self.n_bits,), dtype=np.float32)
        # Convert ExplicitBitVect to numpy array
        for bit_idx in fp.GetOnBits():
            arr[bit_idx] = 1.0
        return arr

    def compute_descriptors(self, mol: Chem.Mol) -> np.ndarray:
        """Compute 10 standardized physicochemical descriptors."""
        values = [
            float(Descriptors.MolWt(mol)),
            float(Descriptors.MolLogP(mol)),
            float(Descriptors.TPSA(mol)),
            float(Descriptors.NumHDonors(mol)),
            float(Descriptors.NumHAcceptors(mol)),
            float(Descriptors.NumRotatableBonds(mol)),
            float(Descriptors.NumAromaticRings(mol)),
            float(Descriptors.HeavyAtomCount(mol)),
            float(Descriptors.FractionCSP3(mol)),
            float(QED.qed(mol)),
        ]
        return np.array(values, dtype=np.float32)

    def featurize_smiles(self, smiles: str) -> tuple[np.ndarray | None, np.ndarray | None]:
        """Featurize a single SMILES string into (morgan_fp, descriptors)."""
        if not smiles or not isinstance(smiles, str):
            return None, None
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None, None
        fp = self.compute_morgan_fp(mol)
        desc = self.compute_descriptors(mol)
        return fp, desc

    def featurize_smiles_list(
        self, smiles_list: Sequence[str]
    ) -> tuple[np.ndarray, np.ndarray, list[int]]:
        """Batch featurize a list of SMILES strings.
        
        Returns:
            X_morgan: (N_valid, n_bits)
            X_desc: (N_valid, 10)
            valid_indices: list of original indices that successfully featurized.
        """
        valid_fps = []
        valid_descs = []
        valid_indices = []

        for idx, smi in enumerate(smiles_list):
            fp, desc = self.featurize_smiles(smi)
            if fp is not None and desc is not None:
                valid_fps.append(fp)
                valid_descs.append(desc)
                valid_indices.append(idx)

        X_morgan = np.vstack(valid_fps) if valid_fps else np.empty((0, self.n_bits), dtype=np.float32)
        X_desc = np.vstack(valid_descs) if valid_descs else np.empty((0, len(DESCRIPTOR_NAMES)), dtype=np.float32)
        return X_morgan, X_desc, valid_indices

    def get_feature_config(self) -> dict:
        """Return metadata configuration dictionary for model persistence."""
        return {
            "morgan_radius": self.radius,
            "morgan_n_bits": self.n_bits,
            "descriptor_names": DESCRIPTOR_NAMES,
            "descriptor_dim": len(DESCRIPTOR_NAMES),
            "hybrid_dim": self.n_bits + len(DESCRIPTOR_NAMES),
        }
