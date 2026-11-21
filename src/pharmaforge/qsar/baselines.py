"""Baseline regression model wrappers for EGFR QSAR affinity prediction.

Implements the four approved model baselines:
- M0: DummyRegressor (mean prediction baseline)
- M1: RandomForestRegressor on Morgan fingerprints (2048 bits)
- M2: RandomForestRegressor on Hybrid representation (Morgan 2048 + scaled descriptors)
- M3: HistGradientBoostingRegressor on 10 physicochemical descriptors
"""

from __future__ import annotations
import numpy as np
from typing import Any
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import RandomForestRegressor, HistGradientBoostingRegressor
from sklearn.preprocessing import StandardScaler

class M0DummyBaseline:
    """M0: Baseline statistical floor predicting the training partition mean."""

    def __init__(self):
        self.model = DummyRegressor(strategy="mean")
        self.mean_value_: float = 0.0

    def fit(self, X: np.ndarray, y: np.ndarray) -> "M0DummyBaseline":
        self.model.fit(X, y)
        self.mean_value_ = float(self.model.constant_[0][0])
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        return self.model.predict(X)

    def predict_with_dispersion(self, X: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        preds = self.predict(X)
        dispersion = np.zeros_like(preds, dtype=np.float32)
        return preds, dispersion


class M1RandomForestMorgan:
    """M1: Random Forest on 2048-bit Morgan circular fingerprints."""

    def __init__(self, n_estimators: int = 100, random_state: int = 42, n_jobs: int = -1):
        self.model = RandomForestRegressor(
            n_estimators=n_estimators,
            max_features="sqrt",
            random_state=random_state,
            n_jobs=n_jobs,
        )

    def fit(self, X_morgan: np.ndarray, y: np.ndarray) -> "M1RandomForestMorgan":
        self.model.fit(X_morgan, y)
        return self

    def predict(self, X_morgan: np.ndarray) -> np.ndarray:
        return self.model.predict(X_morgan)

    def predict_with_dispersion(self, X_morgan: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Compute ensemble mean prediction and standard deviation across tree estimators."""
        all_tree_preds = np.array([tree.predict(X_morgan) for tree in self.model.estimators_])
        mean_preds = np.mean(all_tree_preds, axis=0)
        dispersion = np.std(all_tree_preds, axis=0)
        return mean_preds, dispersion


class M2RandomForestHybrid:
    """M2: Random Forest on Hybrid features (Morgan 2048 + 10 Scaled Descriptors)."""

    def __init__(self, n_estimators: int = 100, random_state: int = 42, n_jobs: int = -1):
        self.scaler = StandardScaler()
        self.model = RandomForestRegressor(
            n_estimators=n_estimators,
            max_features="sqrt",
            random_state=random_state,
            n_jobs=n_jobs,
        )

    def _prepare_features(self, X_morgan: np.ndarray, X_desc: np.ndarray, fit_scaler: bool = False) -> np.ndarray:
        if fit_scaler:
            desc_scaled = self.scaler.fit_transform(X_desc)
        else:
            desc_scaled = self.scaler.transform(X_desc)
        return np.hstack([X_morgan, desc_scaled])

    def fit(self, X_morgan: np.ndarray, X_desc: np.ndarray, y: np.ndarray) -> "M2RandomForestHybrid":
        X_hybrid = self._prepare_features(X_morgan, X_desc, fit_scaler=True)
        self.model.fit(X_hybrid, y)
        return self

    def predict(self, X_morgan: np.ndarray, X_desc: np.ndarray) -> np.ndarray:
        X_hybrid = self._prepare_features(X_morgan, X_desc, fit_scaler=False)
        return self.model.predict(X_hybrid)

    def predict_with_dispersion(self, X_morgan: np.ndarray, X_desc: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        X_hybrid = self._prepare_features(X_morgan, X_desc, fit_scaler=False)
        all_tree_preds = np.array([tree.predict(X_hybrid) for tree in self.model.estimators_])
        mean_preds = np.mean(all_tree_preds, axis=0)
        dispersion = np.std(all_tree_preds, axis=0)
        return mean_preds, dispersion


class M3HistGradientBoostingDesc:
    """M3: HistGradientBoosting on 10 physicochemical descriptors."""

    def __init__(self, max_iter: int = 100, random_state: int = 42):
        self.model = HistGradientBoostingRegressor(
            max_iter=max_iter,
            random_state=random_state,
        )

    def fit(self, X_desc: np.ndarray, y: np.ndarray) -> "M3HistGradientBoostingDesc":
        self.model.fit(X_desc, y)
        return self

    def predict(self, X_desc: np.ndarray) -> np.ndarray:
        return self.model.predict(X_desc)

    def predict_with_dispersion(self, X_desc: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        preds = self.predict(X_desc)
        # Gradient boosting does not have independent tree bags; dispersion is 0.0 heuristic
        dispersion = np.zeros_like(preds, dtype=np.float32)
        return preds, dispersion
