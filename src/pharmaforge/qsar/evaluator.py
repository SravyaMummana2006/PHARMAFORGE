"""QSAR Evaluation engine for metrics, uncertainty curves, applicability domain, and protocol transfer.

Computes:
1. MAE, RMSE, R2, and Spearman correlation.
2. Error-versus-retained-coverage analysis across confidence deciles.
3. Nearest-neighbor Morgan fingerprint Tanimoto applicability domain (T_max).
4. Protocol-transfer evaluation on membrane vesicle assays.
5. Benchmark lead holdout evaluation.
"""

from __future__ import annotations
import numpy as np
from scipy import stats
from sklearn.metrics import mean_absolute_error, root_mean_squared_error, r2_score

def compute_regression_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    """Compute standard regression metrics (MAE, RMSE, R2, Spearman rho)."""
    mae = float(mean_absolute_error(y_true, y_pred))
    rmse = float(root_mean_squared_error(y_true, y_pred))
    
    # R2 can be negative or positive
    r2 = float(r2_score(y_true, y_pred))
    
    # Spearman rank correlation
    if len(y_true) > 1 and np.std(y_pred) > 1e-9 and np.std(y_true) > 1e-9:
        rho, _ = stats.spearmanr(y_true, y_pred)
        spearman_rho = float(rho)
    else:
        spearman_rho = 0.0

    return {
        "mae": round(mae, 4),
        "rmse": round(rmse, 4),
        "r2": round(r2, 4),
        "spearman_rho": round(spearman_rho, 4),
    }

def evaluate_uncertainty_coverage(
    y_true: np.ndarray, y_pred: np.ndarray, dispersion: np.ndarray
) -> dict[str, any]:
    """Evaluate error-versus-retained-coverage curves across coverage percentiles."""
    abs_errors = np.abs(y_true - y_pred)
    
    # Rank correlation between absolute error and dispersion
    if len(abs_errors) > 1 and np.std(dispersion) > 1e-9:
        rho_err, _ = stats.spearmanr(abs_errors, dispersion)
        spearman_error_rho = float(rho_err)
    else:
        spearman_error_rho = 0.0

    # Sort compounds by dispersion ascending (most confident first)
    sorted_order = np.argsort(dispersion)
    n_total = len(y_true)
    
    coverage_tiers = [1.0, 0.9, 0.8, 0.7, 0.5]
    curve_results = {}
    
    for cov in coverage_tiers:
        k = max(1, int(round(n_total * cov)))
        subset_idx = sorted_order[:k]
        sub_mae = float(mean_absolute_error(y_true[subset_idx], y_pred[subset_idx]))
        sub_rmse = float(root_mean_squared_error(y_true[subset_idx], y_pred[subset_idx]))
        curve_results[f"coverage_{int(cov*100)}%"] = {
            "n_retained": k,
            "mae": round(sub_mae, 4),
            "rmse": round(sub_rmse, 4),
            "max_dispersion": round(float(dispersion[sorted_order[k-1]]), 4),
        }

    return {
        "spearman_error_dispersion_rho": round(spearman_error_rho, 4),
        "coverage_curve": curve_results,
    }

def compute_tanimoto_similarity_matrix(X_query: np.ndarray, X_ref: np.ndarray) -> np.ndarray:
    """Compute pairwise Tanimoto similarity between binary Morgan fingerprint matrices.
    
    Formula: T(a, b) = (a . b) / (|a| + |b| - a . b)
    """
    intersection = np.dot(X_query, X_ref.T)  # (N_q, N_ref)
    sum_q = np.sum(X_query, axis=1, keepdims=True)  # (N_q, 1)
    sum_ref = np.sum(X_ref, axis=1, keepdims=True).T  # (1, N_ref)
    union = sum_q + sum_ref - intersection
    
    with np.errstate(divide="ignore", invalid="ignore"):
        tanimoto = np.where(union > 0, intersection / union, 0.0)
    return tanimoto

def compute_max_tanimoto(X_query: np.ndarray, X_ref: np.ndarray) -> np.ndarray:
    """Compute maximum Tanimoto similarity to any reference molecule."""
    # Process in batches to conserve memory
    batch_size = 250
    n_query = X_query.shape[0]
    max_t = np.zeros(n_query, dtype=np.float32)
    
    for i in range(0, n_query, batch_size):
        q_batch = X_query[i : i + batch_size]
        sims = compute_tanimoto_similarity_matrix(q_batch, X_ref)
        max_t[i : i + batch_size] = np.max(sims, axis=1)
        
    return max_t
