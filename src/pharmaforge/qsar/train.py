"""Full training and cross-validation execution engine for EGFR QSAR baselines.

Evaluates M0, M1, M2, and M3 under:
1. Random 5-Fold Cross-Validation.
2. Bemis-Murcko Scaffold 5-Fold Cross-Validation.
3. Document Connected-Component 5-Fold Cross-Validation.
4. Protocol-Transfer Evaluation on Cell-Free Membrane Vesicle assays.
5. Benchmark Holdout Evaluation.
"""

from __future__ import annotations
import os
import json
import csv
import time
import joblib
import numpy as np
from pathlib import Path
from scipy import stats

from pharmaforge.qsar.baselines import (
    M0DummyBaseline,
    M1RandomForestMorgan,
    M2RandomForestHybrid,
    M3HistGradientBoostingDesc,
)
from pharmaforge.qsar.splits import (
    create_random_kfold_splits,
    create_scaffold_kfold_splits,
    create_document_kfold_splits,
)
from pharmaforge.qsar.evaluator import (
    compute_regression_metrics,
    evaluate_uncertainty_coverage,
    compute_max_tanimoto,
)
from pharmaforge.qsar.features import MolecularFeaturizer, DESCRIPTOR_NAMES


def run_cv_experiment(
    model_name: str,
    split_name: str,
    splits: list[tuple[np.ndarray, np.ndarray]],
    X_m: np.ndarray,
    X_d: np.ndarray,
    y: np.ndarray,
    chembl_ids: np.ndarray,
    smiles: np.ndarray,
    scaffolds: np.ndarray,
) -> tuple[dict, list[dict]]:
    """Run 5-Fold CV for a specific baseline model and split protocol."""
    n_samples = len(y)
    oof_preds = np.zeros(n_samples, dtype=np.float32)
    oof_dispersion = np.zeros(n_samples, dtype=np.float32)
    oof_tmax = np.zeros(n_samples, dtype=np.float32)
    oof_folds = np.zeros(n_samples, dtype=np.int32)
    
    fold_metrics = []

    for fold_idx, (train_idx, test_idx) in enumerate(splits):
        # Instantiate fresh model
        if model_name == "M0_Dummy":
            model = M0DummyBaseline()
            model.fit(X_d[train_idx], y[train_idx])
            preds, disp = model.predict_with_dispersion(X_d[test_idx])
        elif model_name == "M1_RF_Morgan":
            model = M1RandomForestMorgan(n_estimators=100, random_state=42 + fold_idx)
            model.fit(X_m[train_idx], y[train_idx])
            preds, disp = model.predict_with_dispersion(X_m[test_idx])
        elif model_name == "M2_RF_Hybrid":
            model = M2RandomForestHybrid(n_estimators=100, random_state=42 + fold_idx)
            model.fit(X_m[train_idx], X_d[train_idx], y[train_idx])
            preds, disp = model.predict_with_dispersion(X_m[test_idx], X_d[test_idx])
        elif model_name == "M3_HGB_Desc":
            model = M3HistGradientBoostingDesc(max_iter=100, random_state=42 + fold_idx)
            model.fit(X_d[train_idx], y[train_idx])
            preds, disp = model.predict_with_dispersion(X_d[test_idx])
        else:
            raise ValueError(f"Unknown model name: {model_name}")

        oof_preds[test_idx] = preds
        oof_dispersion[test_idx] = disp
        oof_folds[test_idx] = fold_idx

        # Compute max Tanimoto to training set for test set
        tmax = compute_max_tanimoto(X_m[test_idx], X_m[train_idx])
        oof_tmax[test_idx] = tmax

        # Fold metrics
        f_m = compute_regression_metrics(y[test_idx], preds)
        fold_metrics.append(f_m)

    # Aggregate Out-of-Fold overall metrics
    overall_metrics = compute_regression_metrics(y, oof_preds)
    uncertainty_eval = evaluate_uncertainty_coverage(y, oof_preds, oof_dispersion)

    # Tanimoto similarity binning error analysis
    sim_bins = {
        "high_sim_ge_0.60": np.where(oof_tmax >= 0.60)[0],
        "mod_sim_0.45_0.60": np.where((oof_tmax >= 0.45) & (oof_tmax < 0.60))[0],
        "low_sim_0.35_0.45": np.where((oof_tmax >= 0.35) & (oof_tmax < 0.45))[0],
        "novel_lt_0.35": np.where(oof_tmax < 0.35)[0],
    }
    binned_metrics = {}
    for bin_name, bin_idx in sim_bins.items():
        if len(bin_idx) > 0:
            bm = compute_regression_metrics(y[bin_idx], oof_preds[bin_idx])
            bm["count"] = len(bin_idx)
            binned_metrics[bin_name] = bm

    summary = {
        "model_name": model_name,
        "split_protocol": split_name,
        "overall_metrics": overall_metrics,
        "fold_metrics": fold_metrics,
        "uncertainty_evaluation": uncertainty_eval,
        "tanimoto_binned_metrics": binned_metrics,
    }

    # Prepare per-molecule prediction rows
    pred_rows = []
    for i in range(n_samples):
        pred_rows.append({
            "molecule_chembl_id": str(chembl_ids[i]),
            "standardized_smiles": str(smiles[i]),
            "bemis_murcko_scaffold": str(scaffolds[i]),
            "fold": int(oof_folds[i]),
            "y_true": round(float(y[i]), 4),
            "y_pred": round(float(oof_preds[i]), 4),
            "abs_error": round(float(abs(y[i] - oof_preds[i])), 4),
            "dispersion_sigma": round(float(oof_dispersion[i]), 4),
            "max_tanimoto_train": round(float(oof_tmax[i]), 4),
        })

    return summary, pred_rows


def load_cached_cv_summary(pred_file: Path, model_name: str, split_name: str) -> dict:
    """Reconstruct experiment summary from saved prediction CSV."""
    y_true_list = []
    y_pred_list = []
    disp_list = []
    tmax_list = []
    folds_dict = {}

    with open(pred_file, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            yt = float(row["y_true"])
            yp = float(row["y_pred"])
            disp = float(row.get("dispersion_sigma", 0.0))
            tmax = float(row.get("max_tanimoto_train", 0.0))
            f_idx = int(row["fold"])
            y_true_list.append(yt)
            y_pred_list.append(yp)
            disp_list.append(disp)
            tmax_list.append(tmax)
            if f_idx not in folds_dict:
                folds_dict[f_idx] = ([], [])
            folds_dict[f_idx][0].append(yt)
            folds_dict[f_idx][1].append(yp)

    y_true = np.array(y_true_list, dtype=np.float32)
    y_pred = np.array(y_pred_list, dtype=np.float32)
    dispersion = np.array(disp_list, dtype=np.float32)
    oof_tmax = np.array(tmax_list, dtype=np.float32)

    overall_metrics = compute_regression_metrics(y_true, y_pred)
    uncertainty_eval = evaluate_uncertainty_coverage(y_true, y_pred, dispersion)

    fold_metrics = []
    for f_idx in sorted(folds_dict.keys()):
        yt_f, yp_f = folds_dict[f_idx]
        fold_metrics.append(compute_regression_metrics(np.array(yt_f), np.array(yp_f)))

    sim_bins = {
        "high_sim_ge_0.60": np.where(oof_tmax >= 0.60)[0],
        "mod_sim_0.45_0.60": np.where((oof_tmax >= 0.45) & (oof_tmax < 0.60))[0],
        "low_sim_0.35_0.45": np.where((oof_tmax >= 0.35) & (oof_tmax < 0.45))[0],
        "novel_lt_0.35": np.where(oof_tmax < 0.35)[0],
    }
    binned_metrics = {}
    for bin_name, bin_idx in sim_bins.items():
        if len(bin_idx) > 0:
            bm = compute_regression_metrics(y_true[bin_idx], y_pred[bin_idx])
            bm["count"] = len(bin_idx)
            binned_metrics[bin_name] = bm

    return {
        "model_name": model_name,
        "split_protocol": split_name,
        "overall_metrics": overall_metrics,
        "fold_metrics": fold_metrics,
        "uncertainty_evaluation": uncertainty_eval,
        "tanimoto_binned_metrics": binned_metrics,
    }


def train_and_evaluate_all():
    artifacts_dir = Path("artifacts/qsar")
    artifacts_dir.mkdir(parents=True, exist_ok=True)

    # Load precomputed data
    npz_path = "data/qsar_precomputed_features.npz"
    if not os.path.exists(npz_path):
        raise FileNotFoundError(f"Cached features {npz_path} not found. Run precompute_features.py first.")

    data = np.load(npz_path)
    X_m_core = data["X_m_core"]
    X_d_core = data["X_d_core"]
    y_core = data["y_core"]
    chembl_ids_core = data["chembl_ids_core"]
    smiles_core = data["smiles_core"]
    scaffolds_core = data["scaffolds_core"]
    docs_core = data["docs_core"]

    X_m_mem = data["X_m_mem"]
    X_d_mem = data["X_d_mem"]
    y_mem = data["y_mem"]
    chembl_ids_mem = data["chembl_ids_mem"]
    smiles_mem = data["smiles_mem"]

    X_m_bench = data["X_m_bench"]
    X_d_bench = data["X_d_bench"]
    y_bench = data["y_bench"]
    chembl_ids_bench = data["chembl_ids_bench"]
    smiles_bench = data["smiles_bench"]

    n_core = len(y_core)
    print(f"Loaded core dataset: {n_core} compounds.")

    # 1. Generate splits
    random_splits = create_random_kfold_splits(n_core, n_splits=5, seed=42)
    scaffold_splits = create_scaffold_kfold_splits(scaffolds_core, n_splits=5, seed=42)
    
    doc_list = [d.split("|") for d in docs_core]
    doc_splits = create_document_kfold_splits(doc_list, n_splits=5, seed=42)

    all_experiments = {}
    models = ["M0_Dummy", "M1_RF_Morgan", "M2_RF_Hybrid", "M3_HGB_Desc"]

    # 2. Run Random 5-Fold CV and Scaffold 5-Fold CV
    for model_name in models:
        for split_name, splits in [("random_5fold", random_splits), ("scaffold_5fold", scaffold_splits)]:
            exp_key = f"{model_name}__{split_name}"
            pred_file = artifacts_dir / f"predictions_{exp_key}.csv"
            if pred_file.exists():
                print(f"Loading cached {exp_key}...")
                summary = load_cached_cv_summary(pred_file, model_name, split_name)
            else:
                print(f"Running {model_name} under {split_name}...")
                t0 = time.perf_counter()
                summary, pred_rows = run_cv_experiment(
                    model_name=model_name,
                    split_name=split_name,
                    splits=splits,
                    X_m=X_m_core,
                    X_d=X_d_core,
                    y=y_core,
                    chembl_ids=chembl_ids_core,
                    smiles=smiles_core,
                    scaffolds=scaffolds_core,
                )
                elapsed = time.perf_counter() - t0
                summary["runtime_seconds"] = round(elapsed, 2)
                with open(pred_file, "w", encoding="utf-8", newline="") as f:
                    writer = csv.DictWriter(f, fieldnames=list(pred_rows[0].keys()))
                    writer.writeheader()
                    writer.writerows(pred_rows)
            all_experiments[exp_key] = summary
            print(f"  {exp_key} -> MAE: {summary['overall_metrics']['mae']}, RMSE: {summary['overall_metrics']['rmse']}, R2: {summary['overall_metrics']['r2']}, rho: {summary['overall_metrics']['spearman_rho']}")

    # 3. Run Document-Disjoint CV for the best model family (M1 and M2)
    for model_name in ["M1_RF_Morgan", "M2_RF_Hybrid"]:
        exp_key = f"{model_name}__document_5fold"
        pred_file = artifacts_dir / f"predictions_{exp_key}.csv"
        if pred_file.exists():
            print(f"Loading cached {exp_key}...")
            summary = load_cached_cv_summary(pred_file, model_name, "document_5fold")
        else:
            print(f"Running {model_name} under document_5fold...")
            t0 = time.perf_counter()
            summary, pred_rows = run_cv_experiment(
                model_name=model_name,
                split_name="document_5fold",
                splits=doc_splits,
                X_m=X_m_core,
                X_d=X_d_core,
                y=y_core,
                chembl_ids=chembl_ids_core,
                smiles=smiles_core,
                scaffolds=scaffolds_core,
            )
            elapsed = time.perf_counter() - t0
            summary["runtime_seconds"] = round(elapsed, 2)
            with open(pred_file, "w", encoding="utf-8", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=list(pred_rows[0].keys()))
                writer.writeheader()
                writer.writerows(pred_rows)
        all_experiments[exp_key] = summary
        print(f"  {exp_key} -> MAE: {summary['overall_metrics']['mae']}, RMSE: {summary['overall_metrics']['rmse']}, R2: {summary['overall_metrics']['r2']}, rho: {summary['overall_metrics']['spearman_rho']}")

    # 4. Train Full Core Models and Evaluate Protocol-Transfer (Membrane) & Benchmarks
    print("\nTraining full production models on all 5,226 core compounds...")
    full_models = {}
    
    # Train M0
    m0_full = M0DummyBaseline().fit(X_d_core, y_core)
    full_models["M0_Dummy"] = m0_full
    
    # Train M1
    m1_full = M1RandomForestMorgan(n_estimators=100, random_state=42).fit(X_m_core, y_core)
    full_models["M1_RF_Morgan"] = m1_full
    
    # Train M2
    m2_full = M2RandomForestHybrid(n_estimators=100, random_state=42).fit(X_m_core, X_d_core, y_core)
    full_models["M2_RF_Hybrid"] = m2_full
    
    # Train M3
    m3_full = M3HistGradientBoostingDesc(max_iter=100, random_state=42).fit(X_d_core, y_core)
    full_models["M3_HGB_Desc"] = m3_full

    # Protocol transfer evaluation on membrane vesicle assays
    transfer_results = {}
    print("\nEvaluating Protocol-Transfer on Cell-Free Membrane Vesicle assays (N=168)...")
    for name, model in full_models.items():
        if name == "M0_Dummy":
            p, disp = model.predict_with_dispersion(X_d_mem)
        elif name == "M1_RF_Morgan":
            p, disp = model.predict_with_dispersion(X_m_mem)
        elif name == "M2_RF_Hybrid":
            p, disp = model.predict_with_dispersion(X_m_mem, X_d_mem)
        elif name == "M3_HGB_Desc":
            p, disp = model.predict_with_dispersion(X_d_mem)
            
        m = compute_regression_metrics(y_mem, p)
        # Residual correlation with LogP (Descriptor index 1 is MolLogP)
        residuals = y_mem - p
        logp_vals = X_d_mem[:, 1]
        logp_corr, _ = stats.spearmanr(logp_vals, residuals)
        m["residual_logp_spearman"] = round(float(logp_corr), 4)
        m["mean_residual_bias"] = round(float(np.mean(residuals)), 4)
        transfer_results[name] = m
        print(f"  {name} on Membrane -> MAE: {m['mae']}, RMSE: {m['rmse']}, R2: {m['r2']}, rho: {m['spearman_rho']}, Bias: {m['mean_residual_bias']}, LogP Corr: {m['residual_logp_spearman']}")

    # Benchmark leads holdout evaluation
    benchmark_eval = {}
    print("\nEvaluating Benchmark Leads Holdout (CHEMBL289959 & CHEMBL290096)...")
    for name, model in full_models.items():
        if name == "M0_Dummy":
            p, disp = model.predict_with_dispersion(X_d_bench)
        elif name == "M1_RF_Morgan":
            p, disp = model.predict_with_dispersion(X_m_bench)
        elif name == "M2_RF_Hybrid":
            p, disp = model.predict_with_dispersion(X_m_bench, X_d_bench)
        elif name == "M3_HGB_Desc":
            p, disp = model.predict_with_dispersion(X_d_bench)
            
        tmax = compute_max_tanimoto(X_m_bench, X_m_core)
        bench_entries = []
        for b_idx in range(len(y_bench)):
            bench_entries.append({
                "molecule_chembl_id": str(chembl_ids_bench[b_idx]),
                "y_true": round(float(y_bench[b_idx]), 4),
                "y_pred": round(float(p[b_idx]), 4),
                "error": round(float(abs(y_bench[b_idx] - p[b_idx])), 4),
                "dispersion_sigma": round(float(disp[b_idx]), 4),
                "max_tanimoto_train": round(float(tmax[b_idx]), 4),
            })
        benchmark_eval[name] = bench_entries
        print(f"  {name} -> {bench_entries[0]['molecule_chembl_id']}: pred={bench_entries[0]['y_pred']} (true={bench_entries[0]['y_true']}), {bench_entries[1]['molecule_chembl_id']}: pred={bench_entries[1]['y_pred']} (true={bench_entries[1]['y_true']})")

    # 5. Save the best model artifact (M2_RF_Hybrid)
    best_model_name = "M2_RF_Hybrid"
    best_model = full_models[best_model_name]
    best_model_path = artifacts_dir / "egfr_qsar_best_model.joblib"
    joblib.dump(best_model, best_model_path)
    print(f"\nSaved best model ({best_model_name}) artifact to {best_model_path}")

    # Save feature configuration
    featurizer = MolecularFeaturizer()
    feature_config = featurizer.get_feature_config()
    with open(artifacts_dir / "feature_config.json", "w", encoding="utf-8") as f:
        json.dump(feature_config, f, indent=2)

    # Save training metadata
    training_metadata = {
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "target": "EGFR (CHEMBL203)",
        "dataset_name": "Core Homogeneous Cell-Free Biochemical Subset",
        "training_samples_count": n_core,
        "features": feature_config,
        "models_evaluated": models,
        "best_model_selected": best_model_name,
        "model_artifact_path": str(best_model_path),
        "membrane_transfer_samples_count": len(y_mem),
        "benchmark_holdout_samples_count": len(y_bench),
    }
    with open(artifacts_dir / "training_metadata.json", "w", encoding="utf-8") as f:
        json.dump(training_metadata, f, indent=2)

    # Save complete metrics summary
    metrics_summary = {
        "cross_validation_experiments": all_experiments,
        "protocol_transfer_membrane": transfer_results,
        "benchmark_holdout_evaluation": benchmark_eval,
    }
    with open(artifacts_dir / "metrics_summary.json", "w", encoding="utf-8") as f:
        json.dump(metrics_summary, f, indent=2)

    print(f"Saved complete metrics summary to {artifacts_dir / 'metrics_summary.json'}")
    return metrics_summary

if __name__ == "__main__":
    train_and_evaluate_all()
