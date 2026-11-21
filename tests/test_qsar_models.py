"""Targeted test suite for QSAR model baselines, feature generation, splitters, and predictor integration."""
import pytest
import numpy as np
from pathlib import Path
from rdkit import Chem

from pharmaforge.qsar.features import MolecularFeaturizer, DESCRIPTOR_NAMES
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
from pharmaforge.qsar.evaluator import compute_regression_metrics, evaluate_uncertainty_coverage
from pharmaforge.qsar.predictor import EGFRQSARPredictor
from pharmaforge.core.types import MoleculeCandidate, LeadOptimizationTask, OptimizationObjective, ObjectiveDirection
from pharmaforge.agents.orchestrator import PharmaForgeOrchestrator


def test_featurizer_dimensions():
    featurizer = MolecularFeaturizer(radius=2, n_bits=2048)
    smiles = "Nc1ncnc2c1c1cccc(Cl)c1n2-c1cccc(Br)c1"
    fp, desc = featurizer.featurize_smiles(smiles)
    
    assert fp is not None
    assert desc is not None
    assert fp.shape == (2048,)
    assert desc.shape == (10,)
    assert len(DESCRIPTOR_NAMES) == 10


def test_baselines_fit_predict():
    X_m = np.random.binomial(1, 0.05, size=(20, 2048)).astype(np.float32)
    X_d = np.random.normal(0, 1, size=(20, 10)).astype(np.float32)
    y = np.random.uniform(5.0, 9.0, size=(20,)).astype(np.float32)

    # M0 Dummy
    m0 = M0DummyBaseline().fit(X_d[:15], y[:15])
    p0, d0 = m0.predict_with_dispersion(X_d[15:])
    assert len(p0) == 5
    assert np.allclose(d0, 0.0)

    # M1 RF Morgan
    m1 = M1RandomForestMorgan(n_estimators=10, random_state=42).fit(X_m[:15], y[:15])
    p1, d1 = m1.predict_with_dispersion(X_m[15:])
    assert len(p1) == 5
    assert len(d1) == 5

    # M2 RF Hybrid
    m2 = M2RandomForestHybrid(n_estimators=10, random_state=42).fit(X_m[:15], X_d[:15], y[:15])
    p2, d2 = m2.predict_with_dispersion(X_m[15:], X_d[15:])
    assert len(p2) == 5
    assert len(d2) == 5

    # M3 HGB Desc
    m3 = M3HistGradientBoostingDesc(max_iter=10, random_state=42).fit(X_d[:15], y[:15])
    p3, d3 = m3.predict_with_dispersion(X_d[15:])
    assert len(p3) == 5


def test_scaffold_disjoint_split():
    scaffolds = ["scaf_A", "scaf_A", "scaf_B", "scaf_B", "scaf_C", "scaf_D", "scaf_E"]
    splits = create_scaffold_kfold_splits(scaffolds, n_splits=3, seed=42)
    
    for train_idx, test_idx in splits:
        train_scafs = set(scaffolds[i] for i in train_idx)
        test_scafs = set(scaffolds[i] for i in test_idx)
        assert train_scafs.isdisjoint(test_scafs), "Scaffolds must be strictly disjoint between train and test!"


def test_document_disjoint_split():
    docs = [["doc_1"], ["doc_1", "doc_2"], ["doc_2"], ["doc_3"], ["doc_4"]]
    splits = create_document_kfold_splits(docs, n_splits=2, seed=42)
    
    for train_idx, test_idx in splits:
        train_docs = set(d for i in train_idx for d in docs[i])
        test_docs = set(d for i in test_idx for d in docs[i])
        assert train_docs.isdisjoint(test_docs), "Documents must be strictly disjoint between train and test!"


def test_qsar_predictor_epistemic_safety():
    predictor = EGFRQSARPredictor()
    if not predictor.is_ready:
        pytest.skip("Model artifact not yet present.")
        
    res = predictor.predict("Nc1ncnc2c1c1cccc(Cl)c1n2-c1cccc(Br)c1")
    assert res["success"] is True
    assert "predicted_pchembl" in res
    assert "heuristic_dispersion_sigma" in res
    assert "applicability_domain" in res
    
    # Must explicitly state non-experimental validation
    epistemic = res["epistemic_metadata"]
    assert epistemic["is_experimentally_validated"] is False
    assert "caveat" in epistemic
    assert "Not an experimental measurement" in epistemic["caveat"]


def test_orchestrator_qsar_integration():
    orch = PharmaForgeOrchestrator()
    assert orch.qsar_predictor is not None
    assert orch.qsar_predictor.is_ready is True
    
    task = LeadOptimizationTask(
        task_id="test_egfr_integration",
        target_name="EGFR",
        initial_smiles="Nc1ncnc2c1c1cccc(Cl)c1n2-c1cccc(Br)c1",
        objectives=[
            OptimizationObjective(name="affinity", direction=ObjectiveDirection.MAXIMIZE, threshold=6.0),
        ],
        max_iterations=1,
    )
    res = orch.run_optimization(task)
    assert res["success"] is True
    assert len(res["pareto_front"]) >= 1
    top_cand = res["pareto_front"][0]
    assert top_cand["affinity_score"] is not None
    assert "qsar_prediction" in top_cand["metadata"]
    assert top_cand["metadata"]["qsar_prediction"]["epistemic_metadata"]["is_experimentally_validated"] is False
