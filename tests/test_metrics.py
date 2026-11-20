"""Unit tests for BenchmarkEvaluator metrics."""
import pytest
from pharmaforge.benchmarks.metrics import BenchmarkEvaluator
from pharmaforge.core.types import MoleculeCandidate, OptimizationObjective, ObjectiveDirection


def test_chemical_validity_calculation():
    """Test validity percentage calculation."""
    assert BenchmarkEvaluator.calculate_chemical_validity_rate(100, 95) == 95.0
    assert BenchmarkEvaluator.calculate_chemical_validity_rate(0, 0) == 0.0


def test_skill_reuse_calculation():
    """Test skill reuse ratio calculation."""
    assert BenchmarkEvaluator.calculate_skill_reuse_ratio(10, 6) == 60.0


def test_pareto_hypervolume_calculation():
    """Test hypervolume metric computation over candidates."""
    candidates = [
        MoleculeCandidate(smiles="CC", qed=0.8, affinity_score=8.5, logp=2.0),
        MoleculeCandidate(smiles="CCC", qed=0.6, affinity_score=7.0, logp=2.5),
    ]
    objectives = [
        OptimizationObjective(name="qed", direction=ObjectiveDirection.MAXIMIZE),
        OptimizationObjective(name="affinity_score", direction=ObjectiveDirection.MAXIMIZE),
    ]
    hv = BenchmarkEvaluator.calculate_pareto_hypervolume(candidates, objectives)
    assert 0.0 < hv <= 1.0
