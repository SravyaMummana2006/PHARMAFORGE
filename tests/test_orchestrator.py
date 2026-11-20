"""Integration tests for the complete closed-loop PharmaForgeOrchestrator."""
import pytest
from pharmaforge.agents.orchestrator import PharmaForgeOrchestrator
from pharmaforge.core.types import LeadOptimizationTask, OptimizationObjective, ObjectiveDirection


@pytest.fixture
def orchestrator():
    return PharmaForgeOrchestrator()


def test_orchestrator_closed_loop_optimization(orchestrator):
    """Test running an end-to-end multi-objective optimization task."""
    task = LeadOptimizationTask(
        task_id="TEST-TASK-01",
        target_name="TestTarget",
        initial_smiles="c1ccccc1O",  # Phenol
        objectives=[
            OptimizationObjective(name="qed", direction=ObjectiveDirection.MAXIMIZE, threshold=0.5),
            OptimizationObjective(name="logp", direction=ObjectiveDirection.MINIMIZE, threshold=3.0),
        ],
        max_iterations=3,
    )

    result = orchestrator.run_optimization(task)
    assert result["success"] is True
    assert result["total_candidates_generated"] >= 1
    assert result["pareto_front_size"] >= 1
    assert "trajectory_summary" in result
    assert result["wall_clock_seconds"] >= 0.0


def test_orchestrator_invalid_initial_lead_handling(orchestrator):
    """Test that invalid initial lead is handled gracefully without crashing."""
    task = LeadOptimizationTask(
        task_id="TEST-FAIL-01",
        target_name="TestTarget",
        initial_smiles="INVALID_SMILES_STRING",
        objectives=[],
        max_iterations=2,
    )
    result = orchestrator.run_optimization(task)
    assert result["success"] is False
    assert "failed verification" in result["error"]
