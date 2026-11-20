"""Unit tests for the ExecutionTraceSelfRepairEngine (ET-SRE)."""
import pytest
from pharmaforge.agents.repair import ExecutionTraceSelfRepairEngine
from pharmaforge.core.exceptions import SkillRepairMaxAttemptsExceededError
from pharmaforge.sandbox.executor import IsolatedCodeExecutor


@pytest.fixture
def repair_engine():
    executor = IsolatedCodeExecutor(timeout_seconds=5.0)
    return ExecutionTraceSelfRepairEngine(executor=executor, max_repair_attempts=3)


def test_repair_missing_module_import(repair_engine):
    """Test repairing code that failed due to a missing standard module import."""
    buggy_code = """
def calculate_hypotenuse(a: float, b: float) -> float:
    # Notice: 'math' is not imported
    return math.sqrt(a**2 + b**2)
"""
    unit_tests = """
def run_tests():
    assert calculate_hypotenuse(3.0, 4.0) == 5.0
    return True
"""
    simulated_traceback = "NameError: name 'math' is not defined"

    repaired_code, changelog = repair_engine.repair_code(
        skill_name="calculate_hypotenuse",
        failing_code=buggy_code,
        unit_tests=unit_tests,
        traceback_str=simulated_traceback,
    )
    assert "import math" in repaired_code
    assert len(changelog) >= 1
    assert "Success" in changelog[-1]


def test_repair_unfixable_code_exceeds_budget(repair_engine):
    """Test that irrecoverably broken code exhausts repair budget and raises error."""
    broken_code = """
def completely_broken():
    raise RuntimeError("Permanent failure")
"""
    unit_tests = """
def run_tests():
    assert completely_broken() == 42
    return True
"""
    with pytest.raises(SkillRepairMaxAttemptsExceededError):
        repair_engine.repair_code(
            skill_name="completely_broken",
            failing_code=broken_code,
            unit_tests=unit_tests,
            traceback_str="RuntimeError: Permanent failure",
        )
