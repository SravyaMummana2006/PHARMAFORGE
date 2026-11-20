"""Unit tests for the Verified Procedural Skill Library (VPSL)."""
import pytest
import tempfile
from pathlib import Path
from pharmaforge.core.exceptions import VerificationFailureError
from pharmaforge.memory.skill_library import VerifiedProceduralSkillLibrary, SkillRetriever


@pytest.fixture
def temp_vpsl():
    with tempfile.TemporaryDirectory() as tmpdir:
        library = VerifiedProceduralSkillLibrary(storage_dir=tmpdir)
        yield library


def test_vpsl_register_valid_skill(temp_vpsl):
    """Test registering a skill that passes unit tests."""
    code = """
def multiply_by_two(val: float) -> float:
    return val * 2.0
"""
    unit_tests = """
def run_tests():
    assert multiply_by_two(5.0) == 10.0
    return True
"""
    skill = temp_vpsl.register_skill(
        name="multiply_by_two",
        description="Multiplies numerical value by 2.0",
        code=code,
        unit_tests=unit_tests,
        tags=["math", "scaling"],
    )
    assert skill.name == "multiply_by_two"
    assert skill.name in temp_vpsl.skills
    assert skill.reliability_score == 1.0


def test_vpsl_reject_failing_unit_tests(temp_vpsl):
    """Test that a skill failing unit tests is rejected and raises VerificationFailureError."""
    code = """
def buggy_add(a: int, b: int) -> int:
    return a - b  # Buggy implementation
"""
    unit_tests = """
def run_tests():
    assert buggy_add(2, 3) == 5  # Will fail
    return True
"""
    with pytest.raises(VerificationFailureError):
        temp_vpsl.register_skill(
            name="buggy_add",
            description="Adds two integers",
            code=code,
            unit_tests=unit_tests,
        )
    assert "buggy_add" not in temp_vpsl.skills


def test_vpsl_semantic_retrieval(temp_vpsl):
    """Test querying skills by semantic similarity."""
    results = temp_vpsl.retrieve_skills("evaluate drug-likeness Lipinski", top_k=2)
    assert len(results) >= 1
    best_skill, sim = results[0]
    assert "lipinski" in best_skill.name or "lipinski" in best_skill.description.lower()
    assert sim > 0.0


def test_vpsl_record_execution_outcome(temp_vpsl):
    """Test updating skill reliability metrics."""
    skill_name = "check_lipinski_rule"
    assert skill_name in temp_vpsl.skills

    temp_vpsl.record_execution_outcome(skill_name, success=True)
    temp_vpsl.record_execution_outcome(skill_name, success=False)
    skill = temp_vpsl.skills[skill_name]
    assert skill.success_count >= 2
    assert skill.failure_count == 1
    assert 0.0 < skill.reliability_score < 1.0
