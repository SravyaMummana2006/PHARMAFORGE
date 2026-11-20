"""Unit tests for the IsolatedCodeExecutor sandbox and AST security policies."""
import pytest
from pharmaforge.core.exceptions import SandboxSecurityViolationError, SkillCompilationError
from pharmaforge.sandbox.executor import IsolatedCodeExecutor


@pytest.fixture
def executor():
    return IsolatedCodeExecutor(timeout_seconds=5.0)


def test_sandbox_successful_execution(executor):
    """Test standard valid Python code execution in sandbox."""
    code = """
def compute_sum(a: int, b: int) -> int:
    return a + b
"""
    trace = executor.execute_function(
        function_code=code,
        function_name="compute_sum",
        input_kwargs={"a": 10, "b": 25},
    )
    assert trace.success is True
    assert trace.return_value == 35
    assert trace.execution_time_seconds > 0.0


def test_sandbox_security_rejection_forbidden_import(executor):
    """Test that unauthorized system module imports are blocked by AST inspection."""
    malicious_code = """
import os

def malicious_call():
    return os.listdir('.')
"""
    trace = executor.execute_function(
        function_code=malicious_code,
        function_name="malicious_call",
        input_kwargs={},
    )
    assert trace.success is False
    assert "SandboxSecurityViolationError" in (trace.exception_type or "")


def test_sandbox_security_rejection_forbidden_call(executor):
    """Test that forbidden built-in calls like eval() are intercepted."""
    code = """
def sneaky_call():
    return eval("2 + 2")
"""
    trace = executor.execute_function(
        function_code=code,
        function_name="sneaky_call",
        input_kwargs={},
    )
    assert trace.success is False
    assert "Forbidden function call: 'eval()'" in (trace.traceback_str or "")


def test_sandbox_timeout_handling(executor):
    """Test that long-running operations are terminated when timeout is reached."""
    slow_code = """
import time

def infinite_loop():
    time.sleep(10)
    return True
"""
    # Use short 1-second timeout
    trace = executor.execute_function(
        function_code=slow_code,
        function_name="infinite_loop",
        input_kwargs={},
        timeout=1.0,
    )
    assert trace.success is False
    assert "Timeout" in (trace.exception_type or "")


def test_sandbox_syntax_error(executor):
    """Test that syntax errors are intercepted gracefully."""
    bad_syntax = """
def broken_syntax(
    return 42
"""
    trace = executor.execute_function(
        function_code=bad_syntax,
        function_name="broken_syntax",
        input_kwargs={},
    )
    assert trace.success is False
    assert "SkillCompilationError" in (trace.exception_type or "")
