"""Isolated execution sandbox for running agent-synthesized code with AST checks and timeouts."""
from __future__ import annotations

import ast
import json
import subprocess
import sys
import tempfile
import time
import uuid
from pathlib import Path
from typing import Any
from pharmaforge.core.exceptions import (
    SandboxExecutionTimeoutError,
    SandboxSecurityViolationError,
    SkillCompilationError,
)
from pharmaforge.core.types import ExecutionTrace

# Disallowed dangerous modules and attributes
FORBIDDEN_MODULES: set[str] = {
    "os", "sys", "subprocess", "socket", "shutil", "ctypes", "winreg",
    "pty", "resource", "builtins.eval", "builtins.exec", "builtins.compile",
    "webbrowser", "http.server", "urllib.request", "ftplib", "smtplib"
}

FORBIDDEN_CALLS: set[str] = {
    "eval", "exec", "compile", "__import__", "globals", "locals"
}


class ASTSecurityInspector(ast.NodeVisitor):
    """Inspects Python AST to reject unauthorized system calls or malicious imports."""

    def __init__(self, allowed_modules: set[str] | None = None):
        self.allowed_modules = allowed_modules
        self.violations: list[str] = []

    def visit_Import(self, node: ast.Import) -> None:
        for alias in node.names:
            base_module = alias.name.split(".")[0]
            if base_module in FORBIDDEN_MODULES:
                self.violations.append(f"Forbidden module import: '{alias.name}'")
            elif self.allowed_modules and base_module not in self.allowed_modules:
                self.violations.append(f"Module not in allowlist: '{alias.name}'")
        self.generic_visit(node)

    def visit_ImportFrom(self, node: ast.ImportFrom) -> None:
        if node.module:
            base_module = node.module.split(".")[0]
            if base_module in FORBIDDEN_MODULES:
                self.violations.append(f"Forbidden module import from: '{node.module}'")
            elif self.allowed_modules and base_module not in self.allowed_modules:
                self.violations.append(f"Module not in allowlist: '{node.module}'")
        self.generic_visit(node)

    def visit_Call(self, node: ast.Call) -> None:
        if isinstance(node.func, ast.Name):
            if node.func.id in FORBIDDEN_CALLS:
                self.violations.append(f"Forbidden function call: '{node.func.id}()'")
        self.generic_visit(node)


class IsolatedCodeExecutor:
    """Safely executes dynamically synthesized scientific code in an isolated subprocess."""

    def __init__(
        self,
        timeout_seconds: float = 15.0,
        allowed_modules: set[str] | None = None,
    ):
        self.timeout_seconds = timeout_seconds
        self.allowed_modules = allowed_modules or {
            "math", "numpy", "scipy", "rdkit", "torch", "collections",
            "itertools", "typing", "pydantic", "re", "json", "time"
        }

    def validate_code_ast(self, code_str: str) -> None:
        """Validate code syntax and check against security policies."""
        try:
            tree = ast.parse(code_str)
        except SyntaxError as e:
            raise SkillCompilationError(code_str, f"SyntaxError at line {e.lineno}: {e.msg}") from e

        inspector = ASTSecurityInspector(allowed_modules=self.allowed_modules)
        inspector.visit(tree)

        if inspector.violations:
            violation_msg = "; ".join(inspector.violations)
            raise SandboxSecurityViolationError(violation_msg)

    def execute_function(
        self,
        function_code: str,
        function_name: str,
        input_kwargs: dict[str, Any],
        timeout: float | None = None,
    ) -> ExecutionTrace:
        """Execute a Python function definition with arguments and capture structured output.
        
        Args:
            function_code: The Python code containing the function definition.
            function_name: The name of the function to invoke.
            input_kwargs: Keyword arguments to pass into the function.
            timeout: Optional override for execution timeout.
            
        Returns:
            ExecutionTrace record with stdout, return value, or traceback.
        """
        execution_id = str(uuid.uuid4())[:8]
        effective_timeout = timeout or self.timeout_seconds
        start_time = time.perf_counter()

        # Step 1: AST Security Inspection
        try:
            self.validate_code_ast(function_code)
        except (SkillCompilationError, SandboxSecurityViolationError) as e:
            return ExecutionTrace(
                execution_id=execution_id,
                tool_or_skill_name=function_name,
                input_kwargs=input_kwargs,
                success=False,
                exception_type=type(e).__name__,
                traceback_str=str(e),
                execution_time_seconds=time.perf_counter() - start_time,
            )

        # Step 2: Build execution harness script
        # Serialize inputs to JSON
        serialized_inputs = json.dumps(input_kwargs)
        runner_script = f"""
import json
import sys
import traceback

{function_code}

def __runner():
    kwargs = json.loads({json.dumps(serialized_inputs)})
    try:
        result = {function_name}(**kwargs)
        print("__RESULT_START__")
        print(json.dumps(result))
        print("__RESULT_END__")
    except Exception as e:
        sys.stderr.write(traceback.format_exc())
        sys.exit(1)

if __name__ == "__main__":
    __runner()
"""

        # Step 3: Run in isolated temporary subprocess
        with tempfile.NamedTemporaryFile("w", suffix=".py", delete=False) as temp_file:
            temp_path = Path(temp_file.name)
            temp_file.write(runner_script)

        try:
            process = subprocess.run(
                [sys.executable, str(temp_path)],
                capture_output=True,
                text=True,
                timeout=effective_timeout,
            )
            elapsed = time.perf_counter() - start_time

            stdout = process.stdout
            stderr = process.stderr

            if process.returncode != 0:
                return ExecutionTrace(
                    execution_id=execution_id,
                    tool_or_skill_name=function_name,
                    input_kwargs=input_kwargs,
                    stdout=stdout,
                    stderr=stderr,
                    success=False,
                    exception_type="SubprocessRuntimeError",
                    traceback_str=stderr or stdout,
                    execution_time_seconds=elapsed,
                )

            # Parse return value
            return_val = None
            if "__RESULT_START__" in stdout and "__RESULT_END__" in stdout:
                raw_json = stdout.split("__RESULT_START__")[1].split("__RESULT_END__")[0].strip()
                try:
                    return_val = json.loads(raw_json)
                except json.JSONDecodeError:
                    return_val = raw_json

            return ExecutionTrace(
                execution_id=execution_id,
                tool_or_skill_name=function_name,
                input_kwargs=input_kwargs,
                stdout=stdout.split("__RESULT_START__")[0],
                stderr=stderr,
                return_value=return_val,
                success=True,
                execution_time_seconds=elapsed,
            )

        except subprocess.TimeoutExpired:
            return ExecutionTrace(
                execution_id=execution_id,
                tool_or_skill_name=function_name,
                input_kwargs=input_kwargs,
                success=False,
                exception_type="SandboxExecutionTimeoutError",
                traceback_str=f"Execution timed out after {effective_timeout}s.",
                execution_time_seconds=effective_timeout,
            )
        finally:
            if temp_path.exists():
                try:
                    temp_path.unlink()
                except OSError:
                    pass
