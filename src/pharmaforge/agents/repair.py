"""Execution-Trace Self-Repair Engine (ET-SRE): Automated program repair for scientific skills."""
from __future__ import annotations

import re
from typing import Any
from pharmaforge.core.exceptions import SkillRepairMaxAttemptsExceededError
from pharmaforge.sandbox.executor import IsolatedCodeExecutor


class ExecutionTraceSelfRepairEngine:
    """Automated Program Repair (APR) engine applying execution traceback reflection.
    
    Diagnoses syntax, runtime, and assertion errors from execution traces, applies
    localized code repairs, and validates against unit tests in the sandbox.
    """

    def __init__(
        self,
        executor: IsolatedCodeExecutor | None = None,
        max_repair_attempts: int = 3,
    ):
        self.executor = executor or IsolatedCodeExecutor()
        self.max_repair_attempts = max_repair_attempts

    def repair_code(
        self,
        skill_name: str,
        failing_code: str,
        unit_tests: str,
        traceback_str: str,
    ) -> tuple[str, list[str]]:
        """Attempt to iteratively repair failing code using execution feedback.
        
        Args:
            skill_name: Name of the failing skill.
            failing_code: Current buggy function implementation.
            unit_tests: The verification test suite.
            traceback_str: Stack trace from the failed run.
            
        Returns:
            Tuple of (repaired_code: str, repair_changelog: list[str])
        """
        current_code = failing_code
        changelog: list[str] = []

        for attempt in range(1, self.max_repair_attempts + 1):
            diagnosis = self._diagnose_fault(current_code, traceback_str)
            patched_code = self._apply_patch(current_code, diagnosis)
            changelog.append(f"Attempt {attempt}: Applied patch for '{diagnosis['issue_type']}'")

            # Test patched code in sandbox
            verification_code = f"""
{patched_code}

# --- UNIT TESTS ---
{unit_tests}
"""
            trace = self.executor.execute_function(
                function_code=verification_code,
                function_name="run_tests",
                input_kwargs={},
            )

            if trace.success:
                changelog.append(f"Success: Unit tests passed on attempt {attempt}.")
                return patched_code, changelog

            # Update traceback for next iteration
            traceback_str = trace.traceback_str or "Test failure"
            current_code = patched_code

        raise SkillRepairMaxAttemptsExceededError(
            skill_name=skill_name,
            attempts=self.max_repair_attempts,
            last_trace=traceback_str,
        )

    def _diagnose_fault(self, code: str, traceback_str: str) -> dict[str, Any]:
        """Analyze traceback and identify the error category and target location."""
        diag: dict[str, Any] = {"issue_type": "unknown", "line": None, "detail": traceback_str}

        if "NameError:" in traceback_str:
            match = re.search(r"NameError: name '(\w+)' is not defined", traceback_str)
            diag["issue_type"] = "undefined_name"
            diag["name"] = match.group(1) if match else ""
        elif "TypeError: 'NoneType' object is not subscriptable" in traceback_str:
            diag["issue_type"] = "none_guard_missing"
        elif "IndexError:" in traceback_str:
            diag["issue_type"] = "index_out_of_bounds"
        elif "AssertionError" in traceback_str:
            diag["issue_type"] = "assertion_failure"
        elif "SyntaxError:" in traceback_str:
            diag["issue_type"] = "syntax_error"

        return diag

    def _apply_patch(self, code: str, diagnosis: dict[str, Any]) -> str:
        """Apply targeted automated patch based on diagnosed fault."""
        issue = diagnosis["issue_type"]

        if issue == "none_guard_missing" or issue == "index_out_of_bounds":
            # Add safe guard clause at beginning of function
            lines = code.splitlines()
            for i, line in enumerate(lines):
                if line.strip().startswith("def "):
                    indent = "    "
                    guard = f"{indent}if not locals().get('smiles', True):\n{indent}    return []"
                    lines.insert(i + 1, guard)
                    return "\n".join(lines)

        elif issue == "undefined_name":
            name = diagnosis.get("name", "")
            if name in {"Chem", "AllChem"}:
                return f"from rdkit import {name}\n" + code
            elif name == "np":
                return "import numpy as np\n" + code
            elif name in {"math", "re", "json", "numpy"}:
                return f"import {name}\n" + code

        elif issue == "assertion_failure":
            # If function returns None instead of expected type, ensure return value exists
            if "return " not in code:
                return code + "\n    return []"

        # General repair heuristic: append fallback return if empty
        if not code.endswith("\n"):
            code += "\n"
        return code
