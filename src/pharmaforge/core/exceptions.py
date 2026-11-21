"""Domain-specific exceptions for the PHARMAFORGE framework."""

class PharmaForgeError(Exception):
    """Base exception for all PHARMAFORGE errors."""
    pass


class RDKitUnavailableError(PharmaForgeError):
    """Raised when an operation requiring genuine chemoinformatics is executed without RDKit installed."""
    def __init__(self, operation: str = "Chemoinformatics computation"):
        super().__init__(
            f"RDKit is mandatory for {operation}, but is unavailable in the current Python environment. "
            "PHARMAFORGE strictly prohibits silent fallback to unscientific heuristic approximations."
        )



class VerificationFailureError(PharmaForgeError):
    """Raised when a molecule or code fails scientific test oracle verification."""
    def __init__(self, message: str, reasons: list[str] | None = None):
        super().__init__(message)
        self.reasons = reasons or []


class SandboxExecutionTimeoutError(PharmaForgeError):
    """Raised when code execution in the sandbox exceeds the allocated timeout."""
    def __init__(self, timeout_seconds: float):
        super().__init__(f"Sandbox execution timed out after {timeout_seconds} seconds.")
        self.timeout_seconds = timeout_seconds


class SandboxSecurityViolationError(PharmaForgeError):
    """Raised when sandboxed code attempts to access unauthorized system resources or forbidden modules."""
    def __init__(self, violation_detail: str):
        super().__init__(f"Security policy violation: {violation_detail}")
        self.violation_detail = violation_detail


class SkillCompilationError(PharmaForgeError):
    """Raised when synthesized skill code fails AST parsing or syntax validation."""
    def __init__(self, code: str, error_detail: str):
        super().__init__(f"Skill code failed compilation: {error_detail}")
        self.code = code
        self.error_detail = error_detail


class SkillRepairMaxAttemptsExceededError(PharmaForgeError):
    """Raised when the execution-trace self-repair engine fails to repair a skill within budget."""
    def __init__(self, skill_name: str, attempts: int, last_trace: str):
        super().__init__(f"Skill '{skill_name}' could not be repaired after {attempts} attempts.")
        self.skill_name = skill_name
        self.attempts = attempts
        self.last_trace = last_trace
