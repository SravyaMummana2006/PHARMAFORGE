"""Base agent definitions and interfaces for PHARMAFORGE."""
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any
from pharmaforge.core.types import LeadOptimizationTask


class BaseScientificAgent(ABC):
    """Abstract base class for all scientific reasoning agents."""

    def __init__(self, name: str):
        self.name = name

    @abstractmethod
    def execute_step(self, task: LeadOptimizationTask, context: dict[str, Any]) -> dict[str, Any]:
        """Execute a single planning, synthesis, or evaluation step."""
        pass
