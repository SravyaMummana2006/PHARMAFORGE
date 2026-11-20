"""Episodic memory buffer for tracking execution trajectories and agent deliberations."""
from __future__ import annotations

import datetime
from typing import Any
from pharmaforge.core.types import ExecutionTrace, MoleculeCandidate


class TrajectoryStep:
    """A single decision-action-observation step in an optimization episode."""

    def __init__(
        self,
        step_number: int,
        thought: str,
        action_name: str,
        action_input: dict[str, Any],
        trace: ExecutionTrace | None = None,
        resulting_candidates: list[MoleculeCandidate] | None = None,
    ):
        self.step_number = step_number
        self.thought = thought
        self.action_name = action_name
        self.action_input = action_input
        self.trace = trace
        self.resulting_candidates = resulting_candidates or []
        self.timestamp = datetime.datetime.now(datetime.timezone.utc)

    def to_dict(self) -> dict[str, Any]:
        return {
            "step_number": self.step_number,
            "thought": self.thought,
            "action_name": self.action_name,
            "action_input": self.action_input,
            "trace": self.trace.model_dump() if self.trace else None,
            "resulting_candidates": [c.model_dump() for c in self.resulting_candidates],
            "timestamp": self.timestamp.isoformat(),
        }


class EpisodicTrajectoryBuffer:
    """Maintains short-term trajectory history for the active task episode."""

    def __init__(self, task_id: str):
        self.task_id = task_id
        self.steps: list[TrajectoryStep] = []
        self.all_generated_molecules: list[MoleculeCandidate] = []
        self.pareto_front: list[MoleculeCandidate] = []

    def record_step(
        self,
        thought: str,
        action_name: str,
        action_input: dict[str, Any],
        trace: ExecutionTrace | None = None,
        candidates: list[MoleculeCandidate] | None = None,
    ) -> TrajectoryStep:
        """Append a new step to the trajectory history."""
        step_num = len(self.steps) + 1
        step = TrajectoryStep(
            step_number=step_num,
            thought=thought,
            action_name=action_name,
            action_input=action_input,
            trace=trace,
            resulting_candidates=candidates,
        )
        self.steps.append(step)
        if candidates:
            self.all_generated_molecules.extend(candidates)
        return step

    def get_successful_steps(self) -> list[TrajectoryStep]:
        """Return all steps whose tool execution was successful."""
        return [s for s in self.steps if s.trace and s.trace.success]

    def get_failed_steps(self) -> list[TrajectoryStep]:
        """Return all steps that encountered execution exceptions."""
        return [s for s in self.steps if s.trace and not s.trace.success]

    def export_summary(self) -> dict[str, Any]:
        """Export serialized summary of the trajectory."""
        return {
            "task_id": self.task_id,
            "total_steps": len(self.steps),
            "successful_steps": len(self.get_successful_steps()),
            "failed_steps": len(self.get_failed_steps()),
            "total_molecules_generated": len(self.all_generated_molecules),
            "steps": [s.to_dict() for s in self.steps],
        }

    def clear(self) -> None:
        """Reset the buffer for a new episode."""
        self.steps.clear()
        self.all_generated_molecules.clear()
        self.pareto_front.clear()
