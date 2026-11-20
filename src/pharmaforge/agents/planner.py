"""Scientific Discovery Planner: Goal decomposition and skill query formulation."""
from __future__ import annotations

from typing import Any
from pharmaforge.core.types import LeadOptimizationTask, MoleculeCandidate, ObjectiveDirection
from pharmaforge.memory.skill_library import VerifiedProceduralSkillLibrary
from pharmaforge.verifier.oracle import ScientificTestOracle


class ScientificDiscoveryPlanner:
    """Decomposes multi-objective pharmaceutical optimization goals into verifiable sub-tasks."""

    def __init__(
        self,
        skill_library: VerifiedProceduralSkillLibrary,
        oracle: ScientificTestOracle | None = None,
    ):
        self.library = skill_library
        self.oracle = oracle or ScientificTestOracle()

    def analyze_lead_deficiencies(
        self,
        task: LeadOptimizationTask,
        current_candidate: MoleculeCandidate,
    ) -> list[dict[str, Any]]:
        """Identify which target objectives are currently violated by the lead compound."""
        deficiencies = []

        for obj in task.objectives:
            val = current_candidate.get_property(obj.name)
            if val is None:
                deficiencies.append({
                    "objective_name": obj.name,
                    "target_direction": obj.direction.value,
                    "threshold": obj.threshold,
                    "observed_value": None,
                    "status": "unmeasured",
                    "priority": 1.0,
                })
            elif not obj.is_satisfied(val):
                gap = abs(val - (obj.threshold or 0.0))
                deficiencies.append({
                    "objective_name": obj.name,
                    "target_direction": obj.direction.value,
                    "threshold": obj.threshold,
                    "observed_value": val,
                    "status": "violated",
                    "priority": obj.weight * (1.0 + gap),
                })

        # Sort deficiencies by priority descending
        deficiencies.sort(key=lambda x: x["priority"], reverse=True)
        return deficiencies

    def plan_next_actions(
        self,
        task: LeadOptimizationTask,
        current_candidate: MoleculeCandidate,
    ) -> list[dict[str, Any]]:
        """Formulate next concrete procedural actions to address priority deficiencies."""
        deficiencies = self.analyze_lead_deficiencies(task, current_candidate)
        actions = []

        if not deficiencies:
            actions.append({
                "action": "terminate",
                "reason": "All multi-objective constraints successfully satisfied.",
            })
            return actions

        top_deficiency = deficiencies[0]
        obj_name = top_deficiency["objective_name"]
        direction = top_deficiency["target_direction"]

        # Formulate query for skill library
        query = f"optimize {obj_name} {direction} molecular lead"
        retrieved = self.library.retrieve_skills(query, top_k=2, threshold=0.15)

        if retrieved:
            best_skill, sim = retrieved[0]
            actions.append({
                "action": "execute_skill",
                "skill_name": best_skill.name,
                "similarity": sim,
                "target_deficiency": top_deficiency,
                "input_kwargs": {"smiles": current_candidate.smiles},
            })
        else:
            actions.append({
                "action": "synthesize_skill",
                "target_property": obj_name,
                "target_deficiency": top_deficiency,
                "input_kwargs": {"smiles": current_candidate.smiles},
            })

        return actions
