"""Main Orchestrator: Closed-loop multi-agent coordinator for pharmaceutical discovery."""
from __future__ import annotations

import time
from typing import Any
from pharmaforge.agents.planner import ScientificDiscoveryPlanner
from pharmaforge.agents.repair import ExecutionTraceSelfRepairEngine
from pharmaforge.agents.synthesizer import SkillSynthesizer
from pharmaforge.core.types import (
    LeadOptimizationTask,
    MoleculeCandidate,
    OptimizationObjective,
    ObjectiveDirection,
)
from pharmaforge.memory.episodic import EpisodicTrajectoryBuffer
from pharmaforge.memory.skill_library import VerifiedProceduralSkillLibrary
from pharmaforge.sandbox.executor import IsolatedCodeExecutor
from pharmaforge.verifier.oracle import ScientificTestOracle


class PharmaForgeOrchestrator:
    """Coordinates Planner, VPSL, Synthesizer, Sandbox, STO-IV, and Self-Repair."""

    def __init__(
        self,
        oracle: ScientificTestOracle | None = None,
        skill_library: VerifiedProceduralSkillLibrary | None = None,
        executor: IsolatedCodeExecutor | None = None,
        repair_engine: ExecutionTraceSelfRepairEngine | None = None,
    ):
        self.oracle = oracle or ScientificTestOracle()
        self.executor = executor or IsolatedCodeExecutor()
        self.library = skill_library or VerifiedProceduralSkillLibrary(
            oracle=self.oracle,
            executor=self.executor,
        )
        self.planner = ScientificDiscoveryPlanner(skill_library=self.library, oracle=self.oracle)
        self.synthesizer = SkillSynthesizer(skill_library=self.library)
        self.repair_engine = repair_engine or ExecutionTraceSelfRepairEngine(executor=self.executor)

    def run_optimization(self, task: LeadOptimizationTask) -> dict[str, Any]:
        """Execute closed-loop multi-objective lead optimization on a target task."""
        start_time = time.perf_counter()
        buffer = EpisodicTrajectoryBuffer(task_id=task.task_id)

        # Step 1: Ingest and verify initial lead molecule
        try:
            initial_candidate = self.oracle.create_candidate(
                smiles=task.initial_smiles,
                transformation="initial_lead",
                affinity_score=6.0,  # Baseline target affinity proxy
            )
        except Exception as e:
            return {
                "success": False,
                "error": f"Initial lead SMILES failed verification: {e}",
                "pareto_front": [],
            }

        candidates_pool: list[MoleculeCandidate] = [initial_candidate]
        current_lead = initial_candidate
        newly_learned_skills: list[str] = []

        # Step 2: Iterative optimization loop
        for iteration in range(1, task.max_iterations + 1):
            actions = self.planner.plan_next_actions(task, current_lead)
            if not actions or actions[0].get("action") == "terminate":
                buffer.record_step(
                    thought="All multi-objective constraints satisfied or no viable actions.",
                    action_name="terminate",
                    action_input={},
                )
                break

            action = actions[0]
            action_type = action["action"]
            executed_skill_name = None
            generated_smiles_list: list[str] = []

            if action_type == "execute_skill":
                executed_skill_name = action["skill_name"]
                skill = self.library.skills.get(executed_skill_name)
                if skill:
                    trace = self.executor.execute_function(
                        function_code=skill.code,
                        function_name=skill.name,
                        input_kwargs=action["input_kwargs"],
                    )
                    # If execution failed, trigger automated program repair
                    if not trace.success:
                        try:
                            repaired_code, changelog = self.repair_engine.repair_code(
                                skill_name=skill.name,
                                failing_code=skill.code,
                                unit_tests=skill.unit_tests,
                                traceback_str=trace.traceback_str or "",
                            )
                            # Update repaired skill in library
                            skill.code = repaired_code
                            self.library.save_to_disk()
                            # Re-execute repaired skill
                            trace = self.executor.execute_function(
                                function_code=repaired_code,
                                function_name=skill.name,
                                input_kwargs=action["input_kwargs"],
                            )
                        except Exception:
                            pass

                    self.library.record_execution_outcome(skill.name, trace.success)
                    if trace.success and isinstance(trace.return_value, list):
                        generated_smiles_list = trace.return_value

            elif action_type == "synthesize_skill":
                target_prop = action["target_property"]
                try:
                    new_skill = self.synthesizer.synthesize_transformation_skill(target_property=target_prop)
                    newly_learned_skills.append(new_skill.name)
                    executed_skill_name = new_skill.name
                    trace = self.executor.execute_function(
                        function_code=new_skill.code,
                        function_name=new_skill.name,
                        input_kwargs=action["input_kwargs"],
                    )
                    self.library.record_execution_outcome(new_skill.name, trace.success)
                    if trace.success and isinstance(trace.return_value, list):
                        generated_smiles_list = trace.return_value
                except Exception as e:
                    trace = None

            # Step 3: Independent Verification (STO-IV) of all generated candidates
            valid_new_candidates: list[MoleculeCandidate] = []
            for smiles in generated_smiles_list:
                report = self.oracle.verify_smiles(smiles)
                if report.is_valid:
                    try:
                        # Compute affinity surrogate score
                        affinity_delta = 0.3 if "F" in smiles or "N" in smiles else -0.1
                        cand = self.oracle.create_candidate(
                            smiles=smiles,
                            parent_smiles=current_lead.smiles,
                            transformation=executed_skill_name or "modification",
                            affinity_score=round((current_lead.affinity_score or 6.0) + affinity_delta, 2),
                        )
                        valid_new_candidates.append(cand)
                    except ValueError:
                        pass

            if valid_new_candidates:
                candidates_pool.extend(valid_new_candidates)
                # Select top candidate as next lead
                current_lead = valid_new_candidates[0]

            buffer.record_step(
                thought=f"Iteration {iteration}: Targeting {action.get('target_property', 'lead improvement')}.",
                action_name=action_type,
                action_input=action,
                trace=trace,
                candidates=valid_new_candidates,
            )

        # Step 4: Extract non-dominated Pareto frontier
        pareto_front = self._compute_pareto_front(candidates_pool, task.objectives)
        elapsed = time.perf_counter() - start_time

        return {
            "success": True,
            "task_id": task.task_id,
            "total_candidates_generated": len(candidates_pool),
            "pareto_front_size": len(pareto_front),
            "pareto_front": [c.model_dump() for c in pareto_front],
            "newly_learned_skills": newly_learned_skills,
            "wall_clock_seconds": round(elapsed, 2),
            "trajectory_summary": buffer.export_summary(),
        }

    def _compute_pareto_front(
        self,
        candidates: list[MoleculeCandidate],
        objectives: list[OptimizationObjective],
    ) -> list[MoleculeCandidate]:
        """Compute the set of non-dominated Pareto-optimal candidates across all objectives."""
        if not candidates:
            return []

        pareto: list[MoleculeCandidate] = []
        for i, c1 in enumerate(candidates):
            dominated = False
            for j, c2 in enumerate(candidates):
                if i == j:
                    continue
                # Check if c2 dominates c1
                c2_better_or_equal = True
                c2_strictly_better = False

                for obj in objectives:
                    v1 = c1.get_property(obj.name)
                    v2 = c2.get_property(obj.name)
                    if v1 is None or v2 is None:
                        continue

                    if obj.direction == ObjectiveDirection.MAXIMIZE:
                        if v2 < v1:
                            c2_better_or_equal = False
                            break
                        if v2 > v1:
                            c2_strictly_better = True
                    else:
                        if v2 > v1:
                            c2_better_or_equal = False
                            break
                        if v2 < v1:
                            c2_strictly_better = True

                if c2_better_or_equal and c2_strictly_better:
                    dominated = True
                    break

            if not dominated and c1 not in pareto:
                pareto.append(c1)

        return pareto
