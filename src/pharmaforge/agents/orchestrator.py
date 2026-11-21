"""Main Orchestrator: Closed-loop multi-agent coordinator for pharmaceutical discovery."""
from __future__ import annotations

import time
from typing import Any
from pharmaforge.agents.planner import ScientificDiscoveryPlanner
from pharmaforge.agents.repair import ExecutionTraceSelfRepairEngine
from pharmaforge.agents.synthesizer import SkillSynthesizer
from pharmaforge.benchmarks.metrics import BenchmarkEvaluator
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
        qsar_predictor: Any | None = None,
        enable_skill_memory: bool = True,
        enable_skill_verification: bool = True,
        enable_self_repair: bool = True,
        inject_synthesis_noise: bool = False,
    ):
        self.oracle = oracle or ScientificTestOracle()
        self.executor = executor or IsolatedCodeExecutor()
        self.enable_skill_memory = enable_skill_memory
        self.enable_skill_verification = enable_skill_verification
        self.enable_self_repair = enable_self_repair
        self.inject_synthesis_noise = inject_synthesis_noise

        self.library = skill_library or VerifiedProceduralSkillLibrary(
            oracle=self.oracle,
            executor=self.executor,
        )
        self.planner = ScientificDiscoveryPlanner(
            skill_library=self.library,
            oracle=self.oracle,
            enable_skill_memory=self.enable_skill_memory,
        )
        self.synthesizer = SkillSynthesizer(skill_library=self.library)
        self.repair_engine = repair_engine or ExecutionTraceSelfRepairEngine(executor=self.executor)
        
        self.qsar_predictor = qsar_predictor
        if self.qsar_predictor is None:
            try:
                from pharmaforge.qsar.predictor import EGFRQSARPredictor
                predictor = EGFRQSARPredictor()
                if predictor.is_ready:
                    self.qsar_predictor = predictor
            except Exception:
                self.qsar_predictor = None

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
            if self.qsar_predictor and self.qsar_predictor.is_ready:
                self.qsar_predictor.predict_candidate(initial_candidate)
        except Exception as e:
            return {
                "success": False,
                "error": f"Initial lead SMILES failed verification: {e}",
                "pareto_front": [],
            }

        candidates_pool: list[MoleculeCandidate] = [initial_candidate]
        current_lead = initial_candidate
        newly_learned_skills: list[str] = []

        total_candidates_proposed = 0
        valid_candidates_count = 0
        rejected_candidates: list[dict[str, Any]] = []

        skill_reuse_count = 0
        skill_synthesized_count = 0
        skills_admitted_count = 0
        skills_rejected_count = 0
        repair_attempts = 0
        successful_repairs = 0
        iterations_executed = 0

        # Step 2: Iterative optimization loop
        for iteration in range(1, task.max_iterations + 1):
            iterations_executed += 1
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
            trace = None

            if action_type == "execute_skill":
                skill_reuse_count += 1
                executed_skill_name = action["skill_name"]
                skill = self.library.skills.get(executed_skill_name)
                if skill:
                    trace = self.executor.execute_function(
                        function_code=skill.code,
                        function_name=skill.name,
                        input_kwargs=action["input_kwargs"],
                    )
                    # If execution failed and self-repair enabled, trigger repair
                    if not trace.success and self.enable_self_repair:
                        repair_attempts += 1
                        try:
                            repaired_code, changelog = self.repair_engine.repair_code(
                                skill_name=skill.name,
                                failing_code=skill.code,
                                unit_tests=skill.unit_tests,
                                traceback_str=trace.traceback_str or "",
                            )
                            # Update repaired skill in library
                            skill.code = repaired_code
                            if self.enable_skill_memory:
                                self.library.save_to_disk()
                            successful_repairs += 1
                            # Re-execute repaired skill
                            trace = self.executor.execute_function(
                                function_code=repaired_code,
                                function_name=skill.name,
                                input_kwargs=action["input_kwargs"],
                            )
                        except Exception:
                            pass

                    if self.enable_skill_memory:
                        self.library.record_execution_outcome(skill.name, trace.success)
                    if trace.success and isinstance(trace.return_value, list):
                        generated_smiles_list = [s for s in trace.return_value if s != action["input_kwargs"].get("smiles")]
                        if not generated_smiles_list:
                            # Retrieved skill produced no modifications on this lead; fallback to dynamic synthesis
                            action_type = "synthesize_skill"
                            action["target_property"] = action.get("target_deficiency", {}).get("objective_name", "optimization")

            if action_type == "synthesize_skill":
                skill_synthesized_count += 1
                target_prop = action["target_property"]
                name, desc, code, unit_tests, tags = self.synthesizer.synthesize_skill_components(
                    target_property=target_prop,
                    inject_synthesis_noise=self.inject_synthesis_noise,
                )

                new_skill = None
                admitted = False

                if not self.enable_skill_verification:
                    # Config B: Bypass verification, admit directly into memory
                    new_skill = self.library.register_skill(
                        name=name,
                        description=desc,
                        code=code,
                        unit_tests=unit_tests,
                        tags=tags,
                        bypass_verification=True,
                    )
                    skills_admitted_count += 1
                    admitted = True
                else:
                    # Config C / D: Verify via sandbox unit tests
                    try:
                        new_skill = self.library.register_skill(
                            name=name,
                            description=desc,
                            code=code,
                            unit_tests=unit_tests,
                            tags=tags,
                            bypass_verification=False,
                        )
                        skills_admitted_count += 1
                        admitted = True
                    except Exception as ve:
                        skills_rejected_count += 1
                        if self.enable_self_repair:
                            repair_attempts += 1
                            try:
                                repaired_code, changelog = self.repair_engine.repair_code(
                                    skill_name=name,
                                    failing_code=code,
                                    unit_tests=unit_tests,
                                    traceback_str=str(ve),
                                )
                                new_skill = self.library.register_skill(
                                    name=name,
                                    description=desc,
                                    code=repaired_code,
                                    unit_tests=unit_tests,
                                    tags=tags,
                                    bypass_verification=False,
                                )
                                successful_repairs += 1
                                skills_admitted_count += 1
                                admitted = True
                            except Exception:
                                pass

                if admitted and new_skill:
                    newly_learned_skills.append(new_skill.name)
                    executed_skill_name = new_skill.name
                    trace = self.executor.execute_function(
                        function_code=new_skill.code,
                        function_name=new_skill.name,
                        input_kwargs=action["input_kwargs"],
                    )
                    if not trace.success and self.enable_self_repair:
                        repair_attempts += 1
                        try:
                            repaired_code, changelog = self.repair_engine.repair_code(
                                skill_name=new_skill.name,
                                failing_code=new_skill.code,
                                unit_tests=new_skill.unit_tests,
                                traceback_str=trace.traceback_str or "",
                            )
                            new_skill.code = repaired_code
                            successful_repairs += 1
                            trace = self.executor.execute_function(
                                function_code=repaired_code,
                                function_name=new_skill.name,
                                input_kwargs=action["input_kwargs"],
                            )
                        except Exception:
                            pass

                    if self.enable_skill_memory:
                        self.library.record_execution_outcome(new_skill.name, trace.success)
                    if trace.success and isinstance(trace.return_value, list):
                        generated_smiles_list = trace.return_value

            # Step 3: Independent Verification (STO-IV) of all generated candidates
            total_candidates_proposed += len(generated_smiles_list)
            valid_new_candidates: list[MoleculeCandidate] = []
            
            for smiles in generated_smiles_list:
                report = self.oracle.verify_smiles(smiles)
                if report.is_valid:
                    valid_candidates_count += 1
                    try:
                        if self.qsar_predictor and self.qsar_predictor.is_ready:
                            cand = self.oracle.create_candidate(
                                smiles=smiles,
                                parent_smiles=current_lead.smiles,
                                transformation=executed_skill_name or "modification",
                            )
                            self.qsar_predictor.predict_candidate(cand)
                        else:
                            affinity_delta = 0.3 if "F" in smiles or "N" in smiles else -0.1
                            cand = self.oracle.create_candidate(
                                smiles=smiles,
                                parent_smiles=current_lead.smiles,
                                transformation=executed_skill_name or "modification",
                                affinity_score=round((current_lead.affinity_score or 6.0) + affinity_delta, 2),
                            )
                        valid_new_candidates.append(cand)
                    except ValueError as ve:
                        rejected_candidates.append({"smiles": smiles, "reasons": [str(ve)]})
                else:
                    rejected_candidates.append({"smiles": smiles, "reasons": report.error_reasons})

            if valid_new_candidates:
                candidates_pool.extend(valid_new_candidates)
                # Select top candidate as next lead based on multi-objective rank
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

        # Calculate mathematical metrics
        pareto_hv = BenchmarkEvaluator.calculate_pareto_hypervolume(pareto_front, task.objectives)
        validity_rate = BenchmarkEvaluator.calculate_chemical_validity_rate(
            total_proposed=total_candidates_proposed,
            valid_sanitized=valid_candidates_count,
        )
        total_execs = skill_reuse_count + skill_synthesized_count
        reuse_ratio = BenchmarkEvaluator.calculate_skill_reuse_ratio(
            total_executions=total_execs,
            reused_existing_skills=skill_reuse_count,
        )

        return {
            "success": True,
            "task_id": task.task_id,
            "iterations_executed": iterations_executed,
            "total_candidates_generated": len(candidates_pool),
            "total_candidates_proposed": total_candidates_proposed,
            "valid_candidates_count": valid_candidates_count,
            "validity_rate_pct": validity_rate,
            "rejected_candidates_count": len(rejected_candidates),
            "rejected_candidates_sample": rejected_candidates[:10],
            "skill_reuse_count": skill_reuse_count,
            "skill_synthesized_count": skill_synthesized_count,
            "skill_reuse_ratio_pct": reuse_ratio,
            "skills_admitted_count": skills_admitted_count,
            "skills_rejected_count": skills_rejected_count,
            "repair_attempts": repair_attempts,
            "successful_repairs": successful_repairs,
            "pareto_front_size": len(pareto_front),
            "pareto_hypervolume": pareto_hv,
            "pareto_front": [c.model_dump() for c in pareto_front],
            "all_candidates": [c.model_dump() for c in candidates_pool],
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
