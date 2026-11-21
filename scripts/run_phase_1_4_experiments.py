"""Run Phase 1.4 Benchmark Lead Optimization and Controlled Ablation Experiments."""
import os
import sys
from pathlib import Path

# Ensure src is in python path
src_dir = Path(__file__).resolve().parent.parent / "src"
if str(src_dir) not in sys.path:
    sys.path.insert(0, str(src_dir))

import json
import time
import shutil
from typing import Any

from pharmaforge.agents.orchestrator import PharmaForgeOrchestrator
from pharmaforge.core.types import LeadOptimizationTask, OptimizationObjective, ObjectiveDirection
from pharmaforge.memory.skill_library import VerifiedProceduralSkillLibrary
from pharmaforge.verifier.oracle import ScientificTestOracle
from pharmaforge.sandbox.executor import IsolatedCodeExecutor
from pharmaforge.qsar.predictor import EGFRQSARPredictor

BENCHMARK_LEADS = {
    "CHEMBL289959": {
        "name": "CHEMBL289959",
        "smiles": "c1ccc(Nc2ncnc3ccccc23)cc1",
        "description": "4-anilinoquinazoline with unsubstituted phenyl (Rewcastle et al. 1996; EGFR IC50 = 344 nM)",
        "experimental_pchembl": 6.463,
    },
    "CHEMBL290096": {
        "name": "CHEMBL290096",
        "smiles": "Brc1cccc(Nc2ncnc3ccccc23)c1",
        "description": "4-(3-bromoanilino)quinazoline with 3-bromophenyl (Rewcastle et al. 1996; EGFR IC50 = 27 nM)",
        "experimental_pchembl": 7.569,
    },
}

def create_task(lead_key: str, max_iterations: int = 5) -> LeadOptimizationTask:
    lead = BENCHMARK_LEADS[lead_key]
    return LeadOptimizationTask(
        task_id=f"EGFR-OPT-{lead_key}",
        target_name="EGFR",
        target_class="Kinase",
        initial_smiles=lead["smiles"],
        objectives=[
            OptimizationObjective(name="affinity", direction=ObjectiveDirection.MAXIMIZE, threshold=7.0, weight=2.0),
            OptimizationObjective(name="qed", direction=ObjectiveDirection.MAXIMIZE, threshold=0.60, weight=1.0),
            OptimizationObjective(name="logp", direction=ObjectiveDirection.MINIMIZE, threshold=3.5, weight=1.0),
            OptimizationObjective(name="sascore", direction=ObjectiveDirection.MINIMIZE, threshold=3.5, weight=1.0),
        ],
        max_iterations=max_iterations,
    )

def run_experiment_suite():
    artifacts_dir = Path("artifacts/benchmark_runs")
    artifacts_dir.mkdir(parents=True, exist_ok=True)
    
    oracle = ScientificTestOracle()
    predictor = EGFRQSARPredictor()
    assert predictor.is_ready, "QSAR predictor model artifact must be loaded!"
    
    ablation_configs = {
        "A_no_memory_no_repair": {
            "name": "A: No persistent skill memory and no repair",
            "enable_skill_memory": False,
            "enable_skill_verification": False,
            "enable_self_repair": False,
            "inject_synthesis_noise": False,
        },
        "B_unverified_memory_no_repair": {
            "name": "B: Persistent but unverified skill memory",
            "enable_skill_memory": True,
            "enable_skill_verification": False,
            "enable_self_repair": False,
            "inject_synthesis_noise": True,
        },
        "C_verified_memory_no_repair": {
            "name": "C: Verified procedural skill memory without repair",
            "enable_skill_memory": True,
            "enable_skill_verification": True,
            "enable_self_repair": False,
            "inject_synthesis_noise": True,
        },
        "D_verified_memory_with_repair": {
            "name": "D: Verified procedural skill memory with execution-feedback repair (Full PHARMAFORGE)",
            "enable_skill_memory": True,
            "enable_skill_verification": True,
            "enable_self_repair": True,
            "inject_synthesis_noise": True,
        },
    }

    all_results = {}

    for lead_key in ("CHEMBL289959", "CHEMBL290096"):
        lead_info = BENCHMARK_LEADS[lead_key]
        print(f"\n==================================================")
        print(f"STARTING OPTIMIZATION FOR {lead_key}: {lead_info['description']}")
        print(f"==================================================")
        all_results[lead_key] = {}

        for config_key, cfg in ablation_configs.items():
            print(f"\n--- Running Configuration {cfg['name']} on {lead_key} ---")
            
            # Isolated storage dir
            cache_dir = Path(f"skills_cache/{config_key}_{lead_key}")
            if cache_dir.exists():
                shutil.rmtree(cache_dir)
            cache_dir.mkdir(parents=True, exist_ok=True)
            
            # Initialize isolated library
            executor = IsolatedCodeExecutor()
            library = VerifiedProceduralSkillLibrary(
                storage_dir=cache_dir,
                oracle=oracle,
                executor=executor,
            )
            
            # Isolate initial skill library: preserve filter rules, remove pre-seeded generative skills
            # so the procedural skill lifecycle (synthesis, verification, repair, reuse) is exercised
            if not cfg["enable_skill_memory"]:
                library.skills.clear()
            else:
                library.skills = {k: v for k, v in library.skills.items() if "filter" in v.tags}

            orchestrator = PharmaForgeOrchestrator(
                oracle=oracle,
                skill_library=library,
                executor=executor,
                qsar_predictor=predictor,
                enable_skill_memory=cfg["enable_skill_memory"],
                enable_skill_verification=cfg["enable_skill_verification"],
                enable_self_repair=cfg["enable_self_repair"],
                inject_synthesis_noise=cfg.get("inject_synthesis_noise", False),
            )

            task = create_task(lead_key, max_iterations=5)
            t0 = time.perf_counter()
            res = orchestrator.run_optimization(task)
            elapsed = time.perf_counter() - t0

            # Attach run metadata
            res["total_duration_seconds"] = round(elapsed, 2)
            res["configuration_key"] = config_key
            res["configuration_name"] = cfg["name"]
            res["lead_id"] = lead_key
            res["initial_smiles"] = lead_info["smiles"]
            res["experimental_pchembl"] = lead_info["experimental_pchembl"]

            print(f"  Completed in {elapsed:.2f}s | Iterations: {res['iterations_executed']}")
            print(f"  Proposed: {res['total_candidates_proposed']}, Valid: {res['valid_candidates_count']} ({res['validity_rate_pct']}%)")
            print(f"  Skill Reuse: {res['skill_reuse_count']}, Admitted: {res['skills_admitted_count']}, Rejected: {res['skills_rejected_count']}")
            print(f"  Repair Attempts: {res['repair_attempts']}, Successful: {res['successful_repairs']}")
            print(f"  Pareto Front Size: {res['pareto_front_size']}, Hypervolume: {res['pareto_hypervolume']}")

            if res["pareto_front"]:
                top = res["pareto_front"][0]
                print(f"  Top Pareto Lead: {top['smiles']} | pIC50={top.get('affinity_score')} | QED={top.get('qed')} | LogP={top.get('logp')} | SAScore={top.get('sascore')}")

            all_results[lead_key][config_key] = res

    # Save full results JSON
    out_file = artifacts_dir / "benchmark_optimization_results.json"
    with open(out_file, "w", encoding="utf-8") as f:
        json.dump(all_results, f, indent=2)
    print(f"\nSaved all benchmark optimization results to {out_file}")

    return all_results

if __name__ == "__main__":
    run_experiment_suite()
