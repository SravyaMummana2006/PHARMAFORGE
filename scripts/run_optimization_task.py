"""CLI script to run live closed-loop lead optimization for a single task."""
import sys
import json
from pathlib import Path
import argparse

root_dir = Path(__file__).resolve().parent.parent
src_dir = root_dir / "src"
if str(src_dir) not in sys.path:
    sys.path.insert(0, str(src_dir))
if str(root_dir) not in sys.path:
    sys.path.insert(0, str(root_dir))

from pharmaforge.verifier.oracle import ScientificTestOracle
from pharmaforge.sandbox.executor import IsolatedCodeExecutor
from pharmaforge.memory.skill_library import VerifiedProceduralSkillLibrary
from pharmaforge.qsar.predictor import EGFRQSARPredictor
from pharmaforge.agents.orchestrator import PharmaForgeOrchestrator
from pharmaforge.benchmarks.tdc_loader import get_standard_benchmark_suite
from scripts.run_phase_1_4_experiments import create_task

def main():
    parser = argparse.ArgumentParser(description="Run single optimization task.")
    parser.add_argument("--task-id", type=str, required=True, help="Task ID (e.g. EGFR-BENCH-CHEMBL289959)")
    parser.add_argument("--iterations", type=int, default=5, help="Max iterations")
    parser.add_argument("--output-json", type=str, default=None, help="Output JSON path")
    args = parser.parse_args()

    # Find task
    suite = get_standard_benchmark_suite()
    task = None
    for t in suite:
        if t.task_id == args.task_id:
            task = t
            task.max_iterations = args.iterations
            break
    
    if task is None:
        if "CHEMBL289959" in args.task_id:
            task = create_task("CHEMBL289959", max_iterations=args.iterations)
        elif "CHEMBL290096" in args.task_id:
            task = create_task("CHEMBL290096", max_iterations=args.iterations)
        else:
            print(json.dumps({"success": False, "error": f"Unknown task ID: {args.task_id}"}))
            sys.exit(1)

    oracle = ScientificTestOracle()
    executor = IsolatedCodeExecutor()
    library = VerifiedProceduralSkillLibrary(storage_dir="skills_cache/live_run", oracle=oracle, executor=executor)
    predictor = EGFRQSARPredictor()

    orchestrator = PharmaForgeOrchestrator(
        oracle=oracle,
        skill_library=library,
        executor=executor,
        qsar_predictor=predictor if predictor.is_ready else None,
        enable_skill_memory=True,
        enable_skill_verification=True,
        enable_self_repair=True,
        inject_synthesis_noise=False,
    )

    result = orchestrator.run_optimization(task)
    serialized = json.dumps(result, indent=2)
    
    if args.output_json:
        with open(args.output_json, "w", encoding="utf-8") as f:
            f.write(serialized)
            
    print(serialized)

if __name__ == "__main__":
    main()
