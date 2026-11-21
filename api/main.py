"""FastAPI application for PHARMAFORGE scientific agent operations."""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field

# Ensure src directory is in sys.path
root_dir = Path(__file__).resolve().parent.parent
src_dir = root_dir / "src"
if str(src_dir) not in sys.path:
    sys.path.insert(0, str(src_dir))

from pharmaforge import __version__
from pharmaforge.core.types import LeadOptimizationTask, OptimizationObjective, ObjectiveDirection

def get_rdkit_python() -> str:
    """Find a Python interpreter capable of executing RDKit and scientific tasks."""
    try:
        import rdkit
        return sys.executable
    except ImportError:
        pass
    env_py = os.environ.get("PHARMAFORGE_PYTHON")
    if env_py and Path(env_py).exists():
        return env_py
    local_py313 = Path.home() / "AppData" / "Local" / "Programs" / "Python" / "Python313" / "python.exe"
    if local_py313.exists():
        return str(local_py313)
    return "python"

RDKIT_PYTHON = get_rdkit_python()

app = FastAPI(
    title="PHARMAFORGE Scientific Agent API",
    description="A Self-Improving Multi-Agent Framework for Scientific Reasoning and Pharmaceutical Discovery",
    version=__version__,
)

# Instantiate core singletons safely if in-process RDKit is available
oracle = None
skill_library = None
qsar_predictor = None
orchestrator = None

try:
    from pharmaforge.verifier.oracle import ScientificTestOracle
    from pharmaforge.memory.skill_library import VerifiedProceduralSkillLibrary
    from pharmaforge.agents.orchestrator import PharmaForgeOrchestrator
    from pharmaforge.qsar.predictor import EGFRQSARPredictor
    
    oracle = ScientificTestOracle()
    if oracle.is_rdkit_backed:
        skill_library = VerifiedProceduralSkillLibrary(oracle=oracle)
        try:
            qsar_predictor = EGFRQSARPredictor()
            if not qsar_predictor.is_ready:
                qsar_predictor = None
        except Exception:
            qsar_predictor = None

        orchestrator = PharmaForgeOrchestrator(
            oracle=oracle,
            skill_library=skill_library,
            qsar_predictor=qsar_predictor,
        )
except Exception:
    oracle = None
    skill_library = None
    qsar_predictor = None
    orchestrator = None


class VerifyRequest(BaseModel):
    smiles: str = Field(..., example="CC(=O)OC1=CC=CC=C1C(=O)O")


class SkillQueryRequest(BaseModel):
    query: str = Field(..., example="substitute hydroxyl with fluorine to improve logP")
    top_k: int = Field(default=3, ge=1, le=10)


class OptimizeRequest(BaseModel):
    task_id: str = Field(default="EGFR-BENCH-CHEMBL289959")
    target_name: str = Field(default="EGFR")
    initial_smiles: str = Field(default="c1ccc(Nc2ncnc3ccccc23)cc1")
    objectives: list[dict[str, Any]] = Field(
        default=[
            {"name": "affinity", "direction": "maximize", "threshold": 7.0, "weight": 2.0},
            {"name": "qed", "direction": "maximize", "threshold": 0.5, "weight": 1.0},
            {"name": "logp", "direction": "minimize", "threshold": 3.5, "weight": 1.0},
            {"name": "sascore", "direction": "minimize", "threshold": 4.5, "weight": 1.0},
        ]
    )
    max_iterations: int = Field(default=1, ge=1, le=10)


@app.get("/")
def root() -> dict[str, str]:
    return {
        "framework": "PHARMAFORGE",
        "version": __version__,
        "status": "operational",
        "oracle_backend": "InProcess-RDKit" if (oracle and oracle.is_rdkit_backed) else "Subprocess-RDKit",
    }


@app.get("/health")
def health_check() -> dict[str, Any]:
    skills_count = len(skill_library.skills) if skill_library else 0
    rdkit_ready = bool(oracle and getattr(oracle, "is_rdkit_backed", False))
    qsar_ready = bool(qsar_predictor and getattr(qsar_predictor, "is_ready", False))
    return {
        "status": "healthy",
        "registered_skills": skills_count,
        "rdkit_available": rdkit_ready,
        "qsar_model_loaded": qsar_ready,
        "version": __version__,
    }


@app.post("/api/v1/verify")
def verify_molecule(payload: VerifyRequest) -> dict[str, Any]:
    """Verify chemical validity and compute physicochemical properties via STO-IV."""
    if oracle and oracle.is_rdkit_backed:
        report = oracle.verify_smiles(payload.smiles)
        return report.model_dump()
    
    script = (
        'import sys, json; sys.path.insert(0, "src"); '
        'from pharmaforge.verifier.oracle import ScientificTestOracle; '
        'o = ScientificTestOracle(); '
        f'r = o.verify_smiles({json.dumps(payload.smiles)}); '
        'print(json.dumps(r.model_dump()))'
    )
    res = subprocess.run([RDKIT_PYTHON, "-c", script], capture_output=True, text=True, cwd=str(root_dir))
    if res.returncode == 0:
        return json.loads(res.stdout.strip())
    raise HTTPException(status_code=500, detail=f"Verification failed: {res.stderr}")


@app.get("/api/v1/skills")
def list_skills() -> list[dict[str, Any]]:
    """List all verified procedural skills currently in the library."""
    if skill_library:
        return [s.model_dump() for s in skill_library.skills.values()]
    
    skills_file = root_dir / "skills_cache" / "skills.json"
    if skills_file.exists():
        with open(skills_file, "r") as f:
            data = json.load(f)
            return list(data.values())
    return []


@app.post("/api/v1/skills/query")
def query_skills(payload: SkillQueryRequest) -> list[dict[str, Any]]:
    """Retrieve procedural skills ranked by semantic relevance."""
    if skill_library:
        results = skill_library.retrieve_skills(payload.query, top_k=payload.top_k)
        return [{"skill": s.model_dump(), "similarity": score} for s, score in results]
    
    skills_file = root_dir / "skills_cache" / "skills.json"
    if skills_file.exists():
        with open(skills_file, "r") as f:
            data = json.load(f)
            results = []
            q_words = payload.query.lower().split()
            for s in data.values():
                score = 0.5 + 0.1 * sum(1 for w in q_words if w in s.get("name", "").lower() or w in s.get("description", "").lower())
                results.append({"skill": s, "similarity": min(round(score, 2), 0.95)})
            results.sort(key=lambda x: x["similarity"], reverse=True)
            return results[:payload.top_k]
    return []


@app.get("/api/v1/benchmarks")
def list_benchmarks() -> list[dict[str, Any]]:
    """Return standardized TDC benchmark tasks."""
    try:
        from pharmaforge.benchmarks.tdc_loader import get_standard_benchmark_suite
        suite = get_standard_benchmark_suite()
        return [task.model_dump() for task in suite]
    except Exception:
        bench_file = root_dir / "artifacts" / "benchmark_runs" / "benchmark_optimization_results.json"
        if bench_file.exists():
            with open(bench_file, "r") as f:
                d = json.load(f)
                return [{"task_id": k, "lead_id": k} for k in d.keys()]
        return []


@app.post("/api/v1/optimize")
def run_optimization(payload: OptimizeRequest) -> dict[str, Any]:
    """Execute closed-loop multi-objective lead optimization."""
    if orchestrator:
        try:
            from pharmaforge.benchmarks.tdc_loader import get_standard_benchmark_suite
            tasks = get_standard_benchmark_suite()
            target_task = next((t for t in tasks if t.task_id == payload.task_id), None)
        except Exception:
            target_task = None

        if target_task:
            target_task.max_iterations = payload.max_iterations
            task = target_task
        else:
            objectives = [
                OptimizationObjective(
                    name=obj["name"],
                    direction=ObjectiveDirection(obj.get("direction", "maximize")),
                    threshold=obj.get("threshold"),
                    weight=obj.get("weight", 1.0),
                )
                for obj in payload.objectives
            ]
            task = LeadOptimizationTask(
                task_id=payload.task_id,
                target_name=payload.target_name,
                initial_smiles=payload.initial_smiles,
                objectives=objectives,
                max_iterations=payload.max_iterations,
            )
        return orchestrator.run_optimization(task)
    
    # Subprocess execution bridge
    runner_script = root_dir / "scripts" / "run_optimization_task.py"
    cmd = [
        RDKIT_PYTHON,
        str(runner_script),
        "--task_id", payload.task_id,
        "--iterations", str(payload.max_iterations)
    ]
    res = subprocess.run(cmd, capture_output=True, text=True, cwd=str(root_dir))
    if res.returncode == 0:
        lines = [line.strip() for line in res.stdout.strip().split("\n") if line.strip()]
        for line in reversed(lines):
            if line.startswith("{") and line.endswith("}"):
                try:
                    return json.loads(line)
                except Exception:
                    pass
    raise HTTPException(status_code=500, detail=f"Optimization runner failed: {res.stderr or res.stdout}")
