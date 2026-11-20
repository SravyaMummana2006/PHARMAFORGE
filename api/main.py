"""FastAPI application for PHARMAFORGE scientific agent operations."""
from __future__ import annotations

from typing import Any
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field

from pharmaforge import __version__
from pharmaforge.agents.orchestrator import PharmaForgeOrchestrator
from pharmaforge.benchmarks.tdc_loader import get_standard_benchmark_suite
from pharmaforge.core.types import LeadOptimizationTask, OptimizationObjective, ObjectiveDirection
from pharmaforge.memory.skill_library import VerifiedProceduralSkillLibrary
from pharmaforge.verifier.oracle import ScientificTestOracle

app = FastAPI(
    title="PHARMAFORGE Scientific Agent API",
    description="A Self-Improving Multi-Agent Framework for Scientific Reasoning and Pharmaceutical Discovery",
    version=__version__,
)

# Instantiate core singletons
oracle = ScientificTestOracle()
skill_library = VerifiedProceduralSkillLibrary(oracle=oracle)
orchestrator = PharmaForgeOrchestrator(oracle=oracle, skill_library=skill_library)


class VerifyRequest(BaseModel):
    smiles: str = Field(..., example="CC(=O)OC1=CC=CC=C1C(=O)O")


class SkillQueryRequest(BaseModel):
    query: str = Field(..., example="substitute hydroxyl with fluorine to improve logP")
    top_k: int = Field(default=3, ge=1, le=10)


class OptimizeRequest(BaseModel):
    task_id: str = Field(default="custom_task_01")
    target_name: str = Field(default="Target_X")
    initial_smiles: str = Field(..., example="c1ccccc1O")
    objectives: list[dict[str, Any]] = Field(
        default=[
            {"name": "qed", "direction": "maximize", "threshold": 0.5},
            {"name": "logp", "direction": "minimize", "threshold": 3.0},
        ]
    )
    max_iterations: int = Field(default=5, ge=1, le=20)


@app.get("/")
def root() -> dict[str, str]:
    return {
        "framework": "PHARMAFORGE",
        "version": __version__,
        "status": "operational",
        "oracle_backend": "RDKit" if oracle.is_rdkit_backed else "Heuristic-Deterministic",
    }


@app.get("/health")
def health_check() -> dict[str, Any]:
    return {
        "status": "healthy",
        "registered_skills": len(skill_library.skills),
        "rdkit_available": oracle.is_rdkit_backed,
    }


@app.post("/api/v1/verify")
def verify_molecule(payload: VerifyRequest) -> dict[str, Any]:
    """Verify chemical validity and compute physicochemical properties via STO-IV."""
    report = oracle.verify_smiles(payload.smiles)
    return report.model_dump()


@app.get("/api/v1/skills")
def list_skills() -> list[dict[str, Any]]:
    """List all verified procedural skills currently in the library."""
    return [s.model_dump() for s in skill_library.skills.values()]


@app.post("/api/v1/skills/query")
def query_skills(payload: SkillQueryRequest) -> list[dict[str, Any]]:
    """Retrieve procedural skills ranked by semantic relevance."""
    results = skill_library.retrieve_skills(payload.query, top_k=payload.top_k)
    return [{"skill": s.model_dump(), "similarity": score} for s, score in results]


@app.get("/api/v1/benchmarks")
def list_benchmarks() -> list[dict[str, Any]]:
    """Return standardized TDC benchmark tasks."""
    suite = get_standard_benchmark_suite()
    return [task.model_dump() for task in suite]


@app.post("/api/v1/optimize")
def run_optimization(payload: OptimizeRequest) -> dict[str, Any]:
    """Execute closed-loop multi-objective lead optimization."""
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

    result = orchestrator.run_optimization(task)
    return result
