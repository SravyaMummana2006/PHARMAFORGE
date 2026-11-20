"""Benchmark task loader: Therapeutics Data Commons (TDC) lead optimization suites."""
from __future__ import annotations

from pharmaforge.core.types import LeadOptimizationTask, OptimizationObjective, ObjectiveDirection


def get_standard_benchmark_suite() -> list[LeadOptimizationTask]:
    """Returns a standardized suite of 4 pharmaceutical lead optimization benchmark tasks."""
    return [
        # Task 1: EGFR Kinase Lead Optimization (Potency + Drug-likeness + SAScore)
        LeadOptimizationTask(
            task_id="TDC-EGFR-01",
            target_name="EGFR",
            target_class="Kinase",
            initial_smiles="COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OC",
            objectives=[
                OptimizationObjective(name="affinity_score", direction=ObjectiveDirection.MAXIMIZE, threshold=7.0, weight=2.0),
                OptimizationObjective(name="qed", direction=ObjectiveDirection.MAXIMIZE, threshold=0.55, weight=1.0),
                OptimizationObjective(name="logp", direction=ObjectiveDirection.MINIMIZE, threshold=3.8, weight=1.0),
                OptimizationObjective(name="sascore", direction=ObjectiveDirection.MINIMIZE, threshold=4.0, weight=1.0),
            ],
            max_iterations=8,
        ),
        # Task 2: Dopamine D2 Receptor (GPCR) Multi-Objective Balancing
        LeadOptimizationTask(
            task_id="TDC-DRD2-02",
            target_name="DRD2",
            target_class="GPCR",
            initial_smiles="O=C(c1ccc(F)cc1)CCCN2CCC(O)(c3ccc(Cl)cc3)CC2",
            objectives=[
                OptimizationObjective(name="affinity_score", direction=ObjectiveDirection.MAXIMIZE, threshold=7.2, weight=2.0),
                OptimizationObjective(name="qed", direction=ObjectiveDirection.MAXIMIZE, threshold=0.60, weight=1.0),
                OptimizationObjective(name="logp", direction=ObjectiveDirection.MINIMIZE, threshold=3.5, weight=1.5),
            ],
            max_iterations=8,
        ),
        # Task 3: Solubility Optimization for Poorly Soluble Scaffold
        LeadOptimizationTask(
            task_id="TDC-AqSol-03",
            target_name="Aqueous_Solubility",
            target_class="ADMET",
            initial_smiles="c1ccc2c(c1)c3ccccc3c4ccccc24",
            objectives=[
                OptimizationObjective(name="logp", direction=ObjectiveDirection.MINIMIZE, threshold=2.5, weight=2.0),
                OptimizationObjective(name="qed", direction=ObjectiveDirection.MAXIMIZE, threshold=0.60, weight=1.0),
                OptimizationObjective(name="sascore", direction=ObjectiveDirection.MINIMIZE, threshold=3.5, weight=1.0),
            ],
            max_iterations=6,
        ),
        # Task 4: Cardiotoxicity (hERG) De-risking
        LeadOptimizationTask(
            task_id="TDC-hERG-04",
            target_name="hERG_DeRisking",
            target_class="Safety",
            initial_smiles="c1ccc(CCN2CCN(CC2)c3ccccc3)cc1",
            objectives=[
                OptimizationObjective(name="affinity_score", direction=ObjectiveDirection.MAXIMIZE, threshold=6.5, weight=1.5),
                OptimizationObjective(name="logp", direction=ObjectiveDirection.MINIMIZE, threshold=3.0, weight=1.0),
                OptimizationObjective(name="qed", direction=ObjectiveDirection.MAXIMIZE, threshold=0.65, weight=1.0),
            ],
            max_iterations=6,
        ),
    ]
