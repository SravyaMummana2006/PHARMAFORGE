"""Standardized scientific benchmark metrics: Pareto Hypervolume, validity rate, and skill reuse."""
from __future__ import annotations

import numpy as np
from typing import Any
from pharmaforge.core.types import MoleculeCandidate, OptimizationObjective, ObjectiveDirection


class BenchmarkEvaluator:
    """Computes rigorous quantitative metrics for agent comparison and ablation studies."""

    @staticmethod
    def calculate_chemical_validity_rate(total_proposed: int, valid_sanitized: int) -> float:
        """Percentage of generated chemical candidates passing RDKit sanitization and valence."""
        if total_proposed == 0:
            return 0.0
        return float(round((valid_sanitized / total_proposed) * 100.0, 2))

    @staticmethod
    def calculate_pareto_hypervolume(
        candidates: list[MoleculeCandidate],
        objectives: list[OptimizationObjective],
        reference_point: dict[str, float] | None = None,
        n_samples: int = 50000,
    ) -> float:
        """Compute the mathematical Hypervolume (HV) indicator of the non-dominated Pareto set.
        
        Evaluates the Lebesgue measure of the union of hypercubes dominated by the
        non-dominated Pareto front points in the normalized [0, 1]^d objective space.
        """
        if not candidates or not objectives:
            return 0.0

        ref = reference_point or {
            "qed": 0.0,
            "affinity": 4.0,
            "affinity_score": 4.0,
            "logp": 6.0,
            "sascore": 10.0,
        }

        # Extract normalized objective points in [0, 1]^d where 1.0 is optimal
        points = []
        for cand in candidates:
            pt = []
            for obj in objectives:
                val = cand.get_property(obj.name)
                if val is None:
                    val = ref.get(obj.name, 0.0)

                # Normalize each dimension to [0, 1] where 1.0 is optimal
                if obj.name.lower() == "qed":
                    # QED is natively bounded in [0, 1]; reference is 0.0
                    ref_val = ref.get("qed", 0.0)
                    norm_val = max(0.0, min(1.0, (val - ref_val) / (1.0 - ref_val + 1e-6)))
                elif obj.direction == ObjectiveDirection.MAXIMIZE:
                    ref_val = ref.get(obj.name, 0.0)
                    norm_val = max(0.0, min(1.0, (val - ref_val) / (10.0 - ref_val + 1e-6)))
                else:
                    ref_val = ref.get(obj.name, 10.0)
                    norm_val = max(0.0, min(1.0, (ref_val - val) / (ref_val + 1e-6)))
                pt.append(norm_val)
            points.append(pt)

        arr = np.array(points, dtype=np.float64)
        if len(arr) == 0:
            return 0.0

        d = len(objectives)
        # Filter points that strictly dominate the reference origin [0.0]^d
        valid_points = arr[np.all(arr > 0.0, axis=1)]
        if len(valid_points) == 0:
            return 0.0

        # Rigorous Monte Carlo integration of the dominated Lebesgue volume in [0, 1]^d
        rng = np.random.default_rng(42)
        samples = rng.uniform(0.0, 1.0, size=(n_samples, d))
        
        # Chunked evaluation to prevent high memory usage
        chunk_size = 5000
        dominated_count = 0
        for i in range(0, n_samples, chunk_size):
            chunk = samples[i : i + chunk_size]
            # sample is dominated if exists point p in valid_points s.t. p >= sample across all d
            is_dom = np.any(np.all(valid_points[None, :, :] >= chunk[:, None, :], axis=2), axis=1)
            dominated_count += int(np.sum(is_dom))

        volume = dominated_count / n_samples
        return float(round(volume, 4))

    @staticmethod
    def calculate_skill_reuse_ratio(
        total_executions: int,
        reused_existing_skills: int,
    ) -> float:
        """Percentage of procedural execution steps fulfilled by previously learned skills."""
        if total_executions == 0:
            return 0.0
        return float(round((reused_existing_skills / total_executions) * 100.0, 2))
