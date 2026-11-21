"""Leakage-safe cross-validation splitters for molecular datasets.

Implements:
1. Random K-Fold CV (compound-level).
2. Bemis-Murcko Scaffold K-Fold CV (scaffold-disjoint).
3. Document Connected-Component K-Fold CV (compound- and document-disjoint).
"""

from __future__ import annotations
import random
from collections import defaultdict
from typing import Sequence
import numpy as np

def create_random_kfold_splits(
    n_samples: int, n_splits: int = 5, seed: int = 42
) -> list[tuple[np.ndarray, np.ndarray]]:
    """Generate reproducible Random K-Fold splits."""
    rng = np.random.default_rng(seed)
    indices = np.arange(n_samples)
    rng.shuffle(indices)
    folds = np.array_split(indices, n_splits)
    splits = []
    for fold_idx in range(n_splits):
        test_idx = folds[fold_idx]
        train_idx = np.concatenate([folds[i] for i in range(n_splits) if i != fold_idx])
        splits.append((train_idx, test_idx))
    return splits

def create_scaffold_kfold_splits(
    scaffolds: Sequence[str], n_splits: int = 5, seed: int = 42
) -> list[tuple[np.ndarray, np.ndarray]]:
    """Generate Bemis-Murcko scaffold-disjoint K-Fold splits via greedy balanced packing.
    
    Guarantees that all molecules with the same Murcko scaffold are strictly confined
    to the same fold (either train or test, never across).
    """
    scaffold_to_indices = defaultdict(list)
    for idx, scaf in enumerate(scaffolds):
        scaffold_to_indices[scaf].append(idx)

    # Sort scaffold clusters by size descending; shuffle ties reproducibly
    rng = random.Random(seed)
    clusters = list(scaffold_to_indices.values())
    rng.shuffle(clusters)
    clusters.sort(key=len, reverse=True)

    fold_indices: list[list[int]] = [[] for _ in range(n_splits)]
    fold_sizes: list[int] = [0 for _ in range(n_splits)]

    # Greedy bin packing into the smallest fold
    for cluster in clusters:
        min_fold = int(np.argmin(fold_sizes))
        fold_indices[min_fold].extend(cluster)
        fold_sizes[min_fold] += len(cluster)

    splits = []
    for fold_idx in range(n_splits):
        test_idx = np.array(fold_indices[fold_idx], dtype=np.int64)
        train_parts = [fold_indices[i] for i in range(n_splits) if i != fold_idx]
        train_idx = np.concatenate(train_parts) if train_parts else np.empty(0, dtype=np.int64)
        splits.append((train_idx, test_idx))
    return splits

def create_document_kfold_splits(
    document_ids_list: Sequence[Sequence[str]], n_splits: int = 5, seed: int = 42
) -> list[tuple[np.ndarray, np.ndarray]]:
    """Generate document- and compound-disjoint splits via bipartite connected components.
    
    Guarantees ZERO shared compounds AND ZERO shared publications across folds.
    """
    adj = defaultdict(set)
    n_samples = len(document_ids_list)
    for idx, docs in enumerate(document_ids_list):
        c_node = f"C_{idx}"
        for d in docs:
            d_node = f"D_{d}"
            adj[c_node].add(d_node)
            adj[d_node].add(c_node)

    visited = set()
    components = []
    for idx in range(n_samples):
        c_node = f"C_{idx}"
        if c_node in visited:
            continue
        comp_compounds = []
        queue = [c_node]
        visited.add(c_node)
        while queue:
            curr = queue.pop(0)
            if curr.startswith("C_"):
                comp_compounds.append(int(curr[2:]))
            for neighbor in adj[curr]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append(neighbor)
        components.append(comp_compounds)

    rng = random.Random(seed)
    rng.shuffle(components)
    components.sort(key=len, reverse=True)

    fold_indices: list[list[int]] = [[] for _ in range(n_splits)]
    fold_sizes: list[int] = [0 for _ in range(n_splits)]

    for comp in components:
        min_fold = int(np.argmin(fold_sizes))
        fold_indices[min_fold].extend(comp)
        fold_sizes[min_fold] += len(comp)

    splits = []
    for fold_idx in range(n_splits):
        test_idx = np.array(fold_indices[fold_idx], dtype=np.int64)
        train_parts = [fold_indices[i] for i in range(n_splits) if i != fold_idx]
        train_idx = np.concatenate(train_parts) if train_parts else np.empty(0, dtype=np.int64)
        splits.append((train_idx, test_idx))
    return splits
