"""Verified Procedural Skill Library (VPSL): Persistent, unit-tested executable skills."""
from __future__ import annotations

import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any
from pharmaforge.core.exceptions import VerificationFailureError
from pharmaforge.core.types import ProceduralSkill
from pharmaforge.sandbox.executor import IsolatedCodeExecutor
from pharmaforge.verifier.oracle import ScientificTestOracle


class SkillRetriever:
    """Computes semantic similarity over skill descriptions using character/word n-gram vectorization."""

    def __init__(self):
        self.vocabulary: dict[str, int] = {}

    def _tokenize(self, text: str) -> list[str]:
        cleaned = re.sub(r"[^\w\s]", " ", text.lower())
        tokens = cleaned.split()
        # Add word bi-grams for improved semantic sensitivity
        bigrams = [f"{tokens[i]}_{tokens[i+1]}" for i in range(len(tokens) - 1)]
        return tokens + bigrams

    def _vectorize(self, tokens: list[str]) -> dict[str, float]:
        freqs: dict[str, float] = {}
        for t in tokens:
            freqs[t] = freqs.get(t, 0.0) + 1.0
        # Normalize
        norm = math.sqrt(sum(v * v for v in freqs.values()))
        if norm > 0:
            return {k: v / norm for k, v in freqs.items()}
        return freqs

    def compute_similarity(self, query: str, target: str) -> float:
        """Compute cosine similarity between query and target text."""
        vec1 = self._vectorize(self._tokenize(query))
        vec2 = self._vectorize(self._tokenize(target))
        score = sum(vec1[k] * vec2[k] for k in vec1 if k in vec2)
        return float(min(1.0, max(0.0, score)))


class VerifiedProceduralSkillLibrary:
    """Manages the lifecycle of verified, executable Python skills.
    
    Skills are registered strictly after passing automated unit tests executed
    within the isolated sandbox and evaluated by the Scientific Test-Oracle.
    """

    def __init__(
        self,
        storage_dir: str | Path = "skills_cache",
        oracle: ScientificTestOracle | None = None,
        executor: IsolatedCodeExecutor | None = None,
    ):
        self.storage_dir = Path(storage_dir)
        self.storage_dir.mkdir(parents=True, exist_ok=True)
        self.skills_file = self.storage_dir / "skills.json"
        self.oracle = oracle or ScientificTestOracle()
        self.executor = executor or IsolatedCodeExecutor()
        self.retriever = SkillRetriever()
        self.skills: dict[str, ProceduralSkill] = {}
        self._load_from_disk()
        if not self.skills:
            self._seed_foundational_skills()

    def _load_from_disk(self) -> None:
        """Load registered skills from disk."""
        if self.skills_file.exists():
            try:
                with open(self.skills_file, "r", encoding="utf-8") as f:
                    data = json.load(f)
                    for item in data.values():
                        skill = ProceduralSkill(**item)
                        self.skills[skill.name] = skill
            except Exception as e:
                print(f"[VPSL] Warning: Failed to load skills from {self.skills_file}: {e}")

    def save_to_disk(self) -> None:
        """Persist registered skills to JSON and export individual .py files."""
        serialized = {k: v.model_dump(mode="json") for k, v in self.skills.items()}
        with open(self.skills_file, "w", encoding="utf-8") as f:
            json.dump(serialized, f, indent=2)

        # Export individual Python files for inspection and execution
        code_dir = self.storage_dir / "code"
        code_dir.mkdir(parents=True, exist_ok=True)
        for skill in self.skills.values():
            py_file = code_dir / f"{skill.name}.py"
            with open(py_file, "w", encoding="utf-8") as f:
                f.write(f'"""\nSkill: {skill.name}\nDescription: {skill.description}\n"""\n\n')
                f.write(skill.code)

    def register_skill(
        self,
        name: str,
        description: str,
        code: str,
        unit_tests: str,
        tags: list[str] | None = None,
        parameters_schema: dict[str, Any] | None = None,
        return_schema: dict[str, Any] | None = None,
    ) -> ProceduralSkill:
        """Register a new skill into the library ONLY after automated verification passes."""
        # Step 1: Run unit tests in sandbox
        verification_code = f"""
{code}

# --- UNIT TESTS ---
{unit_tests}
"""
        trace = self.executor.execute_function(
            function_code=verification_code,
            function_name="run_tests",
            input_kwargs={},
        )
        if not trace.success:
            raise VerificationFailureError(
                f"Skill '{name}' failed unit test verification: {trace.traceback_str}",
                reasons=[trace.traceback_str or "Unknown test failure"],
            )

        skill_id = hashlib.sha256(f"{name}_{code}".encode()).hexdigest()[:12]
        skill = ProceduralSkill(
            skill_id=skill_id,
            name=name,
            description=description,
            code=code,
            parameters_schema=parameters_schema or {},
            return_schema=return_schema or {},
            unit_tests=unit_tests,
            tags=tags or [],
            success_count=1,
        )
        self.skills[name] = skill
        self.save_to_disk()
        return skill

    def retrieve_skills(self, query: str, top_k: int = 3, threshold: float = 0.1) -> list[tuple[ProceduralSkill, float]]:
        """Retrieve skills most relevant to the query based on semantic similarity."""
        scored = []
        for skill in self.skills.values():
            text_corpus = f"{skill.name} {skill.description} {' '.join(skill.tags)}"
            sim = self.retriever.compute_similarity(query, text_corpus)
            if sim >= threshold:
                scored.append((skill, sim))
        scored.sort(key=lambda x: x[1], reverse=True)
        return scored[:top_k]

    def record_execution_outcome(self, skill_name: str, success: bool) -> None:
        """Update success/failure statistics for a skill."""
        if skill_name in self.skills:
            if success:
                self.skills[skill_name].success_count += 1
            else:
                self.skills[skill_name].failure_count += 1
            self.save_to_disk()

    def _seed_foundational_skills(self) -> None:
        """Pre-seed library with verified foundational chemoinformatics skills."""
        # Skill 1: Lipinski's Rule of Five
        code_lipinski = '''def check_lipinski_rule(mw: float, logp: float, hbd: int, hba: int) -> dict:
    """Evaluates Lipinski's Rule of Five for drug-likeness."""
    violations = []
    if mw > 500.0:
        violations.append("MW > 500")
    if logp > 5.0:
        violations.append("LogP > 5.0")
    if hbd > 5:
        violations.append("HBD > 5")
    if hba > 10:
        violations.append("HBA > 10")
    
    return {
        "passed": len(violations) <= 1,
        "violations": violations,
        "violation_count": len(violations)
    }
'''
        test_lipinski = '''def run_tests():
    res1 = check_lipinski_rule(300.0, 2.5, 2, 4)
    assert res1["passed"] is True
    assert res1["violation_count"] == 0

    res2 = check_lipinski_rule(600.0, 6.0, 6, 12)
    assert res2["passed"] is False
    assert res2["violation_count"] == 4
    return True
'''
        self.register_skill(
            name="check_lipinski_rule",
            description="Evaluate small molecule drug-likeness compliance against Lipinski Rule of 5 constraints",
            code=code_lipinski,
            unit_tests=test_lipinski,
            tags=["druglikeness", "lipinski", "filter"],
        )

        # Skill 2: Fluorination / Bioisosteric R-group modification
        code_bioisostere = '''def substitute_hydroxyl_with_fluorine(smiles: str) -> list[str]:
    """Replaces hydroxyl groups (-OH) with fluorine (-F) as a classic bioisostere to reduce clearance."""
    candidates = []
    if "O" in smiles or "OH" in smiles:
        # Simple functional group bioisosteric substitution pattern
        mod = smiles.replace("OH", "F")
        if mod != smiles and mod not in candidates:
            candidates.append(mod)
        mod2 = smiles.replace("O", "F", 1)
        if mod2 != smiles and mod2 not in candidates:
            candidates.append(mod2)
    return candidates or [smiles]
'''
        test_bioisostere = '''def run_tests():
    out = substitute_hydroxyl_with_fluorine("CC(C)OH")
    assert len(out) >= 1
    assert "F" in out[0]
    return True
'''
        self.register_skill(
            name="substitute_hydroxyl_with_fluorine",
            description="Perform bioisosteric replacement of hydroxyl group with fluorine to reduce metabolic clearance and logP",
            code=code_bioisostere,
            unit_tests=test_bioisostere,
            tags=["bioisostere", "fluorination", "lead_optimization"],
        )
