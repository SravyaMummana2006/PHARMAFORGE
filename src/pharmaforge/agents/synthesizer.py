"""Tool & Skill Synthesizer: Generates verified procedural Python functions using RDKit."""
from __future__ import annotations

from typing import Any
from pharmaforge.core.types import ProceduralSkill
from pharmaforge.memory.skill_library import VerifiedProceduralSkillLibrary


class SkillSynthesizer:
    """Synthesizes executable, typed Python procedures using RDKit chemical operations."""

    def __init__(self, skill_library: VerifiedProceduralSkillLibrary):
        self.library = skill_library

    def synthesize_skill_components(
        self,
        target_property: str,
        transformation_type: str = "bioisostere",
        inject_synthesis_noise: bool = False,
    ) -> tuple[str, str, str, str, list[str]]:
        """Return raw (name, description, code, unit_tests, tags) before registration."""
        prop = target_property.lower()

        if prop in ("affinity", "affinity_score"):
            name = "rdkit_aromatic_affinity_optimization"
            desc = "Introduces halogen, methyl, or methoxy substituents onto aromatic rings using RDKit Reaction SMARTS to optimize target binding."
            code = '''def rdkit_aromatic_affinity_optimization(smiles: str) -> list[str]:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []
    
    # Established kinase hinge and pocket optimization reactions
    smarts_patterns = [
        "[c;H1:1]>>[c:1]Br",   # Gatekeeper bromine addition (e.g. CHEMBL290096)
        "[c;H1:1]>>[c:1]Cl",   # Lipophilic pocket chlorine addition
        "[c;H1:1]>>[c:1]F",    # Metabolic stabilization fluorine addition
        "[c;H1:1]>>[c:1]OC",   # Electron-donating methoxy addition (e.g. gefitinib 6,7-OMe)
        "[c;H1:1]>>[c:1]C",    # Methyl pocket contact
    ]
    
    candidates = []
    for smarts in smarts_patterns:
        try:
            rxn = AllChem.ReactionFromSmarts(smarts)
            products = rxn.RunReactants((mol,))
            for prod_tuple in products:
                for prod in prod_tuple:
                    try:
                        Chem.SanitizeMol(prod)
                        s = Chem.MolToSmiles(prod, canonical=True)
                        if s and s != smiles and s not in candidates:
                            candidates.append(s)
                    except Exception:
                        pass
        except Exception:
            continue
            
    return candidates[:10]
'''
            unit_tests = '''def run_tests():
    from rdkit import Chem
    res = rdkit_aromatic_affinity_optimization("c1ccccc1")
    assert len(res) >= 1
    # Verify every generated candidate is valid RDKit molecule
    for s in res:
        m = Chem.MolFromSmiles(s)
        assert m is not None
    return True
'''
            tags = ["affinity", "kinase", "halogenation", "rdkit_rxn"]

        elif prop in ("solubility", "logp"):
            name = "rdkit_polar_solubilizing_modification"
            desc = "Introduces polar nitrogen, hydroxyl, or heterocyclic groups using RDKit to reduce logP and improve aqueous solubility."
            code = '''def rdkit_polar_solubilizing_modification(smiles: str) -> list[str]:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []

    smarts_patterns = [
        "[c;H1:1]>>[n:1]",     # Aromatic CH to N bioisostere (pyridine insertion)
        "[c;H1:1]>>[c:1]O",    # Phenolic hydroxyl
        "[c;H1:1]>>[c:1]N",    # Primary amino substitution
        "[c;H1:1]>>[c:1]C#N",  # Polar cyano dipole
    ]

    candidates = []
    for smarts in smarts_patterns:
        try:
            rxn = AllChem.ReactionFromSmarts(smarts)
            products = rxn.RunReactants((mol,))
            for prod_tuple in products:
                for prod in prod_tuple:
                    try:
                        Chem.SanitizeMol(prod)
                        s = Chem.MolToSmiles(prod, canonical=True)
                        if s and s != smiles and s not in candidates:
                            candidates.append(s)
                    except Exception:
                        pass
        except Exception:
            continue

    return candidates[:10]
'''
            unit_tests = '''def run_tests():
    from rdkit import Chem
    res = rdkit_polar_solubilizing_modification("c1ccccc1")
    assert len(res) >= 1
    for s in res:
        m = Chem.MolFromSmiles(s)
        assert m is not None
    return True
'''
            tags = ["solubility", "logp", "polar_group", "rdkit_rxn"]

        elif prop in ("qed", "druglikeness"):
            name = "rdkit_qed_bioisostere_optimization"
            desc = "Balances lipophilicity, polar surface area, and molecular weight via RDKit bioisosteric functionalizations to maximize QED."
            if inject_synthesis_noise:
                # Simulates realistic LLM synthesis bug: missing Chem import
                code = '''def rdkit_qed_bioisostere_optimization(smiles: str) -> list[str]:
    from rdkit.Chem import AllChem
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []

    smarts_patterns = [
        "[c;H1:1]>>[c:1]F",
        "[c;H1:1]>>[c:1]OC",
    ]

    candidates = []
    for smarts in smarts_patterns:
        try:
            rxn = AllChem.ReactionFromSmarts(smarts)
            products = rxn.RunReactants((mol,))
            for prod_tuple in products:
                for prod in prod_tuple:
                    try:
                        Chem.SanitizeMol(prod)
                        s = Chem.MolToSmiles(prod, canonical=True)
                        if s and s != smiles and s not in candidates:
                            candidates.append(s)
                    except Exception:
                        pass
        except Exception:
            continue

    return candidates[:10]
'''
            else:
                code = '''def rdkit_qed_bioisostere_optimization(smiles: str) -> list[str]:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []

    smarts_patterns = [
        "[c;H1:1]>>[c:1]F",    # Fluorine substitution (improves QED)
        "[c;H1:1]>>[c:1]OC",   # Methoxy substitution
        "[c;H1:1]>>[n:1]",     # Heteroatom insertion
    ]

    candidates = []
    for smarts in smarts_patterns:
        try:
            rxn = AllChem.ReactionFromSmarts(smarts)
            products = rxn.RunReactants((mol,))
            for prod_tuple in products:
                for prod in prod_tuple:
                    try:
                        Chem.SanitizeMol(prod)
                        s = Chem.MolToSmiles(prod, canonical=True)
                        if s and s != smiles and s not in candidates:
                            candidates.append(s)
                    except Exception:
                        pass
        except Exception:
            continue

    return candidates[:10]
'''
            unit_tests = '''def run_tests():
    from rdkit import Chem
    res = rdkit_qed_bioisostere_optimization("c1ccccc1")
    assert len(res) >= 1
    for s in res:
        m = Chem.MolFromSmiles(s)
        assert m is not None
    return True
'''
            tags = ["qed", "druglikeness", "bioisostere", "rdkit_rxn"]

        else:
            name = f"rdkit_substituent_variation_{prop}"
            desc = f"Algorithmic structural variation targeting {prop} optimization via RDKit Reaction SMARTS"
            code = f'''def {name}(smiles: str) -> list[str]:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []

    smarts_patterns = [
        "[c;H1:1]>>[c:1]C",
        "[c;H1:1]>>[c:1]F",
        "[c;H1:1]>>[c:1]Cl",
    ]

    candidates = []
    for smarts in smarts_patterns:
        try:
            rxn = AllChem.ReactionFromSmarts(smarts)
            products = rxn.RunReactants((mol,))
            for prod_tuple in products:
                for prod in prod_tuple:
                    try:
                        Chem.SanitizeMol(prod)
                        s = Chem.MolToSmiles(prod, canonical=True)
                        if s and s != smiles and s not in candidates:
                            candidates.append(s)
                    except Exception:
                        pass
        except Exception:
            continue

    return candidates[:10]
'''
            unit_tests = f'''def run_tests():
    from rdkit import Chem
    res = {name}("c1ccccc1")
    assert len(res) >= 1
    for s in res:
        m = Chem.MolFromSmiles(s)
        assert m is not None
    return True
'''
            tags = [prop, "optimization", "rdkit_rxn"]

        if inject_synthesis_noise:
            # Simulate realistic LLM generation bug: missing Chem import
            code = code.replace("    from rdkit import Chem\n", "")

        return name, desc, code, unit_tests, tags

    def synthesize_transformation_skill(
        self,
        target_property: str,
        transformation_type: str = "bioisostere",
        bypass_verification: bool = False,
        inject_synthesis_noise: bool = False,
    ) -> ProceduralSkill:
        """Synthesize and register a procedural transformation skill tailored for a property objective."""
        name, desc, code, unit_tests, tags = self.synthesize_skill_components(
            target_property=target_property,
            transformation_type=transformation_type,
            inject_synthesis_noise=inject_synthesis_noise,
        )
        return self.library.register_skill(
            name=name,
            description=desc,
            code=code,
            unit_tests=unit_tests,
            tags=tags,
            bypass_verification=bypass_verification,
        )
