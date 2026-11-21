#!/usr/bin/env python3
"""Scientific Environment & RDKit Diagnostic Script for PHARMAFORGE.

Validates the presence, version, and correct execution of all mandatory
cheminformatics dependencies (RDKit, descriptors, QED, FilterCatalog,
AllChem, Morgan fingerprints, and Reaction SMARTS).

Returns:
    Exit code 0 on complete success.
    Exit code 1 on any dependency failure or mathematical inconsistency.
"""
from __future__ import annotations

import sys


def run_diagnostics() -> int:
    print("=" * 70)
    print("PHARMAFORGE Scientific Environment & RDKit Diagnostics")
    print("=" * 70)

    # 1. Python Environment Check
    py_ver = sys.version.split()[0]
    print(f"[1/7] Python Runtime: {py_ver}")
    if sys.version_info < (3, 11):
        print(f"FAIL: Python >= 3.11 is required, but found {py_ver}")
        return 1

    # 2. RDKit Core & Module Availability Check
    print("[2/7] Checking Mandatory RDKit Modules...")
    try:
        import rdkit
        from rdkit import Chem
        from rdkit.Chem import Descriptors, QED, FilterCatalog, AllChem
        print(f"      - rdkit: version {rdkit.__version__}")
        print("      - rdkit.Chem: OK")
        print("      - rdkit.Chem.Descriptors: OK")
        print("      - rdkit.Chem.QED: OK")
        print("      - rdkit.Chem.FilterCatalog: OK")
        print("      - rdkit.Chem.AllChem: OK")
    except ImportError as e:
        print(f"FAIL: Mandatory RDKit module missing: {e}")
        return 1

    # Check Ertl SAScore availability
    try:
        from rdkit.Contrib.SA_Score import sascorer
        print("      - rdkit.Contrib.SA_Score.sascorer: OK (Fragment database loaded)")
    except Exception as e:
        print(f"      - rdkit.Contrib.SA_Score.sascorer: UNAVAILABLE ({e})")

    # 3. Molecular Parsing & Valence Validation Check
    print("[3/7] Validating Molecular Parsing & Rejection Oracles...")
    aspirin_smiles = "CC(=O)Oc1ccccc1C(=O)O"
    mol_valid = Chem.MolFromSmiles(aspirin_smiles)
    if mol_valid is None:
        print("FAIL: Failed to parse valid Aspirin SMILES.")
        return 1
    Chem.SanitizeMol(mol_valid)
    print("      - Valid SMILES parsed and sanitized successfully (Aspirin).")

    # Syntax error rejection test
    bad_syntax = "CC(C(=O)O"
    if Chem.MolFromSmiles(bad_syntax) is not None:
        print("FAIL: Oracle accepted malformed SMILES with unbalanced parentheses.")
        return 1
    print("      - Syntax error rejected as expected.")

    # Valence error rejection test (pentavalent nitrogen)
    pentavalent_n = "CN(C)(C)(C)C"
    mol_bad_valence = Chem.MolFromSmiles(pentavalent_n)
    valence_caught = False
    if mol_bad_valence is not None:
        try:
            Chem.SanitizeMol(mol_bad_valence)
        except Exception:
            valence_caught = True
    else:
        valence_caught = True

    if not valence_caught:
        print("FAIL: Oracle failed to catch invalid explicit valence on pentavalent nitrogen.")
        return 1
    print("      - Invalid valence (pentavalent nitrogen) rejected by SanitizeMol.")

    # 4. Physicochemical Descriptor Calculations
    print("[4/7] Checking Exact Physicochemical Descriptors...")
    mw = Descriptors.MolWt(mol_valid)
    logp = Descriptors.MolLogP(mol_valid)
    tpsa = Descriptors.TPSA(mol_valid)
    qed_score = QED.qed(mol_valid)

    # Reference values for Aspirin (C9H8O4, MW ~180.157, LogP ~1.31)
    if not (180.0 < mw < 180.3):
        print(f"FAIL: Molecular weight calculation incorrect ({mw} != 180.16).")
        return 1
    if not (1.2 < logp < 1.4):
        print(f"FAIL: Wildman-Crippen LogP calculation incorrect ({logp} != 1.31).")
        return 1
    if not (60.0 < tpsa < 65.0):
        print(f"FAIL: TPSA calculation incorrect ({tpsa} != 63.6).")
        return 1
    if not (0.5 < qed_score < 0.6):
        print(f"FAIL: Bickerton QED score incorrect ({qed_score} != 0.54).")
        return 1

    print(f"      - Exact MolWt: {mw:.4f} Da")
    print(f"      - MolLogP: {logp:.4f}")
    print(f"      - TPSA: {tpsa:.4f} Å²")
    print(f"      - Bickerton QED: {qed_score:.4f}")

    # Check Ertl SAScore
    try:
        from rdkit.Contrib.SA_Score import sascorer
        sa = sascorer.calculateScore(mol_valid)
        print(f"      - Ertl SAScore: {sa:.4f}")
    except Exception:
        print("      - Ertl SAScore: Skipped (Contrib not installed)")

    # 5. FilterCatalog Substructure Checking (PAINS & Brenk)
    print("[5/7] Testing RDKit FilterCatalog (PAINS & Brenk)...")
    pains_params = FilterCatalog.FilterCatalogParams()
    pains_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A)
    pains_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B)
    pains_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C)
    pains_catalog = FilterCatalog.FilterCatalog(pains_params)

    # Test Ene-rhodanine derivative (classic PAINS alert: ene_rhod_A)
    rhodanine = Chem.MolFromSmiles("O=C1NC(=S)SC1=Cc1ccccc1")
    matches = [entry.GetDescription() for entry in pains_catalog.GetMatches(rhodanine)]
    if not matches:
        print("FAIL: PAINS filter failed to detect alert in ene-rhodanine test molecule.")
        return 1
    print(f"      - PAINS alert successfully detected on Ene-rhodanine: {matches}")

    # Test Aspirin on PAINS (must be clean)
    aspirin_pains = [entry.GetDescription() for entry in pains_catalog.GetMatches(mol_valid)]
    if aspirin_pains:
        print(f"FAIL: False positive PAINS alert on Aspirin: {aspirin_pains}")
        return 1
    print("      - Non-alert control (Aspirin) clean of PAINS matches.")

    # 6. Morgan Fingerprint Generation Check
    print("[6/7] Testing Morgan Fingerprint (ECFP4) Generation...")
    try:
        from rdkit.Chem import rdFingerprintGenerator
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        fp = gen.GetFingerprint(mol_valid)
    except (ImportError, AttributeError):
        fp = AllChem.GetMorganFingerprintAsBitVect(mol_valid, 2, nBits=2048)

    num_on_bits = fp.GetNumOnBits()
    if num_on_bits <= 0 or num_on_bits > 2048:
        print(f"FAIL: Invalid bit count in Morgan fingerprint: {num_on_bits}")
        return 1
    print(f"      - Morgan Fingerprint (radius=2, 2048 bits) generated: {num_on_bits} active bits.")

    # 7. Reaction SMARTS Execution & Verification Check
    print("[7/7] Validating Chemical Reaction SMARTS Execution (Amide Coupling)...")
    # Benzoic acid + Methylamine -> N-methylbenzamide
    rxn_smarts = "[C:1](=[O:2])[OH].[N:3][C:4]>>[C:1](=[O:2])[N:3][C:4]"
    rxn = AllChem.ReactionFromSmarts(rxn_smarts)
    if rxn is None:
        print("FAIL: Failed to initialize ReactionFromSmarts.")
        return 1
    rxn.Initialize()

    acid = Chem.MolFromSmiles("c1ccccc1C(=O)O")
    amine = Chem.MolFromSmiles("NC")
    if not (rxn.IsMoleculeReactant(acid) and rxn.IsMoleculeReactant(amine)):
        print("FAIL: Reaction reactants did not match reactant templates.")
        return 1

    products = rxn.RunReactants((acid, amine))
    if not products or len(products) == 0:
        print("FAIL: Reaction generated no products.")
        return 1

    prod_mol = products[0][0]
    Chem.SanitizeMol(prod_mol)
    prod_smiles = Chem.MolToSmiles(prod_mol, canonical=True)
    expected_smiles = Chem.MolToSmiles(Chem.MolFromSmiles("c1ccccc1C(=O)NC"), canonical=True)

    if prod_smiles != expected_smiles:
        print(f"FAIL: Product structure mismatch ({prod_smiles} != {expected_smiles}).")
        return 1
    print(f"      - Amide coupling product verified: {prod_smiles} (matches expected N-methylbenzamide).")

    print("=" * 70)
    print("SUCCESS: All scientific environment and RDKit diagnostics passed cleanly.")
    print("=" * 70)
    return 0


if __name__ == "__main__":
    sys.exit(run_diagnostics())
