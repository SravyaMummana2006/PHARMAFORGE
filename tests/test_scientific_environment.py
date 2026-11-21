"""Scientific tests verifying authentic RDKit execution, parsing, descriptors, and reactions."""
from __future__ import annotations

import pytest
from pharmaforge.core.exceptions import RDKitUnavailableError
from pharmaforge.verifier.filters import SubstructureFilter
from pharmaforge.verifier.oracle import ScientificTestOracle


@pytest.fixture
def oracle():
    return ScientificTestOracle()


@pytest.fixture
def sub_filter():
    return SubstructureFilter()


def test_valid_smiles_acceptance_and_sanitization(oracle):
    """Test that authentic drugs pass RDKit parsing, sanitization, and property calculations."""
    aspirin = "CC(=O)Oc1ccccc1C(=O)O"
    caffeine = "Cn1cnc2c1c(=O)n(c(=O)n2C)C"

    rep_asp = oracle.verify_smiles(aspirin)
    assert rep_asp.is_valid is True
    assert rep_asp.passed_sanitization is True
    assert rep_asp.passed_valence is True
    assert rep_asp.passed_pains is True

    rep_caf = oracle.verify_smiles(caffeine)
    assert rep_caf.is_valid is True
    assert rep_caf.passed_sanitization is True
    assert rep_caf.passed_valence is True


def test_invalid_smiles_syntax_rejection(oracle):
    """Test that syntactically malformed SMILES are deterministically rejected."""
    unbalanced = "CC(C(=O)O"
    invalid_char = "C1CC1?"

    rep_unb = oracle.verify_smiles(unbalanced)
    assert rep_unb.is_valid is False
    assert rep_unb.passed_sanitization is False
    assert any("invalid chemical syntax" in r for r in rep_unb.error_reasons)

    rep_inv = oracle.verify_smiles(invalid_char)
    assert rep_inv.is_valid is False


def test_invalid_valence_rejection(oracle):
    """Test that non-physical molecules violating atomic valence are caught and rejected."""
    # Pentavalent nitrogen (neutral nitrogen with 5 single bonds)
    pentavalent_n = "CN(C)(C)(C)C"
    rep = oracle.verify_smiles(pentavalent_n)

    assert rep.is_valid is False
    assert (not rep.passed_sanitization) or (not rep.passed_valence)
    assert any("valence" in r.lower() for r in rep.error_reasons)


def test_real_descriptor_calculations(oracle):
    """Test that calculated properties match authentic chemoinformatics reference values."""
    aspirin = "CC(=O)Oc1ccccc1C(=O)O"
    report = oracle.verify_smiles(aspirin)
    props = report.property_scores

    # Exact chemical constants for Aspirin (C9H8O4)
    assert 180.15 < props["molecular_weight"] < 180.17
    assert 1.30 < props["logp"] < 1.32
    assert 63.5 < props["tpsa"] < 63.7
    assert 0.54 < props["qed"] < 0.56

    # Ertl SAScore (1.58 for Aspirin via pubchem fragment database)
    if oracle.is_sascore_available:
        assert 1.55 < props["sascore"] < 1.62


def test_morgan_fingerprint_generation(oracle):
    """Test generating 2048-bit Morgan fingerprints (ECFP4 equivalent)."""
    aspirin = "CC(=O)Oc1ccccc1C(=O)O"
    fp_bits = oracle.generate_morgan_fingerprint(aspirin, radius=2, n_bits=2048)

    assert len(fp_bits) == 2048
    assert all(b in (0, 1) for b in fp_bits)
    num_on = sum(fp_bits)
    assert 20 <= num_on <= 30  # Aspirin typically has 24-25 active bits


def test_pains_filter_catalog_structural_matching(sub_filter):
    """Test that RDKit's official FilterCatalog accurately identifies authentic PAINS structures."""
    # Ene-rhodanine (classic Baell PAINS alert)
    ene_rhodanine = "O=C1NC(=S)SC1=Cc1ccccc1"
    passed, alerts = sub_filter.check_pains(ene_rhodanine)

    assert passed is False
    assert len(alerts) >= 1
    assert any("ene_rhod" in a for a in alerts)

    # Clean non-alert molecule (Aspirin)
    aspirin = "CC(=O)Oc1ccccc1C(=O)O"
    asp_passed, asp_alerts = sub_filter.check_pains(aspirin)
    assert asp_passed is True
    assert len(asp_alerts) == 0


def test_brenk_filter_catalog_structural_matching(sub_filter):
    """Test that Brenk structural alerts are detected via official FilterCatalog."""
    nitro_aromatic = "c1ccc(cc1)[N+](=O)[O-]"
    passed, alerts = sub_filter.check_brenk(nitro_aromatic)

    assert passed is False
    assert len(alerts) >= 1
    assert any("nitro" in a.lower() for a in alerts)


def test_reaction_smarts_execution(oracle):
    """Test full execution, reactant matching, sanitization, and verification of Reaction SMARTS."""
    # Amide coupling: Carboxylic acid + Amine -> Amide
    rxn_smarts = "[C:1](=[O:2])[OH].[N:3][C:4]>>[C:1](=[O:2])[N:3][C:4]"
    benzoic_acid = "c1ccccc1C(=O)O"
    methylamine = "NC"

    products = oracle.execute_reaction(rxn_smarts, [benzoic_acid, methylamine])
    assert len(products) >= 1
    # Canonical SMILES for N-methylbenzamide is CNC(=O)c1ccccc1
    assert "CNC(=O)c1ccccc1" in products or "c1ccccc1C(=O)NC" in products


def test_reaction_smarts_invalid_reactant_count_raises(oracle):
    """Test that providing wrong number of reactants raises explicit ValueError."""
    rxn_smarts = "[C:1](=[O:2])[OH].[N:3][C:4]>>[C:1](=[O:2])[N:3][C:4]"
    with pytest.raises(ValueError, match="requires 2 reactants"):
        oracle.execute_reaction(rxn_smarts, ["c1ccccc1C(=O)O"])


def test_rdkit_dependency_failure_raises_explicit_error(monkeypatch):
    """Test that absence of RDKit raises RDKitUnavailableError and NOT silent heuristic fallback."""
    import pharmaforge.verifier.oracle as oracle_module

    # Temporarily simulate RDKit being unimportable
    original_init = oracle_module.ScientificTestOracle._init_rdkit

    def mock_init_failure(self):
        raise RDKitUnavailableError("Simulated missing RDKit")

    monkeypatch.setattr(oracle_module.ScientificTestOracle, "_init_rdkit", mock_init_failure)

    with pytest.raises(RDKitUnavailableError, match="RDKit is mandatory"):
        oracle_module.ScientificTestOracle()
