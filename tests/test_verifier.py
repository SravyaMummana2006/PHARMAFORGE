"""Unit tests for the Scientific Test-Oracle (STO-IV) and chemical filters."""
import pytest
from pharmaforge.verifier.oracle import ScientificTestOracle
from pharmaforge.verifier.filters import SubstructureFilter


@pytest.fixture
def oracle():
    return ScientificTestOracle()


@pytest.fixture
def sub_filter():
    return SubstructureFilter()


def test_verify_valid_molecules(oracle):
    """Test that standard valid drug molecules pass verification."""
    # Aspirin
    aspirin_smiles = "CC(=O)Oc1ccccc1C(=O)O"
    report = oracle.verify_smiles(aspirin_smiles)
    assert report.is_valid is True
    assert "molecular_weight" in report.property_scores
    assert report.property_scores["molecular_weight"] > 150.0
    assert report.property_scores["qed"] > 0.35


def test_verify_invalid_syntax(oracle):
    """Test that syntactically malformed SMILES are rejected."""
    bad_smiles = "CC(C(=O)O"  # Unbalanced parenthesis
    report = oracle.verify_smiles(bad_smiles)
    assert report.is_valid is False
    assert len(report.error_reasons) > 0


def test_verify_empty_smiles(oracle):
    """Test rejection of empty input."""
    report = oracle.verify_smiles("")
    assert report.is_valid is False


def test_pains_filter_detection(sub_filter):
    """Test that PAINS toxicophore patterns are detected."""
    # Epoxide (reactive group)
    epoxide_smiles = "CC1OC1"
    passed, alerts = sub_filter.check_pains(epoxide_smiles)
    assert passed is False
    assert len(alerts) >= 1


def test_create_candidate_valid(oracle):
    """Test wrapping validated SMILES into MoleculeCandidate."""
    cand = oracle.create_candidate("c1ccccc1", affinity_score=5.5)
    assert cand.smiles == "c1ccccc1"
    assert cand.affinity_score == 5.5
    assert cand.molecular_weight is not None


def test_create_candidate_invalid_raises(oracle):
    """Test that creating candidate from invalid SMILES raises ValueError."""
    with pytest.raises(ValueError):
        oracle.create_candidate("INVALID_SMILES_@@@@")
