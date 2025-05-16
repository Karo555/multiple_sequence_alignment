import pytest
from src.aligner.models import Sequence


def test_valid_dna_sequence():
    seq = Sequence(identifier="seq1", sequence="ACGTACGT")
    assert seq.id == "seq1"
    assert str(seq.sequence) == "ACGTACGT"
    assert len(seq) == 8


def test_invalid_dna_sequence_raises():
    with pytest.raises(ValueError) as exc:
        Sequence(identifier="bad", sequence="ACGTX")
    assert "Invalid characters" in str(exc.value)


def test_valid_protein_sequence():
    prot = Sequence(identifier="prot1", sequence="ARNDCEQGHI", alphabet="protein")
    assert len(prot) == 10


def test_invalid_protein_sequence_raises():
    with pytest.raises(ValueError):
        Sequence(identifier="badprot", sequence="ARNDXE", alphabet="protein")


def test_unknown_alphabet_raises():
    with pytest.raises(ValueError):
        Sequence(identifier="foo", sequence="ACGT", alphabet="rna")
