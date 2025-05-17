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

def test_unknown_alphabet_raises():
    with pytest.raises(ValueError):
        Sequence(identifier="foo", sequence="ACGT", alphabet="rna")
