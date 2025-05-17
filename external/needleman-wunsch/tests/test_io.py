import pytest
import json
from aligner.models import Sequence
from src.aligner.io import (
    format_multi_report,
    format_report,
    write_matrix,
    write_report,
    read_fasta,
    read_manual,
    create_output_dict,
    write_json,
)
from src.aligner.models import Sequence


def test_read_fasta_single(tmp_path):
    content = ">seq1\nACGTACGT\n"
    path = tmp_path / "single.fasta"
    path.write_text(content)
    records = read_fasta(str(path))
    assert len(records) == 1
    seq = records[0]
    assert isinstance(seq, Sequence)
    assert seq.id == "seq1"
    assert seq.sequence == "ACGTACGT"


def test_read_fasta_multiple(tmp_path):
    content = ">s1\nAAA\n>s2\nTTT\n"
    path = tmp_path / "multi.fasta"
    path.write_text(content)
    records = read_fasta(str(path))
    assert [rec.id for rec in records] == ["s1", "s2"]
    assert [rec.sequence for rec in records] == ["AAA", "TTT"]


def test_read_fasta_missing_file():
    with pytest.raises(FileNotFoundError):
        read_fasta("nonexistent.fasta")


def test_read_fasta_malformed(tmp_path):
    path = tmp_path / "bad.fasta"
    path.write_text("ACGTACGT\n")
    with pytest.raises(ValueError):
        read_fasta(str(path))


def test_read_fasta_empty(tmp_path):
    path = tmp_path / "empty.fasta"
    path.write_text("")
    with pytest.raises(ValueError):
        read_fasta(str(path))


def test_read_manual_valid(monkeypatch):
    inputs = iter(["seqA", "acgt", "seqB", "tttt"])
    monkeypatch.setattr("builtins.input", lambda prompt="": next(inputs))
    s1, s2 = read_manual()
    assert isinstance(s1, Sequence) and isinstance(s2, Sequence)
    assert s1.id == "seqA"
    assert s1.sequence == "ACGT"
    assert s2.id == "seqB"
    assert s2.sequence == "TTTT"


def test_read_manual_invalid(monkeypatch):
    inputs = iter(["seqA", "ACGTX", "seqB", "TTTT"])
    monkeypatch.setattr("builtins.input", lambda prompt="": next(inputs))
    with pytest.raises(ValueError):
        read_manual()


def test_write_report_and_format(tmp_path):
    seq1 = Sequence("s1", "A")
    seq2 = Sequence("s2", "A")
    aligned1, aligned2 = "A", "A"
    match, mismatch, gap = 1, -1, -1
    report = format_report(seq1, seq2, aligned1, aligned2, match, mismatch, gap)
    assert "Match score: 1" in report
    assert "Mismatch score: -1" in report
    assert "Gap penalty: -1" in report
    assert "s1: A" in report
    assert "s2: A" in report
    assert "Alignment length: 1" in report
    assert "Identical positions: 1" in report
    out_file = tmp_path / "report.txt"
    write_report(str(out_file), report)
    assert out_file.read_text() == report


def test_format_multi_report(tmp_path):
    seq1 = Sequence("s1", "GA")
    seq2 = Sequence("s2", "AG")
    alignments = [
        ("GA", "AG"),
        ("-GA", "GA-"),
    ]
    report = format_multi_report(seq1, seq2, alignments, match=1, mismatch=-1, gap=-1)

    assert "Needleman–Wunsch Multi‐Path Alignment Report" in report
    assert "match=1, mismatch=-1, gap=-1" in report
    assert "Sequence 1: s1  GA" in report
    assert "Sequence 2: s2  AG" in report

    assert "Path 1:" in report
    assert "GA" in report.split("Path 1:")[1]
    assert "AG" in report.split("Path 1:")[1]
    assert "Length: 2" in report
    assert "Identical positions: 0 (0.00%)" in report
    assert "Total gaps: 0" in report

    assert "Path 2:" in report
    assert "-GA" in report.split("Path 2:")[1]
    assert "GA-" in report.split("Path 2:")[1]
    assert "Length: 3" in report
    assert "Identical positions: 0 (0.00%)" in report
    assert "Total gaps: 0" in report

    out_file = tmp_path / "multi_report.txt"
    write_report(str(out_file), report)
    text = out_file.read_text()
    assert "Path 1:" in text and "Path 2:" in text


def test_write_matrix(tmp_path):
    matrix = [
        [0, -1, -2],
        [-1, 1, 0],
        [-2, 0, 2],
    ]
    out_file = tmp_path / "matrix.csv"
    write_matrix(str(out_file), matrix)

    lines = out_file.read_text().splitlines()
    assert lines == [
        "0,-1,-2",
        "-1,1,0",
        "-2,0,2",
    ]


def test_create_output_dict_and_write_json(tmp_path):
    seq1 = Sequence("s1", "A")
    seq2 = Sequence("s2", "A")
    matrix = [[0, -1], [-1, 1]]
    alignments = [("A", "A")]

    data = create_output_dict(
        seq1, seq2, matrix, alignments, match=1, mismatch=-1, gap=-1
    )
    assert data["parameters"] == {"match": 1, "mismatch": -1, "gap": -1}
    assert data["sequences"] == {"s1": "A", "s2": "A"}
    assert data["matrix"] == matrix
    assert isinstance(data["alignments"], list)
    assert data["alignments"][0]["aligned_seq1"] == "A"

    out = tmp_path / "out.json"
    write_json(str(out), data)
    loaded = json.loads(out.read_text())
    assert loaded == data
