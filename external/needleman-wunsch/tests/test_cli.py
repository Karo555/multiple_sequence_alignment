import pytest
import sys
import json
import aligner.cli as cli
from pathlib import Path
from src.aligner.cli import parse_args


def test_parse_args_input():
    args = parse_args(
        [
            "--input",
            "file1.fasta",
            "file2.fasta",
            "--match",
            "2",
            "--mismatch",
            "-1",
            "--gap",
            "-3",
        ]
    )
    assert args.input == ["file1.fasta", "file2.fasta"]
    assert not args.manual
    assert args.match == 2
    assert args.mismatch == -1
    assert args.gap == -3
    assert args.output is None
    assert args.plot is None


def test_parse_args_manual():
    args = parse_args(["--manual"])
    assert args.manual
    assert args.input is None
    assert args.output is None
    assert args.plot is None


def test_help_shows_usage(capsys):
    with pytest.raises(SystemExit) as excinfo:
        parse_args(["--help"])
    assert excinfo.value.code == 0
    captured = capsys.readouterr()
    assert "usage:" in captured.out.lower()


def test_cli_matrix_export(tmp_path, monkeypatch):
    project = Path(__file__).parent.parent
    monkeypatch.syspath_prepend(str(project / "src"))

    data = tmp_path / "data"
    data.mkdir()
    (data / "s1.fasta").write_text(">s1\nA\n")
    (data / "s2.fasta").write_text(">s2\nA\n")

    monkeypatch.chdir(tmp_path)

    out_csv = tmp_path / "matrix_out.csv"
    sys.argv = [
        "aligner.cli",
        "--input",
        "data/s1.fasta",
        "data/s2.fasta",
        "--matrix-out",
        str(out_csv),
    ]

    cli.main()

    text = out_csv.read_text().splitlines()
    assert text == ["0,-2", "-2,1"]


def test_cli_json_output(tmp_path, monkeypatch):
    project = Path(__file__).parent.parent
    monkeypatch.syspath_prepend(str(project / "src"))

    data = tmp_path / "data"
    data.mkdir()
    (data / "s1.fasta").write_text(">s1\nA\n")
    (data / "s2.fasta").write_text(">s2\nA\n")

    monkeypatch.chdir(tmp_path)
    reports = tmp_path / "reports"
    reports.mkdir()
    out_json = reports / "out.json"

    sys.argv = [
        "aligner.cli",
        "--input",
        "data/s1.fasta",
        "data/s2.fasta",
        "--json",
        str(out_json),
    ]

    cli.main()

    assert out_json.exists()
    payload = json.loads(out_json.read_text())
    assert payload["parameters"]["match"] == 1
    assert payload["sequences"]["s1"] == "A"
    assert isinstance(payload["matrix"], list)
    assert isinstance(payload["alignments"], list)


def test_cli_html_output(tmp_path, monkeypatch):

    project = Path(__file__).parent.parent
    monkeypatch.syspath_prepend(str(project / "src"))

    data = tmp_path / "data"
    data.mkdir()
    (data / "s1.fasta").write_text(">s1\nA\n")
    (data / "s2.fasta").write_text(">s2\nA\n")

    monkeypatch.chdir(tmp_path)

    out_html = tmp_path / "report.html"
    sys.argv = [
        "aligner.cli",
        "--input",
        "data/s1.fasta",
        "data/s2.fasta",
        "--html",
        str(out_html),
    ]

    cli.main()

    text = out_html.read_text()
    assert "<h1>Needleman" in text
    assert "Sequence 1: s1: A" in text


def test_cli_pdf_output(tmp_path, monkeypatch):
    project = Path(__file__).parent.parent
    monkeypatch.syspath_prepend(str(project / "src"))
    data = tmp_path / "data"
    data.mkdir()
    (data / "s1.fasta").write_text(">s1\nA\n")
    (data / "s2.fasta").write_text(">s2\nA\n")

    monkeypatch.chdir(tmp_path)

    out_pdf = tmp_path / "report.pdf"
    sys.argv = [
        "aligner.cli",
        "--input",
        "data/s1.fasta",
        "data/s2.fasta",
        "--pdf",
        str(out_pdf),
        "--plot",
        "plots/heatmap.png",
    ]

    cli.main()
    assert out_pdf.exists() and out_pdf.stat().st_size > 0
