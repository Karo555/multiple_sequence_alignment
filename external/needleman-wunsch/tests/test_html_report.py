import json
from aligner.models import Sequence
from aligner.io import create_output_dict
from aligner.html_report import format_html_report


def test_format_html_report(tmp_path):
    seq1 = Sequence("s1", "A")
    seq2 = Sequence("s2", "A")
    matrix = [[0, -1], [-1, 1]]
    alignments = [("A", "A")]
    data = create_output_dict(
        seq1, seq2, matrix, alignments, match=1, mismatch=-1, gap=-1
    )

    html = format_html_report(
        seq1, seq2, data["alignments"], data["parameters"], image_path="heatmap.png"
    )
    assert "<h1>Needleman" in html
    assert "Match: 1" in html
    assert "Sequence 1: s1: A" in html
    assert "Path 1" in html
    assert 'src="heatmap.png"' in html
