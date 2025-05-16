import pytest
from src.aligner.models import Sequence
from src.aligner.core import build_score_matrix, trace_all_paths, traceback


def test_build_score_matrix_single_match():
    s1 = Sequence("s1", "A")
    s2 = Sequence("s2", "A")
    mat = build_score_matrix(s1, s2, match=1, mismatch=-1, gap=-1)
    assert len(mat) == 2 and len(mat[0]) == 2
    assert mat[0][0] == 0
    assert mat[0][1] == -1
    assert mat[1][0] == -1
    assert mat[1][1] == 1


def test_build_score_matrix_single_mismatch():
    s1 = Sequence("s1", "A")
    s2 = Sequence("s2", "C")
    mat = build_score_matrix(s1, s2, match=1, mismatch=-1, gap=-1)
    assert mat[1][1] == -1


def test_build_score_matrix_with_gap():
    s1 = Sequence("s1", "A")
    s2 = Sequence("s2", "")
    mat = build_score_matrix(s1, s2, match=1, mismatch=-1, gap=-2)
    assert mat[1][0] == -2


def test_traceback_single_match():
    s1 = Sequence("s1", "A")
    s2 = Sequence("s2", "A")
    mat = build_score_matrix(s1, s2, match=1, mismatch=-1, gap=-1)
    aln1, aln2 = traceback(mat, s1, s2, match=1, mismatch=-1, gap=-1)
    assert aln1 == "A"
    assert aln2 == "A"


def test_traceback_gap():
    s1 = Sequence("s1", "A")
    s2 = Sequence("s2", "")
    mat = build_score_matrix(s1, s2, match=1, mismatch=-1, gap=-2)
    aln1, aln2 = traceback(mat, s1, s2, match=1, mismatch=-1, gap=-2)
    assert aln1 == "A"
    assert aln2 == "-"


def test_trace_all_paths_gap_variants():
    s1 = Sequence("s1", "AG")
    s2 = Sequence("s2", "A")
    mat = build_score_matrix(s1, s2, match=1, mismatch=-1, gap=-1)
    paths = trace_all_paths(mat, s1, s2, match=1, mismatch=-1, gap=-1)
    expected = {("AG", "A-")}
    assert set(paths) == expected
