from typing import List, Tuple
from aligner.models import Sequence


def build_score_matrix(
    seq1: Sequence, seq2: Sequence, match: int, mismatch: int, gap: int
) -> List[List[int]]:
    """
    Build and return the scoring matrix for global alignment using
    the Needleman–Wunsch algorithm.

    :param seq1: First sequence
    :param seq2: Second sequence
    :param match: Score for a match
    :param mismatch: Score for a mismatch
    :param gap: Gap penalty (negative)
    :return: A (len(seq1)+1) x (len(seq2)+1) matrix of scores
    """
    n = len(seq1)
    m = len(seq2)
    matrix: List[List[int]] = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        matrix[i][0] = matrix[i - 1][0] + gap

    for j in range(1, m + 1):
        matrix[0][j] = matrix[0][j - 1] + gap

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            char1 = seq1.sequence[i - 1]
            char2 = seq2.sequence[j - 1]
            if char1 == char2:
                diag = matrix[i - 1][j - 1] + match
            else:
                diag = matrix[i - 1][j - 1] + mismatch
            up = matrix[i - 1][j] + gap
            left = matrix[i][j - 1] + gap
            matrix[i][j] = max(diag, up, left)

    return matrix


def traceback(
    matrix: List[List[int]],
    seq1: Sequence,
    seq2: Sequence,
    match: int,
    mismatch: int,
    gap: int,
) -> Tuple[str, str]:
    """
    Perform a traceback through the scoring matrix to recover one optimal alignment.

    :param matrix: Scoring matrix from build_score_matrix
    :param seq1: First sequence
    :param seq2: Second sequence
    :param match: Score for a match
    :param mismatch: Score for a mismatch
    :param gap: Gap penalty
    :return: A tuple of aligned strings (with '-' for gaps)
    """
    i, j = len(seq1), len(seq2)
    aligned1: List[str] = []
    aligned2: List[str] = []

    while i > 0 or j > 0:
        if i > 0 and j > 0:
            char1 = seq1.sequence[i - 1]
            char2 = seq2.sequence[j - 1]
            if char1 == char2:
                score_diag = matrix[i - 1][j - 1] + match
            else:
                score_diag = matrix[i - 1][j - 1] + mismatch
            if matrix[i][j] == score_diag:
                aligned1.append(char1)
                aligned2.append(char2)
                i -= 1
                j -= 1
                continue
        if i > 0 and matrix[i][j] == matrix[i - 1][j] + gap:
            aligned1.append(seq1.sequence[i - 1])
            aligned2.append("-")
            i -= 1
            continue
        if j > 0 and matrix[i][j] == matrix[i][j - 1] + gap:
            aligned1.append("-")
            aligned2.append(seq2.sequence[j - 1])
            j -= 1
            continue
        break

    aligned1.reverse()
    aligned2.reverse()
    return "".join(aligned1), "".join(aligned2)


def trace_all_paths(
    matrix: List[List[int]],
    seq1: Sequence,
    seq2: Sequence,
    match: int,
    mismatch: int,
    gap: int,
    max_paths: int = 100,
) -> List[Tuple[str, str]]:
    """
    Enumerate all optimal alignment paths through the scoring matrix.

    :param matrix: Scoring matrix from build_score_matrix
    :param seq1: First sequence
    :param seq2: Second sequence
    :param match: Score for a match
    :param mismatch: Score for a mismatch
    :param gap: Gap penalty
    :param max_paths: Maximum number of alignments to return
    :return: List of tuples of aligned strings
    """
    n, m = len(seq1), len(seq2)
    paths: List[Tuple[str, str]] = []

    def recurse(i: int, j: int, a1: List[str], a2: List[str]):
        """
        Recursive function to find all paths in the scoring matrix.
        :param i: Current row index
        :param j: Current column index
        :param a1: Current alignment for seq1
        :param a2: Current alignment for seq2
        """
        if len(paths) >= max_paths:
            return
        if i == 0 and j == 0:
            paths.append(("".join(reversed(a1)), "".join(reversed(a2))))
            return
        if i > 0 and j > 0:
            char1 = seq1.sequence[i - 1]
            char2 = seq2.sequence[j - 1]
            score_diag = matrix[i - 1][j - 1] + (match if char1 == char2 else mismatch)
            if matrix[i][j] == score_diag:
                recurse(i - 1, j - 1, a1 + [char1], a2 + [char2])
        if i > 0 and matrix[i][j] == matrix[i - 1][j] + gap:
            recurse(i - 1, j, a1 + [seq1.sequence[i - 1]], a2 + ["-"])
        if j > 0 and matrix[i][j] == matrix[i][j - 1] + gap:
            recurse(i, j - 1, a1 + ["-"], a2 + [seq2.sequence[j - 1]])

    recurse(n, m, [], [])
    return paths
