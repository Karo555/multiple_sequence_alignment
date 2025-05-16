from aligner.core import needleman_wunsch_alignment

def build_pairwise_score_matrix(sequences, scoring):
    """
    Compute pairwise Needleman-Wunsch alignment scores for all input sequences.
    Returns a symmetric n x n score matrix.
    """
    n = len(sequences)
    score_matrix = [[0] * n for _ in range(n)]

    for i in range(n):
        for j in range(i + 1, n):
            _, _, score = needleman_wunsch_alignment(
                sequences[i],
                sequences[j],
                scoring.match,
                scoring.mismatch,
                scoring.gap
            )
            score_matrix[i][j] = score
            score_matrix[j][i] = score

    return score_matrix
