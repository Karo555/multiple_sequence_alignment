import sys
import os
import inspect
import aligner.core
print(inspect.getmembers(aligner.core, inspect.isfunction))

# Add external/needleman-wunsch/src to sys.path
repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
nw_path = os.path.join(repo_root, 'external', 'needleman-wunsch', 'src')
if nw_path not in sys.path:
    sys.path.append(nw_path)

from aligner.core import needleman_wunsch_alignment
from src.core.build_pairwise_score_matrix import build_pairwise_score_matrix

# Define scoring scheme class (if not already imported)
class ScoringScheme:
    def __init__(self, match=1, mismatch=-1, gap=-2):
        self.match = match
        self.mismatch = mismatch
        self.gap = gap

# Test data
sequences = ["ACGT", "AGCT", "ATGT"]
scoring = ScoringScheme(match=1, mismatch=-1, gap=-2)

# Run test
score_matrix = build_pairwise_score_matrix(sequences, scoring)

# Print matrix
print("Pairwise Score Matrix:")
for row in score_matrix:
    print(row)
