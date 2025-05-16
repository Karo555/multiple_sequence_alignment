import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'external', 'needleman-wunsch', 'src')))
from typing import List
from aligner.core import needleman_wunsch_alignment
from aligner.models import Sequence


class ScoringScheme:
    def __init__(self, match=1, mismatch=-1, gap=-2):
        self.match = match
        self.mismatch = mismatch
        self.gap = gap

    def __repr__(self):
        return f"Scoring(match={self.match}, mismatch={self.mismatch}, gap={self.gap})"


VALID_CHARS = {
    "dna": set("ACGT"),
    "rna": set("ACGU"),
    "protein": set("ACDEFGHIKLMNPQRSTVWY")
}


def normalize_sequences(seq_str: str) -> List[str]:
    """
    Splits and uppercases space-separated sequences from --input.
    """
    return [seq.strip().upper() for seq in seq_str.strip().split() if seq.strip()]


def detect_sequence_type(seqs: List[str]) -> str:
    """
    Auto-detect the type of sequences (dna, rna, protein).
    Returns a string or raises ValueError if ambiguous or invalid.
    """
    matching_types = []

    for seq_type, valid_chars in VALID_CHARS.items():
        if all(set(seq).issubset(valid_chars) for seq in seqs):
            matching_types.append(seq_type)

    if len(matching_types) == 1:
        return matching_types[0]
    elif not matching_types:
        raise ValueError("Sequences contain invalid characters for all known types.")
    else:
        raise ValueError(f"Ambiguous sequence type detected: matches {matching_types}")


def validate_sequences(seqs: List[str], seq_type: str) -> None:
    """
    Validates sequences against the selected or detected sequence type.
    Raises ValueError if any invalid characters are found.
    """
    valid_chars = VALID_CHARS[seq_type]

    for i, seq in enumerate(seqs):
        invalid = set(seq) - valid_chars
        if invalid:
            raise ValueError(
                f"Sequence {i + 1} contains invalid characters for {seq_type.upper()}: {invalid}"
            )


def parse_fasta_file(filepath: str) -> List[str]:
    """
    Parses a FASTA file and returns a list of uppercased sequences.
    Supports multi-line FASTA entries.
    """
    sequences = []
    current_seq = []

    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if current_seq:
                        sequences.append(''.join(current_seq).upper())
                        current_seq = []
                    continue
                else:
                    current_seq.append(line)

            if current_seq:
                sequences.append(''.join(current_seq).upper())

    except FileNotFoundError:
        raise FileNotFoundError(f"FASTA file not found: {filepath}")
    except Exception as e:
        raise RuntimeError(f"Error reading FASTA file: {e}")

    if not sequences:
        raise ValueError("No sequences found in FASTA file.")

    return sequences


def build_pairwise_score_matrix(sequences: List[Sequence], scoring: ScoringScheme) -> List[List[int]]:
    """
    Compute pairwise Needleman-Wunsch alignment scores for all input sequences.
    Returns a symmetric n x n score matrix.
    """
    n = len(sequences)
    score_matrix = [[0] * n for _ in range(n)]

    for i in range(n):
        for j in range(i + 1, n):
            _, _, score = needleman_wunsch_alignment(
                sequences[i].sequence,
                sequences[j].sequence,
                scoring.match,
                scoring.mismatch,
                scoring.gap
            )
            score_matrix[i][j] = score
            score_matrix[j][i] = score

    return score_matrix


def convert_scores_to_distances(score_matrix: List[List[int]]) -> List[List[int]]:
    """
    Converts a score matrix into a distance matrix using:
    distance = max_score - score
    """
    n = len(score_matrix)
    max_score = max(
        score_matrix[i][j]
        for i in range(n)
        for j in range(n)
        if i != j
    )

    distance_matrix = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                distance_matrix[i][j] = max_score - score_matrix[i][j]
            else:
                distance_matrix[i][j] = 0

    return distance_matrix


def find_center_sequence(distance_matrix: List[List[int]]) -> int:
    """
    Finds the index of the sequence with the smallest total distance to others.
    """
    total_distances = [sum(row) for row in distance_matrix]
    center_index = total_distances.index(min(total_distances))
    return center_index


def align_all_to_center(sequences: List[Sequence], center_index: int, scoring: ScoringScheme) -> List[str]:
    """
    Aligns all sequences to the center sequence using Needleman-Wunsch.
    Returns a list of aligned sequences (with gaps), preserving the order of input.
    """
    center_seq = sequences[center_index].sequence
    aligned_sequences = [None] * len(sequences)

    aligned_sequences[center_index] = center_seq

    for i, seq in enumerate(sequences):
        if i == center_index:
            continue
        aligned_center, aligned_other, _ = needleman_wunsch_alignment(
            center_seq,
            seq.sequence,
            scoring.match,
            scoring.mismatch,
            scoring.gap
        )
        aligned_sequences[center_index] = aligned_center  # updated alignment
        aligned_sequences[i] = aligned_other

    return aligned_sequences

def merge_alignments_to_msa(aligned_seqs: List[str], center_index: int) -> List[str]:
    """
    Ensures all sequences align with the final gapped center sequence.
    Inserts gaps as needed so that all sequences share the same gap structure.
    """
    final_center = aligned_seqs[center_index]
    msa = [''] * len(aligned_seqs)

    for pos in range(len(final_center)):
        if final_center[pos] == '-':
            # Insert a gap at this position in all sequences
            for i in range(len(aligned_seqs)):
                if len(aligned_seqs[i]) <= pos or aligned_seqs[i][pos] != '-':
                    msa[i] += '-'
                else:
                    msa[i] += aligned_seqs[i][pos]
        else:
            # Add the actual character from each sequence
            for i in range(len(aligned_seqs)):
                msa[i] += aligned_seqs[i][pos]

    return msa
  
def compute_msa_statistics(aligned_seqs: List[str]) -> dict:
    """
    Computes MSA statistics: identity %, match, mismatch, gap counts.
    Assumes aligned sequences of equal length.
    """
    length = len(aligned_seqs[0])
    num_seqs = len(aligned_seqs)
    
    matches = 0
    mismatches = 0
    gaps = 0

    for pos in range(length):
        column = [seq[pos] for seq in aligned_seqs]
        if all(char == '-' for char in column):
            continue  # skip all-gap columns
        elif '-' in column:
            gaps += 1
        elif all(char == column[0] for char in column):
            matches += 1
        else:
            mismatches += 1

    identity = (matches / length) * 100 if length > 0 else 0.0

    return {
        "alignment_length": length,
        "matches": matches,
        "mismatches": mismatches,
        "gaps": gaps,
        "identity_percent": round(identity, 2)
    }

import os

def save_alignment_output(
    output_path: str,
    aligned_seqs: List[str],
    sequence_ids: List[str],
    scoring: ScoringScheme,
    center_index: int
):
    """
    Saves the final alignment and scoring info to a text file.
    """
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    with open(output_path, "w") as f:
        f.write("# MSA Output\n")
        f.write(f"Scoring: match={scoring.match}, mismatch={scoring.mismatch}, gap={scoring.gap}\n")
        f.write(f"Center Sequence Index: {center_index}\n\n")
        f.write("Alignment:\n")
        for id_, seq in zip(sequence_ids, aligned_seqs):
            f.write(f"{id_}: {seq}\n")
        f.write("\nStatistics:\n")
        stats = compute_msa_statistics(aligned_seqs)
        for key, value in stats.items():
            f.write(f"{key.replace('_', ' ').capitalize()}: {value}\n")