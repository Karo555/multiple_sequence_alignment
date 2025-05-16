from src.aligner.core import needleman_wunsch_alignment

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

def normalize_sequences(seq_str):
    """
    Splits and uppercases space-separated sequences from --input.
    """
    return [seq.strip().upper() for seq in seq_str.strip().split() if seq.strip()]


def detect_sequence_type(seqs):
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


def validate_sequences(seqs, seq_type):
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
        
def parse_fasta_file(filepath):
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

def convert_scores_to_distances(score_matrix):
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
                distance_matrix[i][j] = 0  # distance to self

    return distance_matrix

def find_center_sequence(distance_matrix):
    """
    Finds the index of the sequence with the smallest total distance to others.
    """
    total_distances = [sum(row) for row in distance_matrix]
    center_index = total_distances.index(min(total_distances))
    return center_index
