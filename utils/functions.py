import re

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
