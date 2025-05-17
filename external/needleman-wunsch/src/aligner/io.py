import csv
import json
from typing import Dict, List, Tuple
from aligner.models import Sequence


def read_fasta(path: str, alphabet: str = "dna") -> list[Sequence]:
    """
    Read a FASTA file and return a list of Sequence objects.
    Each sequence is represented by a Sequence object with the following attributes:
    - id: The sequence ID (header line without the leading '>').
    - sequence: The sequence string (all lines after the header line).
    - alphabet: The alphabet used for the sequence (e.g., "dna")
    The function raises a ValueError if the file is empty or if there are no valid sequences.
    Parameters
    ----------
    path : str
        Path to the FASTA file.
    alphabet : str
        The alphabet used for the sequence (default is "dna").
    Returns
    -------
    list[Sequence]
        A list of Sequence objects representing the sequences in the FASTA file.
    Raises
    ------
    ValueError
        If the file is empty or if there are no valid sequences.
    """
    sequences = []
    header = None
    seq_lines: list[str] = []
    with open(path) as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seq_str = "".join(seq_lines)
                    # pass the raw string, not len(seq_str)
                    sequences.append(Sequence(header, seq_str, alphabet))
                header = line[1:].strip()
                seq_lines = []
            else:
                if header is None:
                    raise ValueError("FASTA format error: data before header")
                seq_lines.append(line)
        if header is not None:
            seq_str = "".join(seq_lines)
            sequences.append(Sequence(header, seq_str, alphabet))
    if not sequences:
        raise ValueError(f"No sequences found in FASTA file: {path}")

    return sequences


def read_manual(alphabet: str = "dna") -> tuple[Sequence, Sequence]:
    """
    Read two sequences from user input.
    Parameters
    ----------
    alphabet : str
        The alphabet used for the sequences (default is "dna").
    Returns
    -------
    tuple[Sequence, Sequence]
        A tuple of Sequence objects representing the sequences read from user inpu
    Raises
    ------
    ValueError
        If the input is invalid or if the sequences are empty.
    """

    id1 = input("Sequence 1 ID: ").strip()
    seq1_str = input("Sequence 1: ").strip()
    id2 = input("Sequence 2 ID: ").strip()
    seq2_str = input("Sequence 2: ").strip()
    return (
        Sequence(id1, seq1_str, alphabet),
        Sequence(id2, seq2_str, alphabet),
    )


def write_report(path: str, report: str) -> None:
    """
    Write the alignment report to a text file.
    The report includes alignment parameters, sequences, alignment, and metrics.
    Parameters
    ----------
    path : str
        The path to the file to write the report to.
    report : str
        The alignment report to write to the file.
    Raises
    ------
    IOError
        If there is an error writing to the file.
    """
    with open(path, "w") as f:
        f.write(report)


def format_report(
    seq1: Sequence,
    seq2: Sequence,
    aligned1: str,
    aligned2: str,
    match: int,
    mismatch: int,
    gap: int,
) -> str:
    """
    Format alignment parameters, sequences, alignment, and metrics into a report string.
    Parameters
    ----------
    seq1, seq2
        The original Sequence objects.
    aligned1, aligned2
        The aligned sequences.
    match, mismatch, gap
        The number of matches, mismatches, and gaps in the alignmen
    Returns
    -------
    str
        The formatted alignment report string.
    Raises
    ------
    ValueError
        If the aligned sequences are empty or if the lengths do not match.
    """
    length = len(aligned1)
    identical = sum(1 for a, b in zip(aligned1, aligned2) if a == b)
    percent = (identical / length * 100) if length > 0 else 0.0
    total_gaps = aligned1.count("-") + aligned2.count("-")

    lines = [
        "Parameters:",
        f"  Match score: {match}",
        f"  Mismatch score: {mismatch}",
        f"  Gap penalty: {gap}",
        "",
        "Sequences:",
        f"  {seq1.id}: {seq1.sequence}",
        f"  {seq2.id}: {seq2.sequence}",
        "",
        "Alignment:",
        f"  {aligned1}",
        f"  {aligned2}",
        "",
        "Statistics:",
        f"  Alignment length: {length}",
        f"  Identical positions: {identical}",
        f"  Percentage identity: {percent:.2f}%",
        f"  Total gaps: {total_gaps}",
    ]
    return "\n".join(lines)


def format_multi_report(
    seq1: Sequence,
    seq2: Sequence,
    alignments: List[Tuple[str, str]],
    match: int,
    mismatch: int,
    gap: int,
) -> str:
    """
    Build a multi‐path alignment report.

    Parameters
    ----------
    seq1, seq2
        The original Sequence objects.
    alignments
        List of (aligned_seq1, aligned_seq2) tuples.
    match, mismatch, gap
        Scoring parameters.
    """
    lines: List[str] = []

    lines.append("Needleman–Wunsch Multi‐Path Alignment Report")
    lines.append(f"Parameters: match={match}, mismatch={mismatch}, gap={gap}")
    lines.append(f"Sequence 1: {seq1.id}  {seq1.sequence}")
    lines.append(f"Sequence 2: {seq2.id}  {seq2.sequence}")
    lines.append("")

    for idx, (aln1, aln2) in enumerate(alignments, start=1):
        length = len(aln1)
        matches = sum(1 for a, b in zip(aln1, aln2) if a == b and a != "-")
        identity_pct = matches / length * 100
        gaps = aln1.count("-") + aln2.count("-")

        lines.append(f"Path {idx}:")
        lines.append(aln1)
        lines.append(aln2)
        lines.append(f"Length: {length}")
        lines.append(f"Identical positions: {matches} ({identity_pct:.2f}%)")
        lines.append(f"Total gaps: {gaps}")
        lines.append("")

    return "\n".join(lines)


def write_matrix(path: str, matrix: List[List[int]]) -> None:
    """
    Write the DP score matrix to CSV.
    Each row of the matrix becomes one line of comma-separated values.
    The first row is the header with column indices.
    The first column is the row index.
    Parameters
    ----------
    path : str
        The path to the file to write the matrix to.
    matrix : list[list[int]]
        The DP score matrix to write to the file.
    Raises
    ------
    IOError
        If there is an error writing to the file.
    """
    with open(path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        for row in matrix:
            writer.writerow(row)


def create_output_dict(
    seq1: Sequence,
    seq2: Sequence,
    matrix: List[List[int]],
    alignments: List[Tuple[str, str]],
    match: int,
    mismatch: int,
    gap: int,
) -> Dict:
    """
    Package everything into a serializable dict:
      - sequences (id → raw string)
      - parameters
      - full score matrix
      - list of alignments with stats
    Parameters
    ----------
    seq1, seq2
        The original Sequence objects.
    matrix
        The DP score matrix.
    alignments
        List of (aligned_seq1, aligned_seq2) tuples.
    match, mismatch, gap
        Scoring parameters.
    Returns
    -------
    dict
        A dictionary containing the sequences, parameters, matrix, and alignments.
    Raises
    ------
    ValueError
        If the sequences are empty or if the lengths do not match.
    """
    paths = []
    for aln1, aln2 in alignments:
        length = len(aln1)
        matches = sum(1 for a, b in zip(aln1, aln2) if a == b and a != "-")
        identity_pct = matches / length * 100
        gaps = aln1.count("-") + aln2.count("-")
        paths.append(
            {
                "aligned_seq1": aln1,
                "aligned_seq2": aln2,
                "length": length,
                "matches": matches,
                "identity_pct": identity_pct,
                "gaps": gaps,
            }
        )

    return {
        "sequences": {seq1.id: seq1.sequence, seq2.id: seq2.sequence},
        "parameters": {"match": match, "mismatch": mismatch, "gap": gap},
        "matrix": matrix,
        "alignments": paths,
    }


def write_json(path: str, data: Dict) -> None:
    """
    Write the dict `data` out as pretty JSON.
    Parameters
    ----------
    path
        The path to write the JSON to.
    data
        The data to write.e_path, "w"
    Raises
    ------
    FileNotFoundError
        If the path does not exist.
    """
    with open(path, "w") as f:
        json.dump(data, f, indent=2)
