import sys
import os
import argparse
from utils.functions import *
from aligner.models import Sequence

# Ensure submodule path is added for aligner imports
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'external', 'needleman-wunsch', 'src')))


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Multiple Sequence Alignment (MSA) using the Center Star method"
    )

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--input",
        type=str,
        help="Space-separated sequences passed directly as a single quoted string (e.g., 'ACGT AGCT GTAC')"
    )
    input_group.add_argument(
        "--file",
        type=str,
        help="Path to a FASTA file containing sequences"
    )

    parser.add_argument(
        "--type",
        choices=["dna", "rna"],
        help="Specify the sequence type if known (dna)"
    )

    parser.add_argument(
        "--match",
        type=float,
        default=1,
        help="Score for a match (default: 1)"
    )
    parser.add_argument(
        "--mismatch",
        type=float,
        default=-1,
        help="Penalty for a mismatch (default: -1)"
    )
    parser.add_argument(
        "--gap",
        type=float,
        default=-2,
        help="Penalty for a gap (default: -2)"
    )

    return parser.parse_args()


def main():
    args = parse_arguments()
    scoring = ScoringScheme(match=args.match, mismatch=args.mismatch, gap=args.gap)
    print("Scoring Scheme:", scoring)

    # Step 1: Read sequences
    if args.input:
        raw_sequences = normalize_sequences(args.input)
    else:
        raw_sequences = parse_fasta_file(args.file)

    # Step 2: Detect or validate sequence type
    if args.type:
        sequence_type = args.type
        validate_sequences(raw_sequences, sequence_type)
    else:
        sequence_type = detect_sequence_type(raw_sequences)

    print("Input method:", "Direct Input" if args.input else "FASTA File")
    print("Validated Sequences:", raw_sequences)
    print("Detected type:", sequence_type)

    # Step 3: Create Sequence objects
    sequences = [
        Sequence(f"seq{i+1}", seq, alphabet=sequence_type)
        for i, seq in enumerate(raw_sequences)
    ]

    # Step 4: Build score + distance matrices
    score_matrix = build_pairwise_score_matrix(sequences, scoring)
    distance_matrix = convert_scores_to_distances(score_matrix)
    center_index = find_center_sequence(distance_matrix)

    print(f"Center sequence index: {center_index}")
    print(f"Center sequence: {sequences[center_index].sequence}")

    # Step 5: Align all to center
    aligned_seqs = align_all_to_center(sequences, center_index, scoring)

    print("\nAligned Sequences (center-aligned):")
    for i, aligned in enumerate(aligned_seqs):
        print(f"{sequences[i].id}: {aligned}")

    # Step 6: MSA merging logic
    final_msa = merge_alignments_to_msa(aligned_seqs, center_index)

    print("\nFinal Multiple Sequence Alignment:")
    for i, aligned in enumerate(final_msa):
        print(f"{sequences[i].id}: {aligned}")

    # step 7: compute statistics
    stats = compute_msa_statistics(final_msa)

    print("\nAlignment Statistics:")
    for key, value in stats.items():
        print(f"{key.replace('_', ' ').capitalize()}: {value}")

    output_file = "results/msa_output.txt"
    save_alignment_output(
    output_file,
    final_msa,
    [s.id for s in sequences],
    scoring,
    center_index
    )
    print(f"\nAlignment saved to {output_file}")

    # Step 8: Save MSA to file
    output_file = "results/msa_output.txt"
    save_alignment_output(
    output_file,
    final_msa,
    [s.id for s in sequences],
    scoring,
    center_index
    )
    print(f"\nAlignment saved to {output_file}")

if __name__ == "__main__":
    main()
