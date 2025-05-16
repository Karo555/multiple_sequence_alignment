import argparse
import sys
from utils.functions import normalize_sequences, validate_sequences, detect_sequence_type

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
        choices=["dna", "rna", "protein"],
        help="Specify the sequence type if known (dna, rna, protein)"
    )

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_arguments()
    if args.input:
        sequences = normalize_sequences(args.input)
    else:
        sequences = ["ACGT", "AGGT", "ACCT"]  # placeholder until FASTA parser is added

    if args.type:
        sequence_type = args.type
        validate_sequences(sequences, sequence_type)
    else:
        sequence_type = detect_sequence_type(sequences)

    print(f"Detected type: {sequence_type}")


    print("Input method:", "Direct Input" if args.input else "FASTA File")
    if args.input:
        print("Sequences:", args.input)
    else:
        print("FASTA File Path:", args.file)
    print("Sequence Type:", args.type if args.type else "Auto-detect")