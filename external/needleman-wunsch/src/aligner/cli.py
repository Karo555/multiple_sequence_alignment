import argparse
import os
from aligner.html_report import format_html_report
from typing import Optional
from aligner.plot import plot_matrix
from aligner.pdf_report import write_pdf
from aligner.core import (
    build_score_matrix,
    traceback as single_traceback,
    trace_all_paths,
)
from aligner.io import (
    create_output_dict,
    format_multi_report,
    format_report,
    read_manual,
    read_fasta,
    write_json,
    write_matrix,
    write_report,
)


def parse_args(args=None):
    """
    Parse command-line arguments for the Needleman–Wunsch aligner.
    """
    parser = argparse.ArgumentParser(
        prog="needleman-wunsch",
        description="Global alignment of two sequences using Needleman–Wunsch algorithm",
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--input",
        nargs=2,
        metavar=("FASTA1", "FASTA2"),
        help="Paths to two input FASTA files",
    )
    group.add_argument(
        "--manual",
        action="store_true",
        help="Enter sequences manually via prompts",
    )
    parser.add_argument(
        "--match",
        type=int,
        default=1,
        help="Score for a match (default: 1)",
    )
    parser.add_argument(
        "--mismatch",
        type=int,
        default=-1,
        help="Score for a mismatch (default: -1)",
    )
    parser.add_argument(
        "--gap",
        type=int,
        default=-2,
        help="Penalty for a gap (default: -2)",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Optional custom filename for text output report",
    )

    parser.add_argument(
        "--all-paths",
        action="store_true",
        help="Enumerate and output all optimal alignments",
    )

    parser.add_argument(
        "--plot",
        type=str,
        default=None,
        help="Optional filename for PNG heatmap output",
    )

    parser.add_argument(
        "--matrix-out",
        type=str,
        default=None,
        help="Filename for raw score matrix output (CSV)",
    )

    parser.add_argument(
        "--json",
        dest="json_out",
        type=str,
        default=None,
        help="Filename for structured JSON output",
    )

    parser.add_argument(
        "--html",
        dest="html_out",
        type=str,
        default=None,
        help="Filename for HTML summary report",
    )

    parser.add_argument(
        "--pdf",
        dest="pdf_out",
        type=str,
        default=None,
        help="Filename for PDF summary report",
    )

    parser.add_argument(
        "--alphabet",
        choices=["dna", "protein"],
        default="dna",
        help="Alphabet for sequences (dna or protein)",
    )

    return parser.parse_args(args)


def main():
    """
    Main function to run the Needleman–Wunsch aligner from the command line.
    """
    args = parse_args()

    if args.output:
        os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    if args.plot:
        os.makedirs(os.path.dirname(args.plot) or ".", exist_ok=True)
    if args.matrix_out:
        os.makedirs(os.path.dirname(args.matrix_out) or ".", exist_ok=True)
    if args.json_out:
        os.makedirs(os.path.dirname(args.json_out) or ".", exist_ok=True)
    if args.html_out:
        os.makedirs(os.path.dirname(args.html_out) or ".", exist_ok=True)

    if args.manual:
        seq1, seq2 = read_manual(args.alphabet)
    else:
        recs1 = read_fasta(args.input[0], args.alphabet)
        recs2 = read_fasta(args.input[1], args.alphabet)
        if len(recs1) != 1 or len(recs2) != 1:
            raise ValueError("Each FASTA must contain exactly one record")
        seq1, seq2 = recs1[0], recs2[0]

    matrix = build_score_matrix(seq1, seq2, args.match, args.mismatch, args.gap)

    if args.matrix_out:
        write_matrix(args.matrix_out, matrix)

    if args.all_paths:
        align_list = trace_all_paths(
            matrix, seq1, seq2, args.match, args.mismatch, args.gap
        )
        report_text = format_multi_report(
            seq1, seq2, align_list, args.match, args.mismatch, args.gap
        )
    else:
        aln1, aln2 = single_traceback(
            matrix, seq1, seq2, args.match, args.mismatch, args.gap
        )
        align_list = [(aln1, aln2)]
        report_text = format_report(
            seq1, seq2, aln1, aln2, args.match, args.mismatch, args.gap
        )

    if args.output:
        write_report(args.output, report_text)
    else:
        # Print report to console when no output file is specified
        print(report_text)

    if args.json_out:
        data = create_output_dict(
            seq1, seq2, matrix, align_list, args.match, args.mismatch, args.gap
        )
        write_json(args.json_out, data)

    if args.html_out:
        img_ref = None
        if args.plot:
            img_ref = os.path.relpath(args.plot, start=os.path.dirname(args.html_out))
        data = create_output_dict(
            seq1, seq2, matrix, align_list, args.match, args.mismatch, args.gap
        )
        html = format_html_report(
            seq1, seq2, data["alignments"], data["parameters"], img_ref
        )
        with open(args.html_out, "w") as f:
            f.write(html)

    if args.plot:
        plot_matrix(matrix, args.plot)

    if args.pdf_out:
        data = create_output_dict(
            seq1, seq2, matrix, align_list, args.match, args.mismatch, args.gap
        )
        write_pdf(
            args.pdf_out,
            seq1,
            seq2,
            data["alignments"],
            data["parameters"],
            args.plot,
        )


if __name__ == "__main__":
    main()
