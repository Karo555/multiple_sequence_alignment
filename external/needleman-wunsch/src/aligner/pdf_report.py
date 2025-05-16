from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image as RLImage
from reportlab.platypus import Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.pagesizes import letter
from reportlab.lib import colors
from typing import List, Dict, Optional
from aligner.models import Sequence


def write_pdf(
    path: str,
    seq1: Sequence,
    seq2: Sequence,
    alignments: List[Dict],
    parameters: Dict,
    image_path: Optional[str] = None,
) -> None:
    """
    Create a PDF report with:
      - Parameters
      - Sequences
      - One or more alignment paths (with stats)
      - Optional heatmap image
    Parameters
    ----------
    path
        Path to save the PDF report.
    seq1, seq2
        The original Sequence objects.
    alignments
        List of dictionaries containing alignment information.
    parameters
        Dictionary containing the parameters used to generate the alignmen
        - match
        - mismatch
        - gap
    image_path
        Path to the image to include in the report (optional).
    Returns
    -------
    None
        The function does not return anything. It creates a PDF file at the specified path.
    """
    doc = SimpleDocTemplate(path, pagesize=letter)
    styles = getSampleStyleSheet()
    story = []

    story.append(Paragraph("Needleman–Wunsch Alignment Report", styles["Title"]))
    story.append(Spacer(1, 12))

    param_data = [
        ["Match", parameters["match"]],
        ["Mismatch", parameters["mismatch"]],
        ["Gap penalty", parameters["gap"]],
    ]
    tbl = Table(param_data, colWidths=[100, 50])
    tbl.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.lightgrey),
                ("GRID", (0, 0), (-1, -1), 0.5, colors.grey),
            ]
        )
    )
    story.append(tbl)
    story.append(Spacer(1, 12))

    story.append(
        Paragraph(f"<b>Sequence 1:</b> {seq1.id}: {seq1.sequence}", styles["Normal"])
    )
    story.append(
        Paragraph(f"<b>Sequence 2:</b> {seq2.id}: {seq2.sequence}", styles["Normal"])
    )
    story.append(Spacer(1, 12))

    for idx, path in enumerate(alignments, start=1):
        story.append(Paragraph(f"Path {idx}", styles["Heading2"]))
        story.append(Paragraph(path["aligned_seq1"], styles["Code"]))
        story.append(Paragraph(path["aligned_seq2"], styles["Code"]))
        stats = [
            ["Length", path["length"]],
            ["Matches", f"{path['matches']} ({path['identity_pct']:.2f}%)"],
            ["Gaps", path["gaps"]],
        ]
        t = Table(stats, colWidths=[100, 100])
        t.setStyle(
            TableStyle(
                [
                    ("GRID", (0, 0), (-1, -1), 0.5, colors.grey),
                ]
            )
        )
        story.append(t)
        story.append(Spacer(1, 12))

    if image_path:
        story.append(Paragraph("Score Matrix Heatmap", styles["Heading2"]))
        story.append(RLImage(image_path, width=400, height=400))
        story.append(Spacer(1, 12))

    doc.build(story)
