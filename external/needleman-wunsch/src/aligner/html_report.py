from typing import List, Dict, Optional
from jinja2 import Template
from aligner.models import Sequence

_HTML_TEMPLATE = """<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>Needleman–Wunsch Alignment Report</title>
    <style>
      body { font-family: Arial, sans-serif; margin: 2em; }
      pre { background: #f8f8f8; padding: 1em; }
      ul, li { margin: 0.5em 0; }
      img { max-width: 100%; height: auto; }
    </style>
</head>
<body>
  <h1>Needleman–Wunsch Alignment Report</h1>
  <h2>Parameters</h2>
  <ul>
    <li>Match: {{ parameters.match }}</li>
    <li>Mismatch: {{ parameters.mismatch }}</li>
    <li>Gap penalty: {{ parameters.gap }}</li>
  </ul>
  <h2>Sequences</h2>
  <pre>>
Sequence 1: {{ seq1_id }}: {{ seq1_seq }}
Sequence 2: {{ seq2_id }}: {{ seq2_seq }}</pre>
  <h2>Alignments</h2>
  {% for path in alignments %}
    <h3>Path {{ loop.index }}</h3>
    <pre>{{ path.aligned_seq1 }}
{{ path.aligned_seq2 }}</pre>
    <ul>
      <li>Length: {{ path.length }}</li>
      <li>Matches: {{ path.matches }} ({{ "%.2f"|format(path.identity_pct) }}%)</li>
      <li>Total gaps: {{ path.gaps }}</li>
    </ul>
  {% endfor %}
  {% if image_path %}
    <h2>Score Matrix Heatmap</h2>
    <img src="{{ image_path }}" alt="Score matrix heatmap">
  {% endif %}
</body>
</html>
"""


def format_html_report(
    seq1: Sequence,
    seq2: Sequence,
    alignments: List[Dict],
    parameters: Dict,
    image_path: Optional[str] = None,
) -> str:
    """
    Render an HTML summary report.
    """
    tmpl = Template(_HTML_TEMPLATE)
    return tmpl.render(
        seq1_id=seq1.id,
        seq1_seq=seq1.sequence,
        seq2_id=seq2.id,
        seq2_seq=seq2.sequence,
        parameters=parameters,
        alignments=alignments,
        image_path=image_path,
    )
