# Needleman–Wunsch Aligner
bIOINFORMATICS @WUST

# Installation
`git clone https://github.com/Karo555/needleman-wunsch.git` <br>
`pip install -e ".[html,pdf]"`

# Usage
1. Basic DNA alignment (manual input) <br>
`needleman-wunsch --manual`

You’ll be prompted to enter:
- Sequence 1 ID
- Sequence 1
- Sequence 2 ID
- Sequence 2

2. FASTA file input <br>
`needleman-wunsch --input data/seq1.fasta data/seq2.fasta --output reports/alignment.txt` <br>

--> data/seq1.fasta and data/seq2.fasta each contain exactly one record <br>
--> the text report is saved to /reports <br>

3. Custom scoring. Override match, mismatch, and gap scores. <br>
`needleman-wunsch --input data/seq1.fasta data/seq2.fasta --match 2 --mismatch -1 --gap -2`

4. Enumerate all optimal alignments <br>
`needleman-wunsch --input data/seq1.fasta data/seq2.fasta --all-paths --output reports/all_paths.txt` <br>

--> Lists every equally optimal alignment and write the text report to /reports

5. Export raw DP matrix as CSV <br>
`needleman-wunsch --input data/seq1.fasta data/seq2.fasta --matrix-out reports/matrix.csv` <br>

--> score matrix is saved in .csv at /reports

6. Structured JSON output <br>
`needleman-wunsch --input data/seq1.fasta data/seq2.fasta --json reports/alignment.json`

--> produces a JSON file containing original sequences and IDs, scoring parameters, full DP matrix, all alignments with per-path statistics <br>
--> saved at /reports


7. Heatmap visualization <br>
`needleman-wunsch --input data/seq1.fasta data/seq2.fasta --plot plots/heatmap.png`

--> renders the DP matrix as a heatmap PNG saved at /plots

8. HTML summary report <br>
`needleman-wunsch --input data/seq1.fasta data/seq2.fasta --plot plots/heatmap.png --html reports/alignment.html`

--> generates a standalone HTML page embedding parameters & sequences, all alignments (or just the single optimal one), inline statistics, the heatmap image

--> Open it in your browser via:
`python -m http.server 8000`

then visit `http://localhost:8000/reports/alignment.html`

9. PDF summary report <br>
`needleman-wunsch --input data/seq1.fasta data/seq2.fasta --plot plots/heatmap.png --pdf reports/alignment.pdf`

--> saved to /reports


10. one shot <br>
`needleman-wunsch --input data/seq1.fasta data/seq2.fasta --all-paths --match 2 --mismatch -1 --gap -2 --output reports/alignment.txt --matrix-out reports/matrix.csv --json reports/alignment.json --plot plots/heatmap.png --html reports/alignment.html --pdf reports/alignment.pdf`

reports/alignment.txt (text report)<br>
reports/matrix.csv (raw DP matrix)<br>
reports/alignment.json (structured JSON)<br>
plots/heatmap.png (heatmap image)<br>
reports/alignment.html (HTML summary)<br>
reports/alignment.pdf (PDF summary)<br>


## 📄 License<br>
This project is licensed under the MIT License. See [LICENSE](LICENSE.txt) for details.<br>
