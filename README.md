# MSA - MULTIPLE SEQUENCE ALIGNMENT (center star method)<br>
This tool allows you to perform multiple sequence alignment (MSA) using the Center Star method. It supports DNA sequences, allows you to visualize results, analyze statistics, and export results to text and image files. Both a graphical (GUI) and a console (CLI) interface are available.
## Pre-requirements
Python 3.12+ and [UV](https://github.com/astral-sh/uv) as package manager <br>

## installation 
`git clone https://github.com/Karo555/multiple_sequence_alignment` <br>
`cd multiple_sequence_alignment` <br>

## setup virtual env
`uv sync` <br>
`source .venv/bin/activate` <br>

## install in development mode
`uv pip install -e .` <br>

## Functions
- **Loading sequences** from FASTA files or directly from the text field
- **Automatic detection of sequence type** (DNA)
- **Configurable scoring parameters** (match, mismatch, gap)
- **Alignment visualization** in text and graphical form (color blocks)
- **Calculation of MSA statistics** (identity, number of matches, mismatches, gaps)
- **Exporting results** to text file and PNG image
- **Responsive GUI** (Tkinter)
- **CLI mode** for running from the command line

## Run the GUI
`python -m src.gui` <br>
<br>
--> example of evolutionaly related sequences to use: <br>
ATGCGTAGCTAGCTAGCTAG<br>
ATGCGTAGCTAGCTAGCTAG<br>
ATGCGTAGCTAGCTAGCAGG<br>
ATGCGTAGCTAGCTAGCCAG<br>

--> example of evolutionaly unrelated sequences to use: <br>
GATTACAGATTACAGATTAG<br>
CCGGAATTCGGCCGAAATTC<br>
TACGTCAGCTACGTAGCTGA<br>
TGCATGCTAGCTACGATCG<br>

## Run the CLI
input as argument -> `python src/cli.py --input "ACGT AGCT GTAC"` <br>
path to input data -> `python src/cli.py --file path/to/sequences.fasta` <br>

path to example data stored in the /data directory <br>
--> data/multiple_sequences_dna.fasta

configurable arguments:
- `--input` – space-separated sequences, e.g. `"ACGT AGCT GTAC"`
- `--file` – path to FASTA file
- `--match` – value for match (default 1)
- `--mismatch` – mismatch penalty (default -1)
- `--gap` – gap penalty (default -2)

## Results
- **GUI:** Results can be saved to a text file or the alignment can be exported as a PNG image ("File" menu or buttons below the results).
- **CLI:** Results are automatically saved to `results/msa_output.txt`.

## license
This project is licensed under the MIT License - see the [LICENSE](https://github.com/Karo555/{REPO_NAME}/LICENSE) file for details. <br>
