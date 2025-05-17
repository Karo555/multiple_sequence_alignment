import tkinter as tk
from tkinter import ttk, messagebox
import sys, os
from tkinter import filedialog
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from utils.functions import *
from aligner.models import Sequence

def draw_alignment_blocks(final_msa, seq_objects):
    alignment_canvas.delete("all")
    cell_width = 20
    cell_height = 20
    padding = 10

    colors = {
        "match": "#c8e6c9",
        "mismatch": "#ffcdd2",
        "gap": "#e0e0e0"
    }

    for row_idx, aligned_seq in enumerate(final_msa):
        y = row_idx * cell_height + padding
        alignment_canvas.create_text(5, y + cell_height // 2, text=seq_objects[row_idx].id, anchor="w", font=("Courier", 10, "bold"))

        for col_idx, char in enumerate(aligned_seq):
            column = [seq[col_idx] for seq in final_msa]
            if all(c == '-' for c in column):
                color = colors["gap"]
            elif '-' in column:
                color = colors["gap"]
            elif all(c == column[0] for c in column):
                color = colors["match"]
            else:
                color = colors["mismatch"]

            x = 100 + col_idx * cell_width
            alignment_canvas.create_rectangle(x, y, x + cell_width, y + cell_height, fill=color, outline="gray")
            alignment_canvas.create_text(x + cell_width // 2, y + cell_height // 2, text=char, font=("Courier", 10))

    canvas_width = max(800, 100 + len(final_msa[0]) * cell_width)
    alignment_canvas.configure(scrollregion=(0, 0, canvas_width, len(final_msa) * cell_height + padding))

def load_fasta_file():
    filepaths = filedialog.askopenfilenames(filetypes=[("FASTA files", "*.fasta *.fa"), ("All files", "*")])
    if not filepaths:
        return
    try:
        all_sequences = []
        for filepath in filepaths:
            sequences = parse_fasta_file(filepath)
            all_sequences.extend(sequences)
        for seq in all_sequences:
            sequence_input.insert(tk.END, seq.strip() + "\n")
    except Exception as e:
        messagebox.showerror("File Load Error", str(e))


def save_output_to_file():
    filepath = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text files", "*.txt"), ("All files", "*")])
    if not filepath:
        return
    try:
        content = result_text.get("1.0", tk.END).strip()
        with open(filepath, "w") as f:
            f.write(content)
        messagebox.showinfo("Saved", f"Results saved to {filepath}")
    except Exception as e:
        messagebox.showerror("Save Error", str(e))
        
def run_alignment():
    try:
        raw_input = sequence_input.get("1.0", tk.END).strip()
        sequences = normalize_sequences(raw_input.replace("\n", " "))

        if not sequences:
            raise ValueError("Please enter at least one valid sequence.")

        seq_type = seq_type_var.get()
        if not seq_type:
            sequence_type = detect_sequence_type(sequences)
        else:
            sequence_type = seq_type

        validate_sequences(sequences, sequence_type)

        seq_objects = [
            Sequence(f"seq{i+1}", seq, alphabet=sequence_type)
            for i, seq in enumerate(sequences)
        ]

        scoring = ScoringScheme(
            match=int(match_var.get()),
            mismatch=int(mismatch_var.get()),
            gap=int(gap_var.get())
        )

        score_matrix = build_pairwise_score_matrix(seq_objects, scoring)
        distance_matrix = convert_scores_to_distances(score_matrix)
        center_index = find_center_sequence(distance_matrix)

        aligned = align_all_to_center(seq_objects, center_index, scoring)
        final_msa = merge_alignments_to_msa(aligned, center_index)
        stats = compute_msa_statistics(final_msa)

        # Display results
        result_text.delete("1.0", tk.END)
        result_text.insert(tk.END, f"Scoring: {scoring}\n")
        result_text.insert(tk.END, f"Detected type: {sequence_type}\n")
        result_text.insert(tk.END, f"Center sequence: {seq_objects[center_index].id}\n\n")
        result_text.insert(tk.END, "Alignment:\n")
        for i, aligned_seq in enumerate(final_msa):
            result_text.insert(tk.END, f"{seq_objects[i].id}: {aligned_seq}\n")

        result_text.insert(tk.END, "\nStatistics:\n")
        for k, v in stats.items():
            result_text.insert(tk.END, f"{k.replace('_', ' ').capitalize()}: {v}\n")

        # Draw graphical alignment block
        draw_alignment_blocks(final_msa, seq_objects)

    except Exception as e:
        messagebox.showerror("Error", str(e))


# GUI setup
root = tk.Tk()
root.title("MSA - Center Star Method")
root.geometry("800x650")

mainframe = ttk.Frame(root, padding="10")
mainframe.grid(row=0, column=0, sticky=("N", "W", "E", "S"))


# Sequence input
ttk.Label(mainframe, text="Enter sequences (one per line):").grid(row=0, column=0, sticky="W")
sequence_input = tk.Text(mainframe, height=6, width=80)
sequence_input.grid(row=1, column=0, columnspan=3, pady=5)

style = ttk.Style()
style.configure("Purple.TButton", background="#a259e6", foreground="black", font=("Arial", 12, "bold"))
ttk.Button(mainframe, text="Load from FASTA", command=load_fasta_file, style="Purple.TButton").grid(row=2, column=0, sticky="W", pady=(0, 10))

# Sequence type selection
ttk.Label(mainframe, text="Sequence Type:").grid(row=3, column=0, sticky="E")
seq_type_var = tk.StringVar()
type_dropdown = ttk.Combobox(mainframe, textvariable=seq_type_var, values=["", "dna", "rna", "protein"], state="readonly")
type_dropdown.grid(row=2, column=1, sticky="W")
type_dropdown.set("")

# Scoring inputs
ttk.Label(mainframe, text="Match:").grid(row=3, column=0, sticky="E")
match_var = tk.StringVar(value="1")
ttk.Entry(mainframe, textvariable=match_var, width=5).grid(row=3, column=1, sticky="W")

ttk.Label(mainframe, text="Mismatch:").grid(row=4, column=0, sticky="E")
mismatch_var = tk.StringVar(value="-1")
ttk.Entry(mainframe, textvariable=mismatch_var, width=5).grid(row=4, column=1, sticky="W")

ttk.Label(mainframe, text="Gap:").grid(row=5, column=0, sticky="E")
gap_var = tk.StringVar(value="-2")
ttk.Entry(mainframe, textvariable=gap_var, width=5).grid(row=5, column=1, sticky="W")

# Run button
ttk.Button(mainframe, text="RUN", command=run_alignment, style="Purple.TButton").grid(row=6, column=0, columnspan=2, pady=10)

# Output display
result_text = tk.Text(mainframe, height=20, width=90)
result_text.grid(row=7, column=0, columnspan=3, pady=10)

# Graphical block alignment display pane
canvas_frame = ttk.Frame(mainframe)
canvas_frame.grid(row=8, column=0, columnspan=3, sticky="NSEW", pady=(10, 0))

alignment_canvas = tk.Canvas(canvas_frame, width=780, height=200, bg="white")
alignment_canvas.pack(side="left", fill="both", expand=True)

canvas_scroll = ttk.Scrollbar(canvas_frame, orient="horizontal", command=alignment_canvas.xview)


#output file button
ttk.Button(mainframe, text="Save Output", command=save_output_to_file, style="Purple.TButton").grid(row=9, column=0, sticky="E")

root.mainloop()