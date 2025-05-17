import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import os
import sys
import webbrowser
from pathlib import Path
import threading
from typing import List, Dict, Tuple, Optional
import json

# Add parent directory to path for importing modules
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Import required modules
from utils.functions import (
    parse_fasta_file, normalize_sequences, detect_sequence_type,
    validate_sequences, ScoringScheme, build_pairwise_score_matrix,
    convert_scores_to_distances, find_center_sequence,
    align_all_to_center, merge_alignments_to_msa, compute_msa_statistics
)
from aligner.models import Sequence


class MSAApplication:
    """Main application class for Multiple Sequence Alignment GUI"""
    
    COLORS = {
        "match": "#c8e6c9",      # Light green
        "mismatch": "#ffcdd2",    # Light red
        "gap": "#e0e0e0",         # Light grey
        "primary": "#6200ee",     # Primary color
        "secondary": "#03dac6",   # Secondary color
        "background": "#f5f5f5",  # Background color
        "surface": "#ffffff",     # Surface color
        "error": "#b00020",       # Error color
        "on_primary": "#ffffff",  # Text on primary color
        "on_secondary": "#000000" # Text on secondary color
    }
    
    DEFAULT_SCORING = {
        "match": 1,
        "mismatch": -1,
        "gap": -2
    }
    
    def __init__(self, root):
        """Initialize the application with the root window"""
        self.root = root
        self.root.title("Enhanced MSA - Center Star Method")
        self.root.geometry("1000x800")
        self.root.minsize(800, 600)
        
        # Variables
        self.seq_type_var = tk.StringVar(value="")
        self.match_var = tk.StringVar(value=str(self.DEFAULT_SCORING["match"]))
        self.mismatch_var = tk.StringVar(value=str(self.DEFAULT_SCORING["mismatch"]))
        self.gap_var = tk.StringVar(value=str(self.DEFAULT_SCORING["gap"]))
        self.status_var = tk.StringVar(value="Ready")
        self.alignment_font_size = 12
        self.last_msa = None
        self.last_seq_objects = None
        self.file_paths = []
        self.theme_mode = tk.StringVar(value="light")
        
        # Set up UI
        self._setup_styles()
        self._create_menu()
        self._create_ui()
        
        # Configure window resize event
        self.root.bind("<Configure>", self._on_window_resize)
        
        # Load settings if available
        self._load_settings()
        
    def _setup_styles(self):
        """Set up ttk styles for the application"""
        self.style = ttk.Style()
        
        # Configure common styles
        self.style.configure("TFrame", background=self.COLORS["background"])
        self.style.configure("TLabel", background=self.COLORS["background"], font=("Segoe UI", 10))
        self.style.configure("TButton", font=("Segoe UI", 10))
        
        # Custom button styles
        self.style.configure(
            "Primary.TButton",
            background=self.COLORS["primary"],
            foreground=self.COLORS["on_primary"],
            font=("Segoe UI", 11, "bold")
        )
        
        self.style.configure(
            "Secondary.TButton",
            background=self.COLORS["secondary"],
            foreground=self.COLORS["on_secondary"],
            font=("Segoe UI", 10)
        )
        
        # Status bar style
        self.style.configure(
            "StatusBar.TLabel",
            background="#f0f0f0",
            relief="sunken",
            font=("Segoe UI", 9)
        )
        
        # Configure the notebook style
        self.style.configure("TNotebook", background=self.COLORS["background"])
        self.style.configure("TNotebook.Tab", padding=[10, 5], font=("Segoe UI", 10))
        
    def _create_menu(self):
        """Create the application menu bar"""
        self.menu_bar = tk.Menu(self.root)
        
        # File menu
        file_menu = tk.Menu(self.menu_bar, tearoff=0)
        file_menu.add_command(label="Open FASTA File(s)", command=self._load_fasta_files, accelerator="Ctrl+O")
        file_menu.add_command(label="Save Results", command=self._save_results, accelerator="Ctrl+S")
        file_menu.add_command(label="Export Alignment as Image", command=self._export_alignment_image)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.root.quit, accelerator="Alt+F4")
        self.menu_bar.add_cascade(label="File", menu=file_menu)
        
        # Edit menu
        edit_menu = tk.Menu(self.menu_bar, tearoff=0)
        edit_menu.add_command(label="Clear Input", command=self._clear_input)
        edit_menu.add_command(label="Clear Results", command=self._clear_results)
        edit_menu.add_separator()
        edit_menu.add_command(label="Settings", command=self._show_settings)
        self.menu_bar.add_cascade(label="Edit", menu=edit_menu)
        
        # View menu
        view_menu = tk.Menu(self.menu_bar, tearoff=0)
        
        # Submenu for themes
        theme_menu = tk.Menu(view_menu, tearoff=0)
        theme_menu.add_radiobutton(label="Light Theme", variable=self.theme_mode, 
                                  value="light", command=self._apply_theme)
        theme_menu.add_radiobutton(label="Dark Theme", variable=self.theme_mode, 
                                  value="dark", command=self._apply_theme)
        view_menu.add_cascade(label="Theme", menu=theme_menu)
        
        # Submenu for font size
        font_menu = tk.Menu(view_menu, tearoff=0)
        for size in [8, 10, 12, 14, 16]:
            font_menu.add_command(label=f"{size} pt", 
                                 command=lambda s=size: self._set_font_size(s))
        view_menu.add_cascade(label="Font Size", menu=font_menu)
        
        self.menu_bar.add_cascade(label="View", menu=view_menu)
        
        # Help menu
        help_menu = tk.Menu(self.menu_bar, tearoff=0)
        help_menu.add_command(label="Documentation", command=self._show_documentation)
        help_menu.add_command(label="About", command=self._show_about)
        self.menu_bar.add_cascade(label="Help", menu=help_menu)
        
        self.root.config(menu=self.menu_bar)
        
        # Keyboard shortcuts
        self.root.bind("<Control-o>", lambda e: self._load_fasta_files())
        self.root.bind("<Control-s>", lambda e: self._save_results())
    
    def _create_ui(self):
        """Create the main UI components"""
        # Main container
        self.main_frame = ttk.Frame(self.root, padding="10")
        self.main_frame.grid(row=0, column=0, sticky="nsew")
        
        # Make rows and columns expandable
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        self.main_frame.columnconfigure(0, weight=1)
        
        # Create UI sections
        self._create_input_section()
        self._create_parameters_section()
        self._create_output_section()
        self._create_status_bar()
    
    def _create_input_section(self):
        """Create the sequence input section"""
        # Input frame
        input_frame = ttk.LabelFrame(self.main_frame, text="Sequence Input", padding="5")
        input_frame.grid(row=0, column=0, sticky="nsew", pady=(0, 10))
        input_frame.columnconfigure(0, weight=1)
        input_frame.rowconfigure(1, weight=1)
        
        # File path display
        self.file_path_var = tk.StringVar(value="No files loaded")
        ttk.Label(input_frame, textvariable=self.file_path_var).grid(row=0, column=0, sticky="w", pady=(0, 5))
        
        # Sequence input box with scrollbars
        self.sequence_input = tk.Text(input_frame, height=6, wrap="none", font=("Courier New", 11))
        self.sequence_input.grid(row=1, column=0, sticky="nsew")
        
        # Scrollbars for the input box
        y_scroll = ttk.Scrollbar(input_frame, orient="vertical", command=self.sequence_input.yview)
        y_scroll.grid(row=1, column=1, sticky="ns")
        x_scroll = ttk.Scrollbar(input_frame, orient="horizontal", command=self.sequence_input.xview)
        x_scroll.grid(row=2, column=0, sticky="ew")
        self.sequence_input.configure(yscrollcommand=y_scroll.set, xscrollcommand=x_scroll.set)
        
        # Button frame
        button_frame = ttk.Frame(input_frame, padding="5")
        button_frame.grid(row=3, column=0, columnspan=2, sticky="ew", pady=5)
        button_frame.columnconfigure(0, weight=1)
        button_frame.columnconfigure(1, weight=1)
        
        # Buttons
        ttk.Button(button_frame, text="Load FASTA File(s)", 
                  command=self._load_fasta_files, style="Secondary.TButton").grid(row=0, column=0, sticky="ew", padx=5)
        ttk.Button(button_frame, text="Clear Input", 
                  command=self._clear_input, style="Secondary.TButton").grid(row=0, column=1, sticky="ew", padx=5)
        
        # Add the input frame to the main grid
        self.main_frame.rowconfigure(0, weight=1)
    
    def _create_parameters_section(self):
        """Create the scoring parameters section"""
        # Parameters frame
        param_frame = ttk.LabelFrame(self.main_frame, text="Alignment Parameters", padding="5")
        param_frame.grid(row=1, column=0, sticky="ew", pady=(0, 10))
        
        # Grid configuration
        for i in range(3):
            param_frame.columnconfigure(i, weight=1)
        
        # Sequence type selection
        ttk.Label(param_frame, text="Sequence Type:").grid(row=0, column=0, sticky="e", pady=5, padx=5)
        type_dropdown = ttk.Combobox(param_frame, textvariable=self.seq_type_var, 
                                    values=["", "dna", "rna", "protein"], state="readonly", width=10)
        type_dropdown.grid(row=0, column=1, sticky="w", pady=5)
        ttk.Label(param_frame, text="(Leave empty for auto-detection)").grid(row=0, column=2, sticky="w")
        
        # Scoring parameters
        scoring_frame = ttk.Frame(param_frame)
        scoring_frame.grid(row=1, column=0, columnspan=3, sticky="ew", pady=5)
        
        for i in range(3):
            scoring_frame.columnconfigure(i * 2 + 1, weight=1)
        
        # Match score
        ttk.Label(scoring_frame, text="Match:").grid(row=0, column=0, sticky="e", padx=10)
        ttk.Entry(scoring_frame, textvariable=self.match_var, width=5).grid(row=0, column=1, sticky="w")
        
        # Mismatch score
        ttk.Label(scoring_frame, text="Mismatch:").grid(row=0, column=2, sticky="e", padx=10)
        ttk.Entry(scoring_frame, textvariable=self.mismatch_var, width=5).grid(row=0, column=3, sticky="w")
        
        # Gap score
        ttk.Label(scoring_frame, text="Gap:").grid(row=0, column=4, sticky="e", padx=10)
        ttk.Entry(scoring_frame, textvariable=self.gap_var, width=5).grid(row=0, column=5, sticky="w")
        
        # Run button
        ttk.Button(param_frame, text="RUN ALIGNMENT", 
                  command=self._run_alignment, style="Primary.TButton").grid(
                      row=2, column=0, columnspan=3, sticky="ew", pady=10, padx=50)
    
    def _create_output_section(self):
        """Create the output section with tabbed interface"""
        # Output frame
        output_frame = ttk.LabelFrame(self.main_frame, text="Results", padding="5")
        output_frame.grid(row=2, column=0, sticky="nsew", pady=(0, 5))
        output_frame.columnconfigure(0, weight=1)
        output_frame.rowconfigure(0, weight=1)
        
        # Notebook (tabbed interface)
        self.output_notebook = ttk.Notebook(output_frame)
        self.output_notebook.grid(row=0, column=0, sticky="nsew")
        
        # Text result tab
        text_frame = ttk.Frame(self.output_notebook, padding="5")
        text_frame.columnconfigure(0, weight=1)
        text_frame.rowconfigure(0, weight=1)
        
        self.result_text = tk.Text(text_frame, height=10, width=90, wrap="none", font=("Courier New", 11))
        self.result_text.grid(row=0, column=0, sticky="nsew")
        
        # Scrollbars for text results
        y_scroll = ttk.Scrollbar(text_frame, orient="vertical", command=self.result_text.yview)
        y_scroll.grid(row=0, column=1, sticky="ns")
        x_scroll = ttk.Scrollbar(text_frame, orient="horizontal", command=self.result_text.xview)
        x_scroll.grid(row=1, column=0, sticky="ew")
        self.result_text.configure(yscrollcommand=y_scroll.set, xscrollcommand=x_scroll.set)
        
        # Visual alignment tab
        visual_frame = ttk.Frame(self.output_notebook, padding="5")
        visual_frame.columnconfigure(0, weight=1)
        visual_frame.rowconfigure(0, weight=1)
        
        # Canvas for alignment visualization
        self.alignment_canvas = tk.Canvas(visual_frame, bg="white", highlightthickness=0)
        self.alignment_canvas.grid(row=0, column=0, sticky="nsew")
        
        # Scrollbars for alignment visualization
        y_scroll = ttk.Scrollbar(visual_frame, orient="vertical", command=self.alignment_canvas.yview)
        y_scroll.grid(row=0, column=1, sticky="ns")
        x_scroll = ttk.Scrollbar(visual_frame, orient="horizontal", command=self.alignment_canvas.xview)
        x_scroll.grid(row=1, column=0, sticky="ew")
        self.alignment_canvas.configure(yscrollcommand=y_scroll.set, xscrollcommand=x_scroll.set)
        
        # Statistics tab
        stats_frame = ttk.Frame(self.output_notebook, padding="5")
        stats_frame.columnconfigure(0, weight=1)
        stats_frame.rowconfigure(0, weight=1)
        
        self.stats_text = tk.Text(stats_frame, height=10, width=90, wrap="word", font=("Segoe UI", 11))
        self.stats_text.grid(row=0, column=0, sticky="nsew")
        
        # Scrollbars for stats
        y_scroll = ttk.Scrollbar(stats_frame, orient="vertical", command=self.stats_text.yview)
        y_scroll.grid(row=0, column=1, sticky="ns")
        self.stats_text.configure(yscrollcommand=y_scroll.set)
        
        # Add the frames to the notebook
        self.output_notebook.add(text_frame, text="Text Alignment")
        self.output_notebook.add(visual_frame, text="Visual Alignment")
        self.output_notebook.add(stats_frame, text="Statistics")
        
        # Button frame for actions
        button_frame = ttk.Frame(output_frame)
        button_frame.grid(row=1, column=0, sticky="ew", pady=(5, 0))
        
        # Buttons for saving and exporting
        ttk.Button(button_frame, text="Save Results", 
                  command=self._save_results, style="Secondary.TButton").pack(side="left", padx=5)
        ttk.Button(button_frame, text="Export Alignment", 
                  command=self._export_alignment_image, style="Secondary.TButton").pack(side="left", padx=5)
        ttk.Button(button_frame, text="Clear Results", 
                  command=self._clear_results, style="Secondary.TButton").pack(side="left", padx=5)
        
        # Make the output section expandable
        self.main_frame.rowconfigure(2, weight=3)
    
    def _create_status_bar(self):
        """Create the status bar at the bottom of the window"""
        status_bar = ttk.Label(self.main_frame, textvariable=self.status_var, 
                              relief="sunken", anchor="w", style="StatusBar.TLabel")
        status_bar.grid(row=3, column=0, sticky="ew")
    
    def _on_window_resize(self, event=None):
        """Handle window resize events to update font sizes"""
        # Only handle events from the main window
        if event and event.widget != self.root:
            return
            
        # Update font sizes based on window height
        min_size = 10
        window_height = self.root.winfo_height()
        
        if window_height < 600:
            base_size = min_size
        elif window_height < 800:
            base_size = min_size + 1
        elif window_height < 1000:
            base_size = min_size + 2
        else:
            base_size = min_size + 3
            
        # Update alignment font size
        self.alignment_font_size = base_size
        
        # If we have alignment data, redraw it
        if hasattr(self, "last_msa") and self.last_msa and hasattr(self, "last_seq_objects") and self.last_seq_objects:
            self._draw_alignment_blocks(self.last_msa, self.last_seq_objects)
    
    def _load_fasta_files(self):
        """Open file dialog and load FASTA files"""
        filepaths = filedialog.askopenfilenames(
            title="Select FASTA Files",
            filetypes=[("FASTA files", "*.fasta *.fa"), ("All files", "*")]
        )
        
        if not filepaths:
            return
            
        try:
            # Clear previous input
            self.sequence_input.delete("1.0", tk.END)
            self.file_paths = list(filepaths)
            
            # Update file path display
            if len(filepaths) == 1:
                self.file_path_var.set(f"File: {os.path.basename(filepaths[0])}")
            else:
                self.file_path_var.set(f"Files: {len(filepaths)} files selected")
            
            # Load sequences
            all_sequences = []
            for filepath in filepaths:
                sequences = parse_fasta_file(filepath)
                all_sequences.extend(sequences)
                
            for seq in all_sequences:
                self.sequence_input.insert(tk.END, seq.strip() + "\n")
                
            self.status_var.set(f"Loaded {len(all_sequences)} sequences from {len(filepaths)} file(s)")
        except Exception as e:
            messagebox.showerror("File Load Error", str(e))
            self.status_var.set("Error loading file(s)")
    
    def _save_results(self, event=None):
        """Save the results to a text file"""
        filepath = filedialog.asksaveasfilename(
            title="Save Results",
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*")]
        )
        
        if not filepath:
            return
            
        try:
            content = self.result_text.get("1.0", tk.END).strip()
            with open(filepath, "w") as f:
                f.write(content)
            messagebox.showinfo("Saved", f"Results saved to {filepath}")
            self.status_var.set(f"Results saved to: {os.path.basename(filepath)}")
        except Exception as e:
            messagebox.showerror("Save Error", str(e))
            self.status_var.set("Error saving results")
    
    def _export_alignment_image(self):
        """Export the alignment visualization as an image"""
        if not hasattr(self, "last_msa") or not self.last_msa:
            messagebox.showinfo("No Alignment", "Please run an alignment first.")
            return
            
        filepath = filedialog.asksaveasfilename(
            title="Export Alignment Image",
            defaultextension=".png",
            filetypes=[("PNG images", "*.png"), ("All files", "*")]
        )
        
        if not filepath:
            return
            
        try:
            # Get the canvas dimensions
            x0, y0, x1, y1 = self.alignment_canvas.bbox("all")
            
            # Create a temporary postscript file
            ps_file = filepath + ".ps"
            self.alignment_canvas.postscript(file=ps_file, colormode="color", 
                                          width=x1-x0, height=y1-y0)
            
            # For a more robust solution, you'd use a library like PIL to convert PS to PNG
            # For simplicity, we'll just inform the user
            messagebox.showinfo("Export", 
                             f"Alignment exported as PostScript file: {ps_file}\n"
                             "You can convert it to PNG using an image editor.")
            self.status_var.set(f"Alignment exported to: {os.path.basename(ps_file)}")
        except Exception as e:
            messagebox.showerror("Export Error", str(e))
            self.status_var.set("Error exporting alignment")
    
    def _clear_input(self):
        """Clear the input text area"""
        self.sequence_input.delete("1.0", tk.END)
        self.file_path_var.set("No files loaded")
        self.file_paths = []
        self.status_var.set("Input cleared")
    
    def _clear_results(self):
        """Clear all results"""
        self.result_text.delete("1.0", tk.END)
        self.stats_text.delete("1.0", tk.END)
        self.alignment_canvas.delete("all")
        self.last_msa = None
        self.last_seq_objects = None
        self.status_var.set("Results cleared")
    
    def _show_settings(self):
        """Show settings dialog"""
        settings_window = tk.Toplevel(self.root)
        settings_window.title("Settings")
        settings_window.geometry("400x300")
        settings_window.transient(self.root)
        settings_window.grab_set()
        
        # Center on parent
        settings_window.update_idletasks()
        x = self.root.winfo_x() + (self.root.winfo_width() - settings_window.winfo_width()) // 2
        y = self.root.winfo_y() + (self.root.winfo_height() - settings_window.winfo_height()) // 2
        settings_window.geometry(f"+{x}+{y}")
        
        # Create settings UI
        settings_frame = ttk.Frame(settings_window, padding="10")
        settings_frame.pack(fill=tk.BOTH, expand=True)
        
        # Theme settings
        ttk.Label(settings_frame, text="Theme:").grid(row=0, column=0, sticky="w", pady=5)
        theme_frame = ttk.Frame(settings_frame)
        theme_frame.grid(row=0, column=1, sticky="w", pady=5)
        ttk.Radiobutton(theme_frame, text="Light", variable=self.theme_mode, 
                       value="light").pack(side=tk.LEFT, padx=5)
        ttk.Radiobutton(theme_frame, text="Dark", variable=self.theme_mode, 
                       value="dark").pack(side=tk.LEFT, padx=5)
        
        # Default scoring settings
        ttk.Label(settings_frame, text="Default Scoring:").grid(row=1, column=0, sticky="w", pady=5)
        scoring_frame = ttk.Frame(settings_frame)
        scoring_frame.grid(row=1, column=1, sticky="w", pady=5)
        
        # Match
        ttk.Label(scoring_frame, text="Match:").grid(row=0, column=0, padx=5)
        match_entry = ttk.Entry(scoring_frame, width=5)
        match_entry.grid(row=0, column=1, padx=5)
        match_entry.insert(0, self.DEFAULT_SCORING["match"])
        
        # Mismatch
        ttk.Label(scoring_frame, text="Mismatch:").grid(row=0, column=2, padx=5)
        mismatch_entry = ttk.Entry(scoring_frame, width=5)
        mismatch_entry.grid(row=0, column=3, padx=5)
        mismatch_entry.insert(0, self.DEFAULT_SCORING["mismatch"])
        
        # Gap
        ttk.Label(scoring_frame, text="Gap:").grid(row=0, column=4, padx=5)
        gap_entry = ttk.Entry(scoring_frame, width=5)
        gap_entry.grid(row=0, column=5, padx=5)
        gap_entry.insert(0, self.DEFAULT_SCORING["gap"])
        
        # Font settings
        ttk.Label(settings_frame, text="Base Font Size:").grid(row=2, column=0, sticky="w", pady=5)
        font_frame = ttk.Frame(settings_frame)
        font_frame.grid(row=2, column=1, sticky="w", pady=5)
        font_size = ttk.Spinbox(font_frame, from_=8, to=16, width=5)
        font_size.grid(row=0, column=0)
        font_size.set(self.alignment_font_size)
        
        # Buttons
        button_frame = ttk.Frame(settings_frame)
        button_frame.grid(row=3, column=0, columnspan=2, pady=20)
        
        ttk.Button(button_frame, text="Apply", command=lambda: self._apply_settings(
            match_entry.get(), mismatch_entry.get(), gap_entry.get(), font_size.get()
        )).pack(side=tk.LEFT, padx=10)
        ttk.Button(button_frame, text="Cancel", command=settings_window.destroy).pack(side=tk.LEFT, padx=10)
    
    def _apply_settings(self, match, mismatch, gap, font_size):
        """Apply new settings"""
        try:
            # Update default scoring
            self.DEFAULT_SCORING["match"] = int(match)
            self.DEFAULT_SCORING["mismatch"] = int(mismatch)
            self.DEFAULT_SCORING["gap"] = int(gap)
            
            # Update current scoring if not custom
            self.match_var.set(str(self.DEFAULT_SCORING["match"]))
            self.mismatch_var.set(str(self.DEFAULT_SCORING["mismatch"]))
            self.gap_var.set(str(self.DEFAULT_SCORING["gap"]))
            
            # Update font size
            self.alignment_font_size = int(font_size)
            
            # Apply theme
            self._apply_theme()
            
            # Save settings
            self._save_settings()
            
            # Update visualization if necessary
            if self.last_msa and self.last_seq_objects:
                self._draw_alignment_blocks(self.last_msa, self.last_seq_objects)
                
            self.status_var.set("Settings applied")
        except Exception as e:
            messagebox.showerror("Settings Error", str(e))
    
    def _apply_theme(self):
        """Apply the selected theme to the application"""
        if self.theme_mode.get() == "dark":
            # Dark theme colors
            self.COLORS.update({
                "background": "#121212",  # Dark background
                "surface": "#1e1e1e",     # Dark surface
                "primary": "#bb86fc",     # Purple
                "secondary": "#03dac6",   # Teal
                "on_primary": "#000000",  # Black text on primary
                "on_secondary": "#000000", # Black text on secondary
            })
            
            # Configure styles for dark theme
            self.style.configure("TFrame", background=self.COLORS["background"])
            self.style.configure("TLabel", background=self.COLORS["background"], foreground="white")
            self.style.configure("TButton", background=self.COLORS["surface"])
            self.style.configure("Primary.TButton", background=self.COLORS["primary"])
            self.style.configure("Secondary.TButton", background=self.COLORS["secondary"])
            
            # Configure text widgets for dark mode
            self.sequence_input.config(bg="#1e1e1e", fg="white", insertbackground="white")
            self.result_text.config(bg="#1e1e1e", fg="white", insertbackground="white")
            self.stats_text.config(bg="#1e1e1e", fg="white", insertbackground="white")
            self.alignment_canvas.config(bg="#1e1e1e")
            
            # Status bar dark mode
            self.style.configure("StatusBar.TLabel", background="#1e1e1e", foreground="white")
        else:
            # Light theme colors
            self.COLORS.update({
                "background": "#f5f5f5",  # Light background
                "surface": "#ffffff",     # White surface
                "primary": "#6200ee",     # Purple
                "secondary": "#03dac6",   # Teal
                "on_primary": "#ffffff",  # White text on primary
                "on_secondary": "#000000", # Black text on secondary
            })
            
            # Configure styles for light theme
            self.style.configure("TFrame", background=self.COLORS["background"])
            self.style.configure("TLabel", background=self.COLORS["background"], foreground="black")
            self.style.configure("TButton", background=self.COLORS["surface"])
            self.style.configure("Primary.TButton", background=self.COLORS["primary"])
            self.style.configure("Secondary.TButton", background=self.COLORS["secondary"])
            
            # Configure text widgets for light mode
            self.sequence_input.config(bg="white", fg="black", insertbackground="black")
            self.result_text.config(bg="white", fg="black", insertbackground="black")
            self.stats_text.config(bg="white", fg="black", insertbackground="black")
            self.alignment_canvas.config(bg="white")
            
            # Status bar light mode
            self.style.configure("StatusBar.TLabel", background="#f0f0f0", foreground="black")
    
    def _set_font_size(self, size):
        """Set the font size for all text elements"""
        self.alignment_font_size = size
        
        # Update text widgets
        self.sequence_input.config(font=("Courier New", size))
        self.result_text.config(font=("Courier New", size))
        self.stats_text.config(font=("Segoe UI", size))
        
        # Redraw alignment if available
        if self.last_msa and self.last_seq_objects:
            self._draw_alignment_blocks(self.last_msa, self.last_seq_objects)
            
        self.status_var.set(f"Font size changed to {size}pt")
    
    def _show_documentation(self):
        """Show documentation in a new window"""
        doc_window = tk.Toplevel(self.root)
        doc_window.title("MSA Documentation")
        doc_window.geometry("600x400")
        doc_window.minsize(400, 300)
        
        # Center on parent
        doc_window.update_idletasks()
        x = self.root.winfo_x() + (self.root.winfo_width() - doc_window.winfo_width()) // 2
        y = self.root.winfo_y() + (self.root.winfo_height() - doc_window.winfo_height()) // 2
        doc_window.geometry(f"+{x}+{y}")
        
        # Create documentation text
        doc_frame = ttk.Frame(doc_window, padding="10")
        doc_frame.pack(fill=tk.BOTH, expand=True)
        
        doc_text = tk.Text(doc_frame, wrap="word", padx=10, pady=10)
        doc_text.pack(fill=tk.BOTH, expand=True)
        
        # Documentation content
        doc_content = """# Multiple Sequence Alignment Tool

## Overview
This application performs multiple sequence alignment using the Center Star Method. 
It aligns DNA, RNA, or protein sequences to identify similarities and differences.

## Features
- Load sequences from FASTA files or paste them directly
- Auto-detection of sequence type (DNA, RNA, protein)
- Customizable scoring parameters
- Visual alignment display
- Statistics calculation
- Export results in various formats

## How to Use
1. Enter sequences in the input box or load from FASTA files
2. Set alignment parameters or use defaults
3. Click "RUN ALIGNMENT" to perform the alignment
4. View results in the tabs below
5. Save or export results as needed

## Sequence Types
- DNA: Nucleotide sequences with A, T, G, C
- RNA: Nucleotide sequences with A, U, G, C
- Protein: Amino acid sequences

## Scoring Parameters
- Match: Score for matching characters
- Mismatch: Penalty for mismatched characters
- Gap: Penalty for gaps inserted in the alignment

## Output Explanation
- Text Alignment: Shows the aligned sequences
- Visual Alignment: Color-coded visualization of alignment
- Statistics: Information about alignment quality

## Tips
- Use larger match scores relative to gap penalties for closely related sequences
- For distantly related sequences, use smaller gap penalties
"""
        doc_text.insert("1.0", doc_content)
        doc_text.config(state="disabled")
        
        # Add scrollbar
        scrollbar = ttk.Scrollbar(doc_frame, orient="vertical", command=doc_text.yview)
        scrollbar.pack(side="right", fill="y")
        doc_text.configure(yscrollcommand=scrollbar.set)
    
    def _show_about(self):
        """Show about dialog"""
        messagebox.showinfo(
            "About MSA Tool",
            "Multiple Sequence Alignment Tool\n"
            "Version 2.0\n\n"
            "This application performs multiple sequence alignment using the Center Star Method.\n"
            "It can align DNA, RNA, or protein sequences to identify similarities and differences."
        )
    
    def _save_settings(self):
        """Save settings to a JSON file"""
        settings = {
            "theme": self.theme_mode.get(),
            "font_size": self.alignment_font_size,
            "default_scoring": self.DEFAULT_SCORING
        }
        
        try:
            settings_dir = Path.home() / ".msa_app"
            settings_dir.mkdir(exist_ok=True)
            
            with open(settings_dir / "settings.json", "w") as f:
                json.dump(settings, f)
        except Exception as e:
            print(f"Error saving settings: {e}")
    
    def _load_settings(self):
        """Load settings from a JSON file"""
        try:
            settings_path = Path.home() / ".msa_app" / "settings.json"
            
            if settings_path.exists():
                with open(settings_path, "r") as f:
                    settings = json.load(f)
                
                # Apply loaded settings
                if "theme" in settings:
                    self.theme_mode.set(settings["theme"])
                
                if "font_size" in settings:
                    self.alignment_font_size = settings["font_size"]
                
                if "default_scoring" in settings:
                    self.DEFAULT_SCORING = settings["default_scoring"]
                    self.match_var.set(str(self.DEFAULT_SCORING["match"]))
                    self.mismatch_var.set(str(self.DEFAULT_SCORING["mismatch"]))
                    self.gap_var.set(str(self.DEFAULT_SCORING["gap"]))
                
                # Apply theme
                self._apply_theme()
        except Exception as e:
            print(f"Error loading settings: {e}")
    
    def _run_alignment(self):
        """Run the alignment process"""
        # Update status
        self.status_var.set("Running alignment...")
        self.root.update_idletasks()
        
        # Run in a separate thread to keep UI responsive
        threading.Thread(target=self._perform_alignment, daemon=True).start()
    
    def _perform_alignment(self):
        """Perform the alignment calculation"""
        try:
            # Get input sequences
            raw_input = self.sequence_input.get("1.0", tk.END).strip()
            if not raw_input:
                raise ValueError("Please enter at least one valid sequence.")
            
            # Normalize and validate sequences
            sequences = normalize_sequences(raw_input.replace("\n", " "))
            if not sequences:
                raise ValueError("No valid sequences found.")
            
            # Determine sequence type
            seq_type = self.seq_type_var.get()
            if not seq_type:
                sequence_type = detect_sequence_type(sequences)
                # Update the dropdown with detected type
                self.root.after(0, lambda: self.seq_type_var.set(sequence_type))
            else:
                sequence_type = seq_type
            
            # Validate sequences based on type
            validate_sequences(sequences, sequence_type)
            
            # Create sequence objects
            seq_objects = [
                Sequence(f"seq{i+1}", seq, alphabet=sequence_type)
                for i, seq in enumerate(sequences)
            ]
            
            # Create scoring scheme
            try:
                scoring = ScoringScheme(
                    match=int(self.match_var.get()),
                    mismatch=int(self.mismatch_var.get()),
                    gap=int(self.gap_var.get())
                )
            except ValueError:
                raise ValueError("Scoring parameters must be integers.")
            
            # Build score matrix
            score_matrix = build_pairwise_score_matrix(seq_objects, scoring)
            
            # Find center sequence
            distance_matrix = convert_scores_to_distances(score_matrix)
            center_index = find_center_sequence(distance_matrix)
            
            # Align sequences to center
            aligned = align_all_to_center(seq_objects, center_index, scoring)
            
            # Create MSA
            final_msa = merge_alignments_to_msa(aligned, center_index)
            
            # Compute statistics
            stats = compute_msa_statistics(final_msa)
            
            # Store results for later
            self.last_msa = final_msa
            self.last_seq_objects = seq_objects
            
            # Update UI in the main thread
            self.root.after(0, lambda: self._update_results(final_msa, seq_objects, stats, scoring, sequence_type, center_index))
            
        except Exception as e:
            # Show error in the main thread
            self.root.after(0, lambda: self._show_error(str(e)))
    
    def _update_results(self, final_msa, seq_objects, stats, scoring, sequence_type, center_index):
        """Update the UI with alignment results"""
        # Clear previous results
        self.result_text.delete("1.0", tk.END)
        self.stats_text.delete("1.0", tk.END)
        
        # Display basic information
        self.result_text.insert(tk.END, f"Scoring: Match={scoring.match}, Mismatch={scoring.mismatch}, Gap={scoring.gap}\n")
        self.result_text.insert(tk.END, f"Sequence Type: {sequence_type.upper()}\n")
        self.result_text.insert(tk.END, f"Center Sequence: {seq_objects[center_index].id}\n\n")
        
        # Display alignment
        self.result_text.insert(tk.END, "Alignment:\n")
        for i, aligned_seq in enumerate(final_msa):
            self.result_text.insert(tk.END, f"{seq_objects[i].id}: {aligned_seq}\n")
        
        # Display statistics
        self.stats_text.insert(tk.END, "Alignment Statistics:\n\n")
        
        for k, v in stats.items():
            # Format the key for better readability
            key_name = k.replace('_', ' ').title()
            self.stats_text.insert(tk.END, f"{key_name}: {v}\n")
        
        # Add some explanation for the statistics
        self.stats_text.insert(tk.END, "\nExplanation:\n\n")
        self.stats_text.insert(tk.END, "• Identity: Percentage of positions with identical residues in all sequences\n")
        self.stats_text.insert(tk.END, "• Conservation: Measure of biochemical similarity at each position\n")
        self.stats_text.insert(tk.END, "• Gaps: Percentage of gap positions in the alignment\n")
        self.stats_text.insert(tk.END, "• Length: Length of the aligned sequences (including gaps)\n")
        self.stats_text.insert(tk.END, "• Sum of Pairs Score: Overall quality measure for the alignment\n")
        
        # Draw the visual alignment
        self._draw_alignment_blocks(final_msa, seq_objects)
        
        # Switch to the text alignment tab
        self.output_notebook.select(0)
        
        # Update status
        self.status_var.set(f"Alignment complete: {len(seq_objects)} sequences aligned")
    
    def _draw_alignment_blocks(self, final_msa, seq_objects):
        """Draw the alignment visualization on the canvas"""
        # Clear canvas
        self.alignment_canvas.delete("all")
        
        # Get the number of sequences and alignment length
        num_sequences = len(final_msa)
        alignment_length = len(final_msa[0]) if final_msa else 0
        
        if not alignment_length:
            return
        
        # Configure cell dimensions based on font size
        cell_width = max(16, self.alignment_font_size * 1.5)
        cell_height = max(16, self.alignment_font_size * 2)
        padding = 10
        id_width = 80  # Width for sequence IDs
        
        # Configure font
        font = ("Courier New", self.alignment_font_size)
        
        # Calculate canvas dimensions
        canvas_width = id_width + alignment_length * cell_width + padding * 2
        canvas_height = num_sequences * cell_height + padding * 2
        
        # Configure scrolling region
        self.alignment_canvas.configure(scrollregion=(0, 0, canvas_width, canvas_height))
        
        # Draw column numbers (every 10 positions)
        for i in range(0, alignment_length, 10):
            col_num = i + 1  # 1-based indexing for display
            x = id_width + i * cell_width + cell_width // 2
            self.alignment_canvas.create_text(x, padding // 2, 
                                           text=str(col_num), font=font, anchor="s")
        
        # Calculate column conservation for coloring
        column_conservation = []
        for col_idx in range(alignment_length):
            column = [seq[col_idx] for seq in final_msa]
            
            # Skip all-gap columns
            if all(c == '-' for c in column):
                column_conservation.append("gap")
                continue
            
            # Check if all non-gap characters are the same
            non_gaps = [c for c in column if c != '-']
            if all(c == non_gaps[0] for c in non_gaps):
                column_conservation.append("match")
            else:
                column_conservation.append("mismatch")
        
        # Draw the alignment blocks
        for row_idx, aligned_seq in enumerate(final_msa):
            # Y position for this sequence
            y = padding + row_idx * cell_height
            
            # Draw sequence ID
            seq_id = seq_objects[row_idx].id if len(seq_objects[row_idx].id) <= 10 else seq_objects[row_idx].id[:10]
            self.alignment_canvas.create_text(padding, y + cell_height // 2, 
                                           text=seq_id, font=font, anchor="w")
            
            # Draw sequence characters
            for col_idx, char in enumerate(aligned_seq):
                # X position for this character
                x = id_width + col_idx * cell_width
                
                # Determine cell color
                conservation = column_conservation[col_idx]
                if char == '-':
                    color = self.COLORS["gap"]
                elif conservation == "match":
                    color = self.COLORS["match"]
                else:
                    color = self.COLORS["mismatch"]
                
                # Draw cell background
                self.alignment_canvas.create_rectangle(
                    x, y, x + cell_width, y + cell_height,
                    fill=color, outline="gray"
                )
                
                # Draw character
                self.alignment_canvas.create_text(
                    x + cell_width // 2, y + cell_height // 2,
                    text=char, font=font
                )
        
        # Draw a ruler at the bottom
        ruler_y = padding + num_sequences * cell_height + 5
        self.alignment_canvas.create_line(
            id_width, ruler_y, id_width + alignment_length * cell_width, ruler_y,
            width=2
        )
        
        for i in range(0, alignment_length + 1, 10):
            tick_x = id_width + i * cell_width
            self.alignment_canvas.create_line(
                tick_x, ruler_y, tick_x, ruler_y + 5,
                width=2
            )
            if i > 0:  # Don't show 0
                self.alignment_canvas.create_text(
                    tick_x, ruler_y + 15,
                    text=str(i), font=("Arial", 8)
                )
    
    def _show_error(self, message):
        """Show error message and update status"""
        messagebox.showerror("Error", message)
        self.status_var.set("Error: " + message[:50] + ("..." if len(message) > 50 else ""))

def main():
    """Main function to start the application"""
    root = tk.Tk()
    app = MSAApplication(root)
    root.mainloop()


if __name__ == "__main__":
    main()