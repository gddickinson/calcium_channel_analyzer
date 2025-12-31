#!/usr/bin/env python3
"""
Advanced Calcium Channel Bioinformatics Analysis Platform
==========================================================

A comprehensive tool for identifying, analyzing, and characterizing calcium channels
across species using state-of-the-art bioinformatics approaches.

Features:
- Multi-species sequence analysis with modern BLAST alternatives
- Domain detection (EF-hand, SPRY, CNBD, GLR, TPC, etc.)
- Phylogenetic analysis and tree construction
- AlphaFold3 structure prediction integration
- Transmembrane domain prediction
- Interactive visualizations
- Comprehensive reporting

Author: George Dickinson
Institution: UC Irvine, Department of Neurobiology & Behavior
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import threading
import queue
from pathlib import Path
import json
import re
from datetime import datetime
from typing import List, Dict, Tuple, Optional
import webbrowser

# Bioinformatics libraries
from Bio import SeqIO, Entrez, AlignIO, Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align import PairwiseAligner

# Data analysis
import pandas as pd
import numpy as np

# Visualization
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import matplotlib.patches as mpatches
import seaborn as sns

# Set styles
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 10


class CalciumChannelAnalyzer:
    """Main application class for calcium channel analysis."""

    # Calcium channel conserved domain patterns (based on literature and Pfam)
    DOMAIN_PATTERNS = {
        'EF_hand': {
            'pattern': r'D.{1,3}D.{1,3}[DNS].{1,3}[EDNST]',
            'description': 'EF-hand calcium-binding motif',
            'min_length': 12,
            'channels': ['All voltage-gated', 'RyR', 'IP3R', 'TPC']
        },
        'CNBD': {
            'pattern': r'[RK].{1,5}[LIVMF].{1,5}G.{1,5}[LIVMF].{1,5}[LIVMF]',
            'description': 'Cyclic nucleotide-binding domain',
            'min_length': 20,
            'channels': ['CNGC']
        },
        'Pore_loop': {
            'pattern': r'[TGSA][LIVMF][LIVMF]G[LIVMF][GAS]',
            'description': 'Ion channel selectivity filter',
            'min_length': 6,
            'channels': ['All Ca2+ channels']
        },
        'SPRY': {
            'pattern': r'[LIVMF].{2,3}[LIVMF].{2,3}[LIVMF].{10,15}[LIVMF].{2,3}[LIVMF]',
            'description': 'SPRY domain (found in RyR)',
            'min_length': 25,
            'channels': ['RyR']
        },
        'IP3_binding': {
            'pattern': r'[KR].{1,3}[LIVMF].{1,3}[LIVMF].{1,3}[KR]',
            'description': 'IP3 receptor binding domain',
            'min_length': 15,
            'channels': ['IP3R']
        },
        'Voltage_sensor': {
            'pattern': r'[RK]{2,4}.{3,5}[RK]{2,4}',
            'description': 'Voltage-sensing domain S4',
            'min_length': 10,
            'channels': ['Voltage-gated Ca2+ channels', 'TPC']
        },
        'TM_helix': {
            'pattern': r'[AILVMFGW]{7,}',
            'description': 'Potential transmembrane helix',
            'min_length': 18,
            'channels': ['All membrane channels']
        },
        'Ca_selectivity': {
            'pattern': r'[DE].{2,4}[DE].{2,4}[DE]',
            'description': 'Calcium selectivity filter',
            'min_length': 8,
            'channels': ['L-type', 'P/Q-type', 'N-type']
        },
        'Calmodulin_binding': {
            'pattern': r'[ILVF]{3,5}[RK]{1,2}[ILVF]{3,5}[RK]{1,2}',
            'description': 'Calmodulin-binding IQ motif',
            'min_length': 15,
            'channels': ['CNGC', 'L-type Ca2+']
        },
        'GLR_ligand': {
            'pattern': r'[LIVMF].{1,3}[DE].{1,3}[LIVMF].{1,3}[DE]',
            'description': 'Glutamate receptor-like ligand binding',
            'min_length': 12,
            'channels': ['GLR']
        }
    }

    # Known calcium channel families
    CHANNEL_FAMILIES = {
        'Voltage-gated': {
            'Cav1': ['Cav1.1', 'Cav1.2', 'Cav1.3', 'Cav1.4'],  # L-type
            'Cav2': ['Cav2.1', 'Cav2.2', 'Cav2.3'],  # P/Q, N, R-type
            'Cav3': ['Cav3.1', 'Cav3.2', 'Cav3.3']   # T-type
        },
        'Ligand-gated': {
            'RyR': ['RyR1', 'RyR2', 'RyR3'],
            'IP3R': ['IP3R1', 'IP3R2', 'IP3R3']
        },
        'Plant-specific': {
            'CNGC': ['CNGC1-20'],  # Arabidopsis has 20 CNGCs
            'GLR': ['GLR1.1-3.7'],  # Glutamate receptor-like
            'TPC': ['TPC1'],         # Two-pore channel
            'MSL': ['MSL1-10'],      # Mechanosensitive
            'MCA': ['MCA1-2']        # Mid1-complementing activity
        },
        'Other': {
            'TRP': ['TRPC', 'TRPV', 'TRPM', 'TRPA', 'TRPML', 'TRPP'],
            'P2X': ['P2X1-7'],
            'Orai': ['Orai1-3']
        }
    }

    def __init__(self, root):
        self.root = root
        self.root.title("Calcium Channel Bioinformatics Analyzer v2.0")
        self.root.geometry("1400x900")

        # Set email for Entrez (required for NCBI API)
        Entrez.email = "george.dickinson@uci.edu"

        # Data storage
        self.sequences = []
        self.alignment = None
        self.tree = None
        self.blast_results = []
        self.domain_results = []
        self.current_figure = None

        # Thread management
        self.task_queue = queue.Queue()
        self.result_queue = queue.Queue()

        # Setup UI
        self.setup_ui()

        # Start result checker
        self.check_results()

    def setup_ui(self):
        """Create the user interface."""
        # Create menu bar
        menubar = tk.Menu(self.root)
        self.root.config(menu=menubar)

        # File menu
        file_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="File", menu=file_menu)
        file_menu.add_command(label="Load FASTA", command=self.load_fasta)
        file_menu.add_command(label="Load from GenBank", command=self.load_genbank)
        file_menu.add_separator()
        file_menu.add_command(label="Export Results", command=self.export_results)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.root.quit)

        # Analysis menu
        analysis_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Analysis", menu=analysis_menu)
        analysis_menu.add_command(label="BLAST Search", command=self.run_blast_dialog)
        analysis_menu.add_command(label="Detect Domains", command=self.detect_domains)
        analysis_menu.add_command(label="Multiple Alignment", command=self.run_alignment)
        analysis_menu.add_command(label="Phylogenetic Tree", command=self.build_tree)
        analysis_menu.add_command(label="Classify Channels", command=self.classify_channels)

        # Tools menu
        tools_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Tools", menu=tools_menu)
        tools_menu.add_command(label="AlphaFold Prediction", command=self.alphafold_predict)
        tools_menu.add_command(label="Transmembrane Prediction", command=self.predict_tm_domains)
        tools_menu.add_command(label="Hydropathy Plot", command=self.plot_hydropathy)

        # Help menu
        help_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Help", menu=help_menu)
        help_menu.add_command(label="User Guide", command=self.show_help)
        help_menu.add_command(label="About", command=self.show_about)

        # Main container with notebook
        main_container = ttk.Frame(self.root)
        main_container.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Create notebook (tabs)
        self.notebook = ttk.Notebook(main_container)
        self.notebook.pack(fill=tk.BOTH, expand=True)

        # Tab 1: Sequence Manager
        self.setup_sequence_tab()

        # Tab 2: Domain Analysis
        self.setup_domain_tab()

        # Tab 3: Phylogenetic Analysis
        self.setup_phylo_tab()

        # Tab 4: Structure Prediction
        self.setup_structure_tab()

        # Tab 5: Results Summary
        self.setup_results_tab()

        # Status bar
        self.status_bar = ttk.Label(self.root, text="Ready", relief=tk.SUNKEN)
        self.status_bar.pack(side=tk.BOTTOM, fill=tk.X)

        # Progress bar
        self.progress = ttk.Progressbar(self.root, mode='indeterminate')

    def setup_sequence_tab(self):
        """Setup sequence management tab."""
        seq_frame = ttk.Frame(self.notebook)
        self.notebook.add(seq_frame, text="Sequences")

        # Top panel with controls
        control_frame = ttk.LabelFrame(seq_frame, text="Sequence Management")
        control_frame.pack(fill=tk.X, padx=5, pady=5)

        ttk.Button(control_frame, text="Load FASTA",
                  command=self.load_fasta).grid(row=0, column=0, padx=5, pady=5)
        ttk.Button(control_frame, text="Load GenBank",
                  command=self.load_genbank).grid(row=0, column=1, padx=5, pady=5)
        ttk.Button(control_frame, text="Add Manual Entry",
                  command=self.add_manual_sequence).grid(row=0, column=2, padx=5, pady=5)
        ttk.Button(control_frame, text="Clear All",
                  command=self.clear_sequences).grid(row=0, column=3, padx=5, pady=5)

        # Sequence info label
        self.seq_count_label = ttk.Label(control_frame, text="Sequences loaded: 0")
        self.seq_count_label.grid(row=0, column=4, padx=20, pady=5)

        # Sequence list with scrollbar
        list_frame = ttk.Frame(seq_frame)
        list_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Treeview for sequences
        columns = ('ID', 'Description', 'Length', 'Family')
        self.seq_tree = ttk.Treeview(list_frame, columns=columns, show='tree headings')

        # Configure columns
        self.seq_tree.heading('#0', text='#')
        self.seq_tree.heading('ID', text='Accession/ID')
        self.seq_tree.heading('Description', text='Description')
        self.seq_tree.heading('Length', text='Length (aa)')
        self.seq_tree.heading('Family', text='Predicted Family')

        self.seq_tree.column('#0', width=50)
        self.seq_tree.column('ID', width=150)
        self.seq_tree.column('Description', width=400)
        self.seq_tree.column('Length', width=100)
        self.seq_tree.column('Family', width=150)

        # Scrollbars
        vsb = ttk.Scrollbar(list_frame, orient="vertical", command=self.seq_tree.yview)
        hsb = ttk.Scrollbar(list_frame, orient="horizontal", command=self.seq_tree.xview)
        self.seq_tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)

        self.seq_tree.grid(row=0, column=0, sticky='nsew')
        vsb.grid(row=0, column=1, sticky='ns')
        hsb.grid(row=1, column=0, sticky='ew')

        list_frame.grid_rowconfigure(0, weight=1)
        list_frame.grid_columnconfigure(0, weight=1)

        # Sequence details panel
        detail_frame = ttk.LabelFrame(seq_frame, text="Sequence Details")
        detail_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        self.seq_details = scrolledtext.ScrolledText(detail_frame, height=8, wrap=tk.WORD)
        self.seq_details.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Bind selection event
        self.seq_tree.bind('<<TreeviewSelect>>', self.on_sequence_select)

    def setup_domain_tab(self):
        """Setup domain analysis tab."""
        domain_frame = ttk.Frame(self.notebook)
        self.notebook.add(domain_frame, text="Domain Analysis")

        # Control panel
        control_frame = ttk.LabelFrame(domain_frame, text="Domain Detection")
        control_frame.pack(fill=tk.X, padx=5, pady=5)

        ttk.Label(control_frame, text="Select domains to search:").grid(row=0, column=0,
                                                                         padx=5, pady=5, sticky='w')

        # Checkboxes for domain types
        self.domain_vars = {}
        col = 0
        row = 1
        for domain_name, domain_info in self.DOMAIN_PATTERNS.items():
            var = tk.BooleanVar(value=True)
            self.domain_vars[domain_name] = var
            ttk.Checkbutton(control_frame, text=f"{domain_name} - {domain_info['description']}",
                           variable=var).grid(row=row, column=col, sticky='w', padx=5)
            row += 1
            if row > 5:
                row = 1
                col += 1

        ttk.Button(control_frame, text="Run Domain Detection",
                  command=self.detect_domains).grid(row=7, column=0, columnspan=3, pady=10)

        # Results display
        results_frame = ttk.LabelFrame(domain_frame, text="Domain Detection Results")
        results_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Treeview for domain results
        columns = ('Sequence', 'Domain', 'Position', 'Sequence_Match', 'Score')
        self.domain_tree = ttk.Treeview(results_frame, columns=columns, show='headings')

        for col in columns:
            self.domain_tree.heading(col, text=col.replace('_', ' '))
            self.domain_tree.column(col, width=150)

        vsb = ttk.Scrollbar(results_frame, orient="vertical", command=self.domain_tree.yview)
        self.domain_tree.configure(yscrollcommand=vsb.set)

        self.domain_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        vsb.pack(side=tk.RIGHT, fill=tk.Y)

        # Visualization frame
        viz_frame = ttk.LabelFrame(domain_frame, text="Domain Map Visualization")
        viz_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        ttk.Button(viz_frame, text="Generate Domain Map",
                  command=self.plot_domain_map).pack(pady=5)

        self.domain_canvas_frame = ttk.Frame(viz_frame)
        self.domain_canvas_frame.pack(fill=tk.BOTH, expand=True)

    def setup_phylo_tab(self):
        """Setup phylogenetic analysis tab."""
        phylo_frame = ttk.Frame(self.notebook)
        self.notebook.add(phylo_frame, text="Phylogenetics")

        # Control panel
        control_frame = ttk.LabelFrame(phylo_frame, text="Phylogenetic Analysis")
        control_frame.pack(fill=tk.X, padx=5, pady=5)

        ttk.Label(control_frame, text="Algorithm:").grid(row=0, column=0, padx=5, pady=5)
        self.phylo_method = ttk.Combobox(control_frame,
                                         values=['UPGMA', 'Neighbor Joining'],
                                         state='readonly')
        self.phylo_method.set('Neighbor Joining')
        self.phylo_method.grid(row=0, column=1, padx=5, pady=5)

        ttk.Button(control_frame, text="Build Phylogenetic Tree",
                  command=self.build_tree).grid(row=0, column=2, padx=5, pady=5)

        ttk.Button(control_frame, text="Export Tree (Newick)",
                  command=self.export_tree).grid(row=0, column=3, padx=5, pady=5)

        # Tree visualization
        self.phylo_canvas_frame = ttk.Frame(phylo_frame)
        self.phylo_canvas_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

    def setup_structure_tab(self):
        """Setup structure prediction tab."""
        struct_frame = ttk.Frame(self.notebook)
        self.notebook.add(struct_frame, text="Structure Prediction")

        # Info panel
        info_frame = ttk.LabelFrame(struct_frame, text="AlphaFold3 Integration")
        info_frame.pack(fill=tk.X, padx=5, pady=5)

        info_text = """
        AlphaFold3 Protein Structure Prediction:

        • AlphaFold3 can predict calcium channel structures with unprecedented accuracy
        • Predicts protein-protein interactions, protein-ligand complexes
        • Useful for understanding channel gating mechanisms and drug binding
        • Submit sequences to AlphaFold Server for non-commercial research

        Note: This tool will help you prepare sequences and link to AlphaFold Server.
        Results can be imported back for analysis.
        """
        ttk.Label(info_frame, text=info_text, justify=tk.LEFT).pack(padx=10, pady=10)

        # Sequence selection
        select_frame = ttk.LabelFrame(struct_frame, text="Select Sequence for Prediction")
        select_frame.pack(fill=tk.X, padx=5, pady=5)

        ttk.Label(select_frame, text="Sequence:").grid(row=0, column=0, padx=5, pady=5)
        self.struct_seq_combo = ttk.Combobox(select_frame, width=50)
        self.struct_seq_combo.grid(row=0, column=1, padx=5, pady=5)

        ttk.Button(select_frame, text="Open AlphaFold Server",
                  command=self.open_alphafold).grid(row=0, column=2, padx=5, pady=5)

        ttk.Button(select_frame, text="Export FASTA for AlphaFold",
                  command=self.export_for_alphafold).grid(row=0, column=3, padx=5, pady=5)

        # Transmembrane prediction
        tm_frame = ttk.LabelFrame(struct_frame, text="Transmembrane Domain Prediction")
        tm_frame.pack(fill=tk.X, padx=5, pady=5)

        ttk.Button(tm_frame, text="Predict TM Domains (Kyte-Doolittle)",
                  command=self.predict_tm_domains).pack(pady=5)

        # Results display
        self.struct_results = scrolledtext.ScrolledText(struct_frame, height=15, wrap=tk.WORD)
        self.struct_results.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

    def setup_results_tab(self):
        """Setup results summary tab."""
        results_frame = ttk.Frame(self.notebook)
        self.notebook.add(results_frame, text="Results Summary")

        # Control panel
        control_frame = ttk.LabelFrame(results_frame, text="Export Options")
        control_frame.pack(fill=tk.X, padx=5, pady=5)

        ttk.Button(control_frame, text="Export Full Report (HTML)",
                  command=self.export_html_report).grid(row=0, column=0, padx=5, pady=5)
        ttk.Button(control_frame, text="Export to Excel",
                  command=self.export_excel).grid(row=0, column=1, padx=5, pady=5)
        ttk.Button(control_frame, text="Export Sequences (FASTA)",
                  command=self.export_fasta).grid(row=0, column=2, padx=5, pady=5)

        # Summary text
        self.summary_text = scrolledtext.ScrolledText(results_frame, wrap=tk.WORD)
        self.summary_text.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

    # ==================== Core Functions ====================

    def update_status(self, message):
        """Update status bar."""
        self.status_bar.config(text=message)
        self.root.update_idletasks()

    def show_progress(self):
        """Show progress bar."""
        self.progress.pack(side=tk.BOTTOM, fill=tk.X, before=self.status_bar)
        self.progress.start()

    def hide_progress(self):
        """Hide progress bar."""
        self.progress.stop()
        self.progress.pack_forget()

    def load_fasta(self):
        """Load sequences from FASTA file."""
        filename = filedialog.askopenfilename(
            title="Select FASTA file",
            filetypes=[("FASTA files", "*.fasta *.fa *.faa"), ("All files", "*.*")]
        )

        if filename:
            try:
                self.update_status(f"Loading {filename}...")
                records = list(SeqIO.parse(filename, "fasta"))
                self.sequences.extend(records)
                self.update_sequence_list()
                self.update_status(f"Loaded {len(records)} sequences from {Path(filename).name}")
                messagebox.showinfo("Success", f"Loaded {len(records)} sequences")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load file: {str(e)}")
                self.update_status("Error loading file")

    def load_genbank(self):
        """Load sequence from GenBank by accession."""
        dialog = tk.Toplevel(self.root)
        dialog.title("Load from GenBank")
        dialog.geometry("400x150")

        ttk.Label(dialog, text="Enter GenBank Accession(s):").pack(pady=10)

        entry = ttk.Entry(dialog, width=50)
        entry.pack(pady=5)
        entry.focus()

        ttk.Label(dialog, text="(separate multiple accessions with commas)").pack()

        def fetch():
            accessions = [acc.strip() for acc in entry.get().split(',')]
            dialog.destroy()

            self.show_progress()
            self.update_status(f"Fetching {len(accessions)} sequence(s) from GenBank...")

            def fetch_thread():
                results = []
                for acc in accessions:
                    try:
                        handle = Entrez.efetch(db="protein", id=acc, rettype="fasta", retmode="text")
                        record = SeqIO.read(handle, "fasta")
                        results.append(record)
                        handle.close()
                    except Exception as e:
                        results.append(f"Error fetching {acc}: {str(e)}")
                self.result_queue.put(('genbank', results))

            thread = threading.Thread(target=fetch_thread, daemon=True)
            thread.start()

        ttk.Button(dialog, text="Fetch", command=fetch).pack(pady=10)
        ttk.Button(dialog, text="Cancel", command=dialog.destroy).pack()

    def add_manual_sequence(self):
        """Add sequence manually."""
        dialog = tk.Toplevel(self.root)
        dialog.title("Add Manual Sequence")
        dialog.geometry("600x400")

        ttk.Label(dialog, text="Sequence ID:").pack(pady=5)
        id_entry = ttk.Entry(dialog, width=50)
        id_entry.pack(pady=5)

        ttk.Label(dialog, text="Description:").pack(pady=5)
        desc_entry = ttk.Entry(dialog, width=50)
        desc_entry.pack(pady=5)

        ttk.Label(dialog, text="Sequence (amino acids):").pack(pady=5)
        seq_text = scrolledtext.ScrolledText(dialog, width=60, height=10)
        seq_text.pack(pady=5, padx=10)

        def add():
            seq_id = id_entry.get().strip()
            description = desc_entry.get().strip()
            sequence = seq_text.get("1.0", tk.END).strip().upper()
            # Remove whitespace and numbers
            sequence = ''.join(c for c in sequence if c.isalpha())

            if not seq_id or not sequence:
                messagebox.showwarning("Warning", "ID and sequence are required")
                return

            try:
                record = SeqRecord(Seq(sequence), id=seq_id, description=description)
                self.sequences.append(record)
                self.update_sequence_list()
                dialog.destroy()
                messagebox.showinfo("Success", "Sequence added")
            except Exception as e:
                messagebox.showerror("Error", f"Invalid sequence: {str(e)}")

        ttk.Button(dialog, text="Add", command=add).pack(pady=5)
        ttk.Button(dialog, text="Cancel", command=dialog.destroy).pack()

    def clear_sequences(self):
        """Clear all loaded sequences."""
        if messagebox.askyesno("Confirm", "Clear all sequences?"):
            self.sequences = []
            self.update_sequence_list()
            self.update_status("Sequences cleared")

    def update_sequence_list(self):
        """Update the sequence treeview."""
        # Clear existing items
        for item in self.seq_tree.get_children():
            self.seq_tree.delete(item)

        # Add sequences
        for idx, record in enumerate(self.sequences, 1):
            # Try to predict family based on description
            family = self.predict_family_from_description(record.description)

            self.seq_tree.insert('', tk.END, text=str(idx),
                               values=(record.id,
                                      record.description[:60],
                                      len(record.seq),
                                      family))

        # Update count
        self.seq_count_label.config(text=f"Sequences loaded: {len(self.sequences)}")

        # Update structure prediction combo
        seq_list = [f"{r.id} - {r.description[:40]}" for r in self.sequences]
        self.struct_seq_combo['values'] = seq_list
        if seq_list:
            self.struct_seq_combo.current(0)

    def predict_family_from_description(self, description: str) -> str:
        """Predict channel family from description."""
        desc_lower = description.lower()

        # Check for known patterns
        if any(term in desc_lower for term in ['cav1', 'l-type', 'ltcc']):
            return 'Cav1 (L-type)'
        elif any(term in desc_lower for term in ['cav2', 'p/q-type', 'n-type']):
            return 'Cav2 (P/Q/N-type)'
        elif any(term in desc_lower for term in ['cav3', 't-type']):
            return 'Cav3 (T-type)'
        elif 'ryanodine' in desc_lower or 'ryr' in desc_lower:
            return 'RyR'
        elif 'ip3' in desc_lower or 'insp3' in desc_lower:
            return 'IP3R'
        elif 'cngc' in desc_lower or 'cyclic nucleotide' in desc_lower:
            return 'CNGC'
        elif 'glr' in desc_lower or 'glutamate receptor' in desc_lower:
            return 'GLR'
        elif 'tpc' in desc_lower or 'two-pore' in desc_lower:
            return 'TPC'
        elif 'trp' in desc_lower:
            return 'TRP'
        elif 'orai' in desc_lower:
            return 'Orai'

        return 'Unknown'

    def on_sequence_select(self, event):
        """Handle sequence selection."""
        selection = self.seq_tree.selection()
        if selection:
            item = self.seq_tree.item(selection[0])
            idx = int(item['text']) - 1

            if 0 <= idx < len(self.sequences):
                record = self.sequences[idx]
                details = f"ID: {record.id}\n"
                details += f"Description: {record.description}\n"
                details += f"Length: {len(record.seq)} amino acids\n\n"
                details += f"Sequence:\n{record.seq}\n\n"

                # Calculate composition
                seq_str = str(record.seq)
                details += "Amino Acid Composition:\n"
                for aa in sorted(set(seq_str)):
                    count = seq_str.count(aa)
                    percent = (count / len(seq_str)) * 100
                    details += f"{aa}: {count} ({percent:.1f}%)\n"

                self.seq_details.delete(1.0, tk.END)
                self.seq_details.insert(1.0, details)

    def detect_domains(self):
        """Detect conserved domains in sequences."""
        if not self.sequences:
            messagebox.showwarning("Warning", "No sequences loaded")
            return

        # Get selected domains
        selected_domains = [name for name, var in self.domain_vars.items() if var.get()]

        if not selected_domains:
            messagebox.showwarning("Warning", "No domains selected")
            return

        self.show_progress()
        self.update_status("Detecting domains...")

        def detect_thread():
            results = []
            for record in self.sequences:
                seq_str = str(record.seq)

                for domain_name in selected_domains:
                    domain_info = self.DOMAIN_PATTERNS[domain_name]
                    pattern = domain_info['pattern']

                    # Find all matches
                    for match in re.finditer(pattern, seq_str):
                        position = f"{match.start()+1}-{match.end()}"
                        matched_seq = match.group()
                        score = self.calculate_domain_score(matched_seq, domain_name)

                        results.append({
                            'sequence': record.id,
                            'domain': domain_name,
                            'position': position,
                            'match': matched_seq,
                            'score': f"{score:.2f}",
                            'description': domain_info['description']
                        })

            self.result_queue.put(('domains', results))

        thread = threading.Thread(target=detect_thread, daemon=True)
        thread.start()

    def calculate_domain_score(self, sequence: str, domain_name: str) -> float:
        """Calculate a simple conservation score for domain match."""
        # This is a simplified scoring - could be enhanced with position-specific scoring
        domain_info = self.DOMAIN_PATTERNS[domain_name]
        min_length = domain_info['min_length']

        # Length score
        length_score = min(len(sequence) / min_length, 1.0) * 50

        # Hydrophobicity score for TM domains
        if domain_name == 'TM_helix':
            hydrophobic = 'AILVMFGW'
            hydro_count = sum(1 for aa in sequence if aa in hydrophobic)
            hydro_score = (hydro_count / len(sequence)) * 50
        else:
            hydro_score = 25

        return length_score + hydro_score

    def plot_domain_map(self):
        """Visualize domain architecture."""
        if not self.domain_results:
            messagebox.showinfo("Info", "No domain results to plot. Run domain detection first.")
            return

        # Clear previous plot
        for widget in self.domain_canvas_frame.winfo_children():
            widget.destroy()

        # Group results by sequence
        seq_domains = {}
        for result in self.domain_results:
            seq_id = result['sequence']
            if seq_id not in seq_domains:
                seq_domains[seq_id] = []
            seq_domains[seq_id].append(result)

        # Create figure
        n_seqs = len(seq_domains)
        fig, axes = plt.subplots(n_seqs, 1, figsize=(14, max(6, n_seqs * 1.5)))
        if n_seqs == 1:
            axes = [axes]

        # Define colors for different domain types
        domain_colors = {
            'EF_hand': '#FF6B6B',
            'CNBD': '#4ECDC4',
            'Pore_loop': '#45B7D1',
            'SPRY': '#FFA07A',
            'IP3_binding': '#98D8C8',
            'Voltage_sensor': '#F7DC6F',
            'TM_helix': '#BB8FCE',
            'Ca_selectivity': '#85C1E2',
            'Calmodulin_binding': '#F8B88B',
            'GLR_ligand': '#FAD7A0'
        }

        for idx, (seq_id, domains) in enumerate(seq_domains.items()):
            ax = axes[idx]

            # Get sequence length
            seq_record = next(r for r in self.sequences if r.id == seq_id)
            seq_length = len(seq_record.seq)

            # Draw sequence as a line
            ax.plot([0, seq_length], [0.5, 0.5], 'k-', linewidth=2, zorder=1)

            # Draw domains
            for domain in domains:
                start, end = map(int, domain['position'].split('-'))
                domain_type = domain['domain']
                color = domain_colors.get(domain_type, '#CCCCCC')

                # Draw rectangle for domain
                rect = mpatches.Rectangle((start, 0.3), end-start, 0.4,
                                         facecolor=color, edgecolor='black',
                                         linewidth=1, alpha=0.7, zorder=2)
                ax.add_patch(rect)

                # Add label if space allows
                if (end - start) > seq_length * 0.05:
                    ax.text((start + end) / 2, 0.5, domain_type,
                           ha='center', va='center', fontsize=8, zorder=3)

            ax.set_xlim(-seq_length*0.02, seq_length*1.02)
            ax.set_ylim(0, 1)
            ax.set_xlabel('Position (aa)')
            ax.set_ylabel(seq_id, fontsize=10)
            ax.set_yticks([])
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(False)

        # Add legend
        legend_elements = [mpatches.Patch(facecolor=color, edgecolor='black',
                                         label=domain, alpha=0.7)
                          for domain, color in domain_colors.items()]
        fig.legend(handles=legend_elements, loc='upper right',
                  bbox_to_anchor=(0.98, 0.98), ncol=2, fontsize=9)

        plt.tight_layout()

        # Embed in tkinter
        canvas = FigureCanvasTkAgg(fig, master=self.domain_canvas_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Add toolbar
        toolbar = NavigationToolbar2Tk(canvas, self.domain_canvas_frame)
        toolbar.update()

        self.current_figure = fig
        self.update_status("Domain map generated")

    def run_alignment(self):
        """Run multiple sequence alignment."""
        if len(self.sequences) < 2:
            messagebox.showwarning("Warning", "Need at least 2 sequences for alignment")
            return

        messagebox.showinfo("Info",
                           "For production use, external alignment tools like MUSCLE, "
                           "MAFFT, or Clustal Omega are recommended. "
                           "This demo uses pairwise alignment for illustration.")

        self.show_progress()
        self.update_status("Performing alignment...")

        def align_thread():
            # Simple pairwise alignment demo
            # In production, would use Bio.Align.Applications for MUSCLE/MAFFT
            results = []

            # Create aligner
            aligner = PairwiseAligner()
            aligner.mode = 'global'
            aligner.match_score = 1
            aligner.mismatch_score = 0
            aligner.open_gap_score = -0.5
            aligner.extend_gap_score = -0.1

            for i in range(len(self.sequences)):
                for j in range(i+1, len(self.sequences)):
                    # Limit sequence length for demo performance
                    seq1 = str(self.sequences[i].seq)[:500]
                    seq2 = str(self.sequences[j].seq)[:500]

                    # Perform alignment
                    alignments = aligner.align(seq1, seq2)

                    if alignments:
                        best_alignment = alignments[0]
                        results.append({
                            'seq1': self.sequences[i].id,
                            'seq2': self.sequences[j].id,
                            'score': best_alignment.score,
                            'alignment': best_alignment
                        })

            self.result_queue.put(('alignment', results))

        thread = threading.Thread(target=align_thread, daemon=True)
        thread.start()

    def build_tree(self):
        """Build phylogenetic tree."""
        if len(self.sequences) < 3:
            messagebox.showwarning("Warning", "Need at least 3 sequences for tree construction")
            return

        self.show_progress()
        self.update_status("Building phylogenetic tree...")

        def tree_thread():
            try:
                # Create simple alignment (in production, use proper MSA)
                sequences = [str(s.seq) for s in self.sequences]
                ids = [s.id for s in self.sequences]

                # Ensure equal length for demo (pad shorter sequences)
                max_len = max(len(s) for s in sequences)
                aligned_seqs = [s + '-' * (max_len - len(s)) for s in sequences]

                # Create alignment object
                records = [SeqRecord(Seq(seq), id=id_)
                          for seq, id_ in zip(aligned_seqs, ids)]
                alignment = MultipleSeqAlignment(records)

                # Calculate distance matrix
                calculator = DistanceCalculator('identity')
                dm = calculator.get_distance(alignment)

                # Construct tree
                constructor = DistanceTreeConstructor(calculator)
                method = self.phylo_method.get()

                if method == 'UPGMA':
                    tree = constructor.upgma(dm)
                else:
                    tree = constructor.nj(dm)

                self.result_queue.put(('tree', (tree, alignment)))

            except Exception as e:
                self.result_queue.put(('error', f"Tree construction failed: {str(e)}"))

        thread = threading.Thread(target=tree_thread, daemon=True)
        thread.start()

    def plot_tree(self, tree):
        """Plot phylogenetic tree."""
        # Clear previous plot
        for widget in self.phylo_canvas_frame.winfo_children():
            widget.destroy()

        fig, ax = plt.subplots(figsize=(12, 8))
        Phylo.draw(tree, axes=ax, do_show=False)
        plt.title("Phylogenetic Tree of Calcium Channel Sequences", fontsize=14, fontweight='bold')
        plt.tight_layout()

        # Embed in tkinter
        canvas = FigureCanvasTkAgg(fig, master=self.phylo_canvas_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Add toolbar
        toolbar = NavigationToolbar2Tk(canvas, self.phylo_canvas_frame)
        toolbar.update()

        self.current_figure = fig
        self.update_status("Phylogenetic tree generated")

    def predict_tm_domains(self):
        """Predict transmembrane domains using Kyte-Doolittle hydropathy."""
        if not self.sequences:
            messagebox.showwarning("Warning", "No sequences loaded")
            return

        # Select sequence
        if not hasattr(self, 'struct_seq_combo') or not self.struct_seq_combo.get():
            if self.sequences:
                idx = 0
            else:
                return
        else:
            idx = self.struct_seq_combo.current()

        record = self.sequences[idx]
        seq = str(record.seq)

        # Kyte-Doolittle hydropathy scale
        kd = {'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8, 'G': -0.4,
              'H': -3.2, 'I': 4.5, 'K': -3.9, 'L': 3.8, 'M': 1.9, 'N': -3.5,
              'P': -1.6, 'Q': -3.5, 'R': -4.5, 'S': -0.8, 'T': -0.7, 'V': 4.2,
              'W': -0.9, 'Y': -1.3}

        # Calculate hydropathy with window
        window = 19
        hydropathy = []
        positions = []

        for i in range(len(seq) - window + 1):
            window_seq = seq[i:i+window]
            score = sum(kd.get(aa, 0) for aa in window_seq) / window
            hydropathy.append(score)
            positions.append(i + window//2)

        # Find TM regions (hydropathy > 1.6 typically indicates TM)
        tm_threshold = 1.6
        in_tm = False
        tm_regions = []
        start = 0

        for i, score in enumerate(hydropathy):
            if score > tm_threshold and not in_tm:
                start = positions[i]
                in_tm = True
            elif score <= tm_threshold and in_tm:
                tm_regions.append((start, positions[i]))
                in_tm = False

        if in_tm:
            tm_regions.append((start, positions[-1]))

        # Display results
        results = f"Transmembrane Domain Prediction for {record.id}\n"
        results += f"{'='*60}\n\n"
        results += f"Sequence length: {len(seq)} aa\n"
        results += f"Predicted TM domains: {len(tm_regions)}\n\n"

        for i, (start, end) in enumerate(tm_regions, 1):
            results += f"TM{i}: {start}-{end} ({end-start} aa)\n"
            results += f"Sequence: {seq[start:end]}\n\n"

        # Plot
        fig, ax = plt.subplots(figsize=(12, 6))
        ax.plot(positions, hydropathy, 'b-', linewidth=2)
        ax.axhline(y=tm_threshold, color='r', linestyle='--', label=f'TM threshold ({tm_threshold})')
        ax.axhline(y=0, color='gray', linestyle='-', alpha=0.3)
        ax.fill_between(positions, hydropathy, tm_threshold,
                        where=np.array(hydropathy) >= tm_threshold,
                        color='red', alpha=0.3, label='TM regions')

        ax.set_xlabel('Position', fontsize=12)
        ax.set_ylabel('Hydropathy Score', fontsize=12)
        ax.set_title(f'Kyte-Doolittle Hydropathy Plot - {record.id}', fontsize=14, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.tight_layout()

        # Show plot
        plot_window = tk.Toplevel(self.root)
        plot_window.title("Hydropathy Plot")
        plot_window.geometry("900x600")

        canvas = FigureCanvasTkAgg(fig, master=plot_window)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2Tk(canvas, plot_window)
        toolbar.update()

        # Show text results
        self.struct_results.delete(1.0, tk.END)
        self.struct_results.insert(1.0, results)

        self.update_status("TM prediction complete")

    def plot_hydropathy(self):
        """Alias for predict_tm_domains."""
        self.predict_tm_domains()

    def classify_channels(self):
        """Classify loaded sequences into calcium channel families."""
        if not self.sequences:
            messagebox.showwarning("Warning", "No sequences loaded")
            return

        self.show_progress()
        self.update_status("Classifying channel sequences...")

        def classify_thread():
            classifications = []

            for record in self.sequences:
                seq_str = str(record.seq)
                desc = record.description.lower()

                # Score for each family
                family_scores = {}

                # Check voltage-gated markers
                has_voltage_sensor = bool(re.search(self.DOMAIN_PATTERNS['Voltage_sensor']['pattern'], seq_str))
                has_ca_selectivity = bool(re.search(self.DOMAIN_PATTERNS['Ca_selectivity']['pattern'], seq_str))

                # Check ligand-gated markers
                has_cnbd = bool(re.search(self.DOMAIN_PATTERNS['CNBD']['pattern'], seq_str))
                has_spry = bool(re.search(self.DOMAIN_PATTERNS['SPRY']['pattern'], seq_str))
                has_ip3_binding = bool(re.search(self.DOMAIN_PATTERNS['IP3_binding']['pattern'], seq_str))

                # Classification logic
                classification = "Unknown"
                confidence = "Low"

                # Check description first
                if any(term in desc for term in ['cav1', 'l-type']):
                    classification = "Voltage-gated (Cav1/L-type)"
                    confidence = "High"
                elif any(term in desc for term in ['cav2', 'p/q', 'n-type']):
                    classification = "Voltage-gated (Cav2/P-Q-N-type)"
                    confidence = "High"
                elif any(term in desc for term in ['cav3', 't-type']):
                    classification = "Voltage-gated (Cav3/T-type)"
                    confidence = "High"
                elif 'ryanodine' in desc or 'ryr' in desc:
                    classification = "Ligand-gated (RyR)"
                    confidence = "High"
                elif 'ip3' in desc:
                    classification = "Ligand-gated (IP3R)"
                    confidence = "High"
                elif 'cngc' in desc:
                    classification = "Plant CNGC"
                    confidence = "High"
                elif 'tpc' in desc or 'two-pore' in desc:
                    classification = "Two-pore channel (TPC)"
                    confidence = "High"

                # If not in description, use domains
                elif has_voltage_sensor and has_ca_selectivity:
                    classification = "Voltage-gated Ca2+ channel"
                    confidence = "Medium"
                elif has_spry:
                    classification = "Ligand-gated (likely RyR)"
                    confidence = "Medium"
                elif has_ip3_binding:
                    classification = "Ligand-gated (likely IP3R)"
                    confidence = "Medium"
                elif has_cnbd:
                    classification = "CNGC-like channel"
                    confidence = "Medium"

                classifications.append({
                    'id': record.id,
                    'description': record.description[:60],
                    'classification': classification,
                    'confidence': confidence,
                    'length': len(seq_str),
                    'markers': {
                        'voltage_sensor': has_voltage_sensor,
                        'ca_selectivity': has_ca_selectivity,
                        'cnbd': has_cnbd,
                        'spry': has_spry,
                        'ip3_binding': has_ip3_binding
                    }
                })

            self.result_queue.put(('classification', classifications))

        thread = threading.Thread(target=classify_thread, daemon=True)
        thread.start()

    def run_blast_dialog(self):
        """Show BLAST search dialog."""
        if not self.sequences:
            messagebox.showwarning("Warning", "No sequences loaded")
            return

        dialog = tk.Toplevel(self.root)
        dialog.title("BLAST Search")
        dialog.geometry("500x300")

        info = """BLAST Search Options:

        Note: For production use, consider faster alternatives:
        • DIAMOND (100x faster, good sensitivity)
        • MMseqs2 (fast with low error rates)
        • BLAST+ (local installation recommended)

        This demo uses NCBI BLAST web service (slower).
        For large-scale analyses, use local BLAST databases.
        """

        ttk.Label(dialog, text=info, justify=tk.LEFT).pack(padx=10, pady=10)

        ttk.Label(dialog, text="Select sequence:").pack()
        seq_combo = ttk.Combobox(dialog, values=[s.id for s in self.sequences], state='readonly')
        seq_combo.pack(pady=5)
        if self.sequences:
            seq_combo.current(0)

        ttk.Label(dialog, text="Database:").pack()
        db_combo = ttk.Combobox(dialog, values=['nr', 'swissprot', 'pdb'], state='readonly')
        db_combo.set('nr')
        db_combo.pack(pady=5)

        ttk.Label(dialog, text="Max hits:").pack()
        hits_entry = ttk.Entry(dialog)
        hits_entry.insert(0, "10")
        hits_entry.pack(pady=5)

        def run_blast():
            idx = seq_combo.current()
            if idx < 0:
                messagebox.showwarning("Warning", "Select a sequence")
                return

            database = db_combo.get()
            max_hits = int(hits_entry.get())
            dialog.destroy()

            messagebox.showinfo("BLAST",
                              "BLAST search started. This may take several minutes.\n"
                              "Check status bar for progress.")

            self.show_progress()
            self.update_status("Running BLAST search...")

            def blast_thread():
                try:
                    record = self.sequences[idx]
                    result_handle = NCBIWWW.qblast("blastp", database, record.seq, hitlist_size=max_hits)
                    blast_records = list(NCBIXML.parse(result_handle))
                    self.result_queue.put(('blast', blast_records))
                except Exception as e:
                    self.result_queue.put(('error', f"BLAST failed: {str(e)}"))

            thread = threading.Thread(target=blast_thread, daemon=True)
            thread.start()

        ttk.Button(dialog, text="Run BLAST", command=run_blast).pack(pady=10)
        ttk.Button(dialog, text="Cancel", command=dialog.destroy).pack()

    def alphafold_predict(self):
        """Prepare sequence for AlphaFold prediction."""
        if not self.sequences:
            messagebox.showwarning("Warning", "No sequences loaded")
            return

        # Create dialog
        dialog = tk.Toplevel(self.root)
        dialog.title("AlphaFold Structure Prediction")
        dialog.geometry("600x400")

        info_text = """
        AlphaFold3 Structure Prediction Workflow:

        1. Select a sequence below
        2. Click 'Export FASTA' to save sequence file
        3. Visit AlphaFold Server (opens in browser)
        4. Upload your FASTA file
        5. Download predicted structure (PDB file)

        AlphaFold3 can predict:
        • Protein structure (near-atomic accuracy)
        • Protein-protein interactions
        • Protein-ligand complexes
        • Multi-domain assemblies

        Results typically available in minutes to hours.
        """

        ttk.Label(dialog, text=info_text, justify=tk.LEFT).pack(padx=10, pady=10)

        ttk.Label(dialog, text="Select Sequence:").pack(pady=5)
        seq_combo = ttk.Combobox(dialog,
                                 values=[f"{s.id} - {s.description[:40]}" for s in self.sequences],
                                 width=60,
                                 state='readonly')
        seq_combo.pack(pady=5)
        if self.sequences:
            seq_combo.current(0)

        button_frame = ttk.Frame(dialog)
        button_frame.pack(pady=20)

        def export_seq():
            idx = seq_combo.current()
            if idx < 0:
                messagebox.showwarning("Warning", "Please select a sequence")
                return

            record = self.sequences[idx]
            filename = filedialog.asksaveasfilename(
                defaultextension=".fasta",
                filetypes=[("FASTA files", "*.fasta"), ("All files", "*.*")],
                initialfile=f"{record.id}_alphafold.fasta"
            )

            if filename:
                with open(filename, 'w') as f:
                    SeqIO.write([record], f, "fasta")
                messagebox.showinfo("Success",
                                  f"Sequence exported to:\n{filename}\n\n"
                                  "Now opening AlphaFold Server...")
                webbrowser.open("https://alphafoldserver.com/")

        ttk.Button(button_frame, text="Export FASTA & Open AlphaFold Server",
                  command=export_seq).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Close",
                  command=dialog.destroy).pack(side=tk.LEFT, padx=5)

    def open_alphafold(self):
        """Open AlphaFold Server in browser."""
        webbrowser.open("https://alphafoldserver.com/")
        messagebox.showinfo("AlphaFold Server",
                           "AlphaFold Server opened in browser.\n"
                           "Use 'Export FASTA for AlphaFold' to prepare your sequence.")

    def export_for_alphafold(self):
        """Export selected sequence for AlphaFold."""
        if not self.sequences:
            messagebox.showwarning("Warning", "No sequences loaded")
            return

        idx = self.struct_seq_combo.current()
        if idx < 0:
            idx = 0

        record = self.sequences[idx]

        filename = filedialog.asksaveasfilename(
            defaultextension=".fasta",
            filetypes=[("FASTA files", "*.fasta"), ("All files", "*.*")],
            initialfile=f"{record.id}_alphafold.fasta"
        )

        if filename:
            with open(filename, 'w') as f:
                SeqIO.write([record], f, "fasta")
            messagebox.showinfo("Success", f"Sequence exported to {filename}\n"
                              "Upload this file to AlphaFold Server for structure prediction.")

    def export_tree(self):
        """Export phylogenetic tree in Newick format."""
        if not self.tree:
            messagebox.showwarning("Warning", "No tree generated. Build tree first.")
            return

        filename = filedialog.asksaveasfilename(
            defaultextension=".nwk",
            filetypes=[("Newick files", "*.nwk"), ("All files", "*.*")]
        )

        if filename:
            Phylo.write(self.tree, filename, "newick")
            messagebox.showinfo("Success", f"Tree exported to {filename}")

    def export_fasta(self):
        """Export sequences to FASTA file."""
        if not self.sequences:
            messagebox.showwarning("Warning", "No sequences to export")
            return

        filename = filedialog.asksaveasfilename(
            defaultextension=".fasta",
            filetypes=[("FASTA files", "*.fasta"), ("All files", "*.*")]
        )

        if filename:
            with open(filename, 'w') as f:
                SeqIO.write(self.sequences, f, "fasta")
            messagebox.showinfo("Success", f"Exported {len(self.sequences)} sequences")

    def export_excel(self):
        """Export results to Excel."""
        if not self.domain_results and not self.sequences:
            messagebox.showwarning("Warning", "No data to export")
            return

        filename = filedialog.asksaveasfilename(
            defaultextension=".xlsx",
            filetypes=[("Excel files", "*.xlsx"), ("All files", "*.*")]
        )

        if filename:
            try:
                with pd.ExcelWriter(filename, engine='openpyxl') as writer:
                    # Sequences sheet
                    if self.sequences:
                        seq_data = []
                        for record in self.sequences:
                            seq_data.append({
                                'ID': record.id,
                                'Description': record.description,
                                'Length': len(record.seq),
                                'Sequence': str(record.seq)
                            })
                        pd.DataFrame(seq_data).to_excel(writer, sheet_name='Sequences', index=False)

                    # Domains sheet
                    if self.domain_results:
                        pd.DataFrame(self.domain_results).to_excel(writer, sheet_name='Domains', index=False)

                messagebox.showinfo("Success", f"Data exported to {filename}")
            except Exception as e:
                messagebox.showerror("Error", f"Export failed: {str(e)}")

    def export_html_report(self):
        """Export comprehensive HTML report."""
        if not self.sequences:
            messagebox.showwarning("Warning", "No data to export")
            return

        filename = filedialog.asksaveasfilename(
            defaultextension=".html",
            filetypes=[("HTML files", "*.html"), ("All files", "*.*")]
        )

        if filename:
            html = self.generate_html_report()
            with open(filename, 'w') as f:
                f.write(html)
            messagebox.showinfo("Success", f"Report exported to {filename}")

    def generate_html_report(self) -> str:
        """Generate comprehensive HTML report."""
        html = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>Calcium Channel Analysis Report</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 20px; }}
                h1 {{ color: #2c3e50; }}
                h2 {{ color: #34495e; border-bottom: 2px solid #3498db; }}
                table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
                th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
                th {{ background-color: #3498db; color: white; }}
                .summary {{ background-color: #ecf0f1; padding: 15px; border-radius: 5px; }}
            </style>
        </head>
        <body>
            <h1>Calcium Channel Bioinformatics Analysis Report</h1>
            <p><strong>Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>

            <div class="summary">
                <h2>Summary</h2>
                <p><strong>Total Sequences:</strong> {len(self.sequences)}</p>
                <p><strong>Domains Detected:</strong> {len(self.domain_results)}</p>
            </div>

            <h2>Sequences</h2>
            <table>
                <tr>
                    <th>#</th>
                    <th>ID</th>
                    <th>Description</th>
                    <th>Length (aa)</th>
                </tr>
        """

        for idx, record in enumerate(self.sequences, 1):
            html += f"""
                <tr>
                    <td>{idx}</td>
                    <td>{record.id}</td>
                    <td>{record.description}</td>
                    <td>{len(record.seq)}</td>
                </tr>
            """

        html += """
            </table>
        """

        if self.domain_results:
            html += """
            <h2>Domain Analysis</h2>
            <table>
                <tr>
                    <th>Sequence</th>
                    <th>Domain</th>
                    <th>Position</th>
                    <th>Match</th>
                    <th>Score</th>
                </tr>
            """

            for result in self.domain_results:
                html += f"""
                <tr>
                    <td>{result['sequence']}</td>
                    <td>{result['domain']}</td>
                    <td>{result['position']}</td>
                    <td><code>{result['match']}</code></td>
                    <td>{result['score']}</td>
                </tr>
                """

            html += """
            </table>
            """

        html += """
            <h2>Methods</h2>
            <p>Analysis performed using the Calcium Channel Bioinformatics Analyzer.</p>
            <p>Domain detection based on conserved motif patterns from literature and Pfam database.</p>

            <h2>References</h2>
            <ul>
                <li>AlphaFold 3: Accurate structure prediction of biomolecular interactions</li>
                <li>Plant calcium-permeable channels (Swarbreck et al., 2013)</li>
                <li>Voltage-gated calcium channels (Catterall et al., 2020)</li>
            </ul>
        </body>
        </html>
        """

        return html

    def export_results(self):
        """General export function."""
        self.export_excel()

    def show_help(self):
        """Show user guide."""
        help_window = tk.Toplevel(self.root)
        help_window.title("User Guide")
        help_window.geometry("700x600")

        help_text = scrolledtext.ScrolledText(help_window, wrap=tk.WORD)
        help_text.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        guide = """
        CALCIUM CHANNEL BIOINFORMATICS ANALYZER - USER GUIDE
        =====================================================

        1. LOADING SEQUENCES
           - Load FASTA: Import sequences from FASTA files
           - Load GenBank: Fetch sequences from NCBI by accession
           - Add Manual Entry: Type/paste sequences directly

        2. DOMAIN ANALYSIS
           - Select domain types to search
           - Run detection to find conserved motifs
           - Visualize domain architecture
           - Supported domains:
             * EF-hand (Ca2+ binding)
             * CNBD (cyclic nucleotide binding)
             * Voltage sensors
             * Pore loops
             * SPRY domains (RyR-specific)
             * And more...

        3. PHYLOGENETIC ANALYSIS
           - Build trees using UPGMA or Neighbor-Joining
           - Export trees in Newick format
           - Requires ≥3 sequences

        4. STRUCTURE PREDICTION
           - Export sequences for AlphaFold3
           - Predict transmembrane domains
           - Visualize hydropathy plots

        5. CLASSIFICATION
           - Automatic channel family prediction
           - Based on:
             * Sequence descriptions
             * Domain patterns
             * Conserved motifs

        6. BLAST SEARCHES
           - Search NCBI databases
           - Note: Use DIAMOND/MMseqs2 for large-scale work

        7. EXPORT OPTIONS
           - HTML reports
           - Excel spreadsheets
           - FASTA files
           - Newick trees

        TIPS FOR OPTIMAL RESULTS:
        -------------------------
        • Start with well-annotated sequences from SwissProt
        • For plant channels, look for CNGC, GLR, TPC families
        • For animal channels, focus on Cav, RyR, IP3R families
        • Use AlphaFold3 for structure-function insights
        • Compare results across multiple species

        ADVANCED WORKFLOWS:
        ------------------
        1. Identify novel calcium channels:
           - BLAST search with known channels
           - Domain detection on hits
           - Classify by family
           - Predict structure with AlphaFold

        2. Evolutionary analysis:
           - Collect orthologs across species
           - Build phylogenetic tree
           - Map domain evolution

        3. Structure-function:
           - Identify key domains
           - Predict TM topology
           - Model with AlphaFold3
           - Analyze ligand binding sites

        For questions or issues, contact:
        George Dickinson
        UC Irvine, Department of Neurobiology & Behavior
        """

        help_text.insert(1.0, guide)
        help_text.config(state=tk.DISABLED)

    def show_about(self):
        """Show about dialog."""
        about_text = """
        Calcium Channel Bioinformatics Analyzer
        Version 2.0

        Developed by George Dickinson
        UC Irvine, Department of Neurobiology & Behavior

        A comprehensive platform for identifying and analyzing
        calcium channel proteins across species.

        Features:
        • Multi-species sequence analysis
        • Advanced domain detection
        • Phylogenetic reconstruction
        • AlphaFold3 integration
        • Interactive visualizations

        Built with:
        • Biopython
        • Matplotlib/Seaborn
        • Pandas/NumPy
        • Tkinter GUI

        For the latest bioinformatics tools and databases:
        • AlphaFold: https://alphafold.ebi.ac.uk/
        • NCBI: https://www.ncbi.nlm.nih.gov/
        • UniProt: https://www.uniprot.org/
        • Pfam: https://pfam.xfam.org/

        © 2024
        """

        messagebox.showinfo("About", about_text)

    def check_results(self):
        """Check for results from background threads."""
        try:
            while True:
                result_type, data = self.result_queue.get_nowait()

                self.hide_progress()

                if result_type == 'error':
                    messagebox.showerror("Error", data)
                    self.update_status("Error occurred")

                elif result_type == 'genbank':
                    success = [r for r in data if isinstance(r, SeqRecord)]
                    errors = [r for r in data if isinstance(r, str)]

                    self.sequences.extend(success)
                    self.update_sequence_list()

                    msg = f"Fetched {len(success)} sequence(s)"
                    if errors:
                        msg += f"\n\nErrors:\n" + "\n".join(errors)
                    messagebox.showinfo("GenBank Fetch", msg)
                    self.update_status(f"Loaded {len(success)} sequences from GenBank")

                elif result_type == 'domains':
                    self.domain_results = data

                    # Update treeview
                    for item in self.domain_tree.get_children():
                        self.domain_tree.delete(item)

                    for result in data:
                        self.domain_tree.insert('', tk.END, values=(
                            result['sequence'],
                            result['domain'],
                            result['position'],
                            result['match'],
                            result['score']
                        ))

                    self.update_status(f"Found {len(data)} domain matches")
                    messagebox.showinfo("Domain Detection",
                                      f"Found {len(data)} domain matches across {len(self.sequences)} sequences")

                elif result_type == 'alignment':
                    msg = f"Pairwise alignments completed\n\n"
                    for result in data[:5]:  # Show first 5
                        msg += f"{result['seq1']} vs {result['seq2']}: Score {result['score']:.1f}\n"
                    if len(data) > 5:
                        msg += f"\n... and {len(data)-5} more"

                    messagebox.showinfo("Alignment Results", msg)
                    self.update_status("Alignment complete")

                elif result_type == 'tree':
                    tree, alignment = data
                    self.tree = tree
                    self.alignment = alignment
                    self.plot_tree(tree)
                    self.update_status("Phylogenetic tree constructed")

                elif result_type == 'classification':
                    # Display classification results
                    results = "CALCIUM CHANNEL CLASSIFICATION RESULTS\n"
                    results += "="*60 + "\n\n"

                    for item in data:
                        results += f"Sequence: {item['id']}\n"
                        results += f"Description: {item['description']}\n"
                        results += f"Classification: {item['classification']}\n"
                        results += f"Confidence: {item['confidence']}\n"
                        results += f"Length: {item['length']} aa\n"
                        results += f"Detected markers:\n"
                        for marker, present in item['markers'].items():
                            results += f"  {marker}: {'Yes' if present else 'No'}\n"
                        results += "\n" + "-"*60 + "\n\n"

                    self.summary_text.delete(1.0, tk.END)
                    self.summary_text.insert(1.0, results)
                    self.notebook.select(4)  # Switch to results tab

                    self.update_status("Classification complete")
                    messagebox.showinfo("Classification", f"Classified {len(data)} sequences")

                elif result_type == 'blast':
                    blast_records = data
                    results = "BLAST SEARCH RESULTS\n"
                    results += "="*60 + "\n\n"

                    for record in blast_records:
                        for alignment in record.alignments[:10]:  # Top 10 hits
                            for hsp in alignment.hsps:
                                results += f"Hit: {alignment.title}\n"
                                results += f"E-value: {hsp.expect}\n"
                                results += f"Identity: {hsp.identities}/{hsp.align_length}\n"
                                results += f"Score: {hsp.score}\n\n"

                    self.summary_text.delete(1.0, tk.END)
                    self.summary_text.insert(1.0, results)
                    self.notebook.select(4)

                    self.update_status("BLAST search complete")
                    messagebox.showinfo("BLAST", "Search complete. See Results Summary tab.")

        except queue.Empty:
            pass

        # Schedule next check
        self.root.after(100, self.check_results)


def main():
    """Main entry point."""
    root = tk.Tk()
    app = CalciumChannelAnalyzer(root)
    root.mainloop()


if __name__ == "__main__":
    main()
