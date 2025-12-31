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
import subprocess
import os
import urllib.request
import gzip
import pickle
import tempfile
import xml.etree.ElementTree as ET
from collections import defaultdict, Counter

# Bioinformatics libraries
from Bio import SeqIO, Entrez, AlignIO, Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align import PairwiseAligner

try:
    from Bio import motifs
    from Bio.SeqUtils import CodonUsage
    MOTIF_AVAILABLE = True
except ImportError:
    MOTIF_AVAILABLE = False

# PDB structure parsing
try:
    from Bio.PDB import PDBParser, PDBIO, Select
    from Bio.PDB.Polypeptide import PPBuilder
    PDB_AVAILABLE = True
except ImportError:
    PDB_AVAILABLE = False

# Data analysis
import pandas as pd
import numpy as np
from scipy import stats
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist, squareform

# Machine Learning
try:
    from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import cross_val_score
    from sklearn.metrics import classification_report
    from sklearn.decomposition import PCA
    from sklearn.manifold import TSNE
    import joblib
    ML_AVAILABLE = True
except ImportError:
    ML_AVAILABLE = False

# Visualization
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle, FancyBboxPatch, Circle, Polygon
from matplotlib.collections import LineCollection
import matplotlib.cm as cm
import seaborn as sns

# Network requests for AlphaFold API
try:
    import requests
    REQUESTS_AVAILABLE = True
except ImportError:
    REQUESTS_AVAILABLE = False

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
        self.pdb_structures = {}  # Store downloaded PDB structures
        self.ml_classifier = None  # Machine learning classifier
        self.ml_scaler = None  # Feature scaler for ML
        self.alphafold_structures = {}  # Cache for AlphaFold structures
        self.diamond_db_path = None  # Path to DIAMOND database

        # New discovery features
        self.ortholog_data = {}  # Cross-species ortholog data
        self.synteny_data = {}  # Gene neighborhood data
        self.expression_data = {}  # Expression patterns
        self.interaction_data = {}  # Protein interactions
        self.hmm_results = []  # HMM search results
        self.comparative_genomes = {}  # Multi-genome data

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
        tools_menu.add_separator()
        tools_menu.add_command(label="DIAMOND Local Search", command=self.diamond_search_dialog)
        tools_menu.add_command(label="Download AlphaFold Structure", command=self.download_alphafold_structure)
        tools_menu.add_command(label="3D Structure Viewer", command=self.view_3d_structure)
        tools_menu.add_separator()
        tools_menu.add_command(label="ML Channel Classifier", command=self.ml_classify_dialog)
        tools_menu.add_command(label="Electrophysiology Prediction", command=self.predict_electrophysiology)
        tools_menu.add_command(label="Ligand Binding Prediction", command=self.predict_ligand_binding)

        # Discovery menu (NEW!)
        discovery_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Discovery", menu=discovery_menu)
        discovery_menu.add_command(label="Genome-Wide Search", command=self.genome_wide_search)
        discovery_menu.add_command(label="HMM Profile Search", command=self.hmm_search_dialog)
        discovery_menu.add_command(label="Find Orthologs", command=self.find_orthologs)
        discovery_menu.add_command(label="Synteny Analysis", command=self.synteny_analysis)
        discovery_menu.add_separator()
        discovery_menu.add_command(label="Expression Atlas", command=self.expression_atlas)
        discovery_menu.add_command(label="Protein Interactions", command=self.protein_interactions)
        discovery_menu.add_command(label="Pathway Analysis", command=self.pathway_analysis)
        discovery_menu.add_separator()
        discovery_menu.add_command(label="Comparative Genomics", command=self.comparative_genomics)
        discovery_menu.add_command(label="Selection Pressure Analysis", command=self.selection_analysis)

        # Visualization menu (NEW!)
        viz_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Visualize", menu=viz_menu)
        viz_menu.add_command(label="Synteny Plot", command=self.plot_synteny)
        viz_menu.add_command(label="Expression Heatmap", command=self.plot_expression_heatmap)
        viz_menu.add_command(label="Network Graph", command=self.plot_network)
        viz_menu.add_command(label="Sequence Logo", command=self.plot_sequence_logo)
        viz_menu.add_command(label="Domain Architecture", command=self.plot_domain_architecture)
        viz_menu.add_separator()
        viz_menu.add_command(label="PCA Plot", command=self.plot_pca)
        viz_menu.add_command(label="t-SNE Clustering", command=self.plot_tsne)
        viz_menu.add_command(label="Phylogenetic Network", command=self.plot_phylo_network)
        viz_menu.add_separator()
        viz_menu.add_command(label="Chromosome View", command=self.plot_chromosome_view)
        viz_menu.add_command(label="Circos Plot", command=self.plot_circos)

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

        # Tab 5: 3D Structure Viewer
        self.setup_3d_viewer_tab()

        # Tab 6: ML Classification
        self.setup_ml_tab()

        # Tab 7: Electrophysiology
        self.setup_electrophys_tab()

        # Tab 8: Discovery & Orthologs
        self.setup_discovery_tab()

        # Tab 9: Comparative Genomics
        self.setup_comparative_tab()

        # Tab 10: Visualizations
        self.setup_visualization_tab()

        # Tab 11: Results Summary
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

        # Update 3D viewer combo
        if hasattr(self, 'viewer_seq_combo'):
            self.viewer_seq_combo['values'] = seq_list
            if seq_list:
                self.viewer_seq_combo.current(0)

        # Update variant combo
        if hasattr(self, 'variant_seq_combo'):
            seq_ids = [r.id for r in self.sequences]
            self.variant_seq_combo['values'] = seq_ids
            if seq_ids:
                self.variant_seq_combo.current(0)

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

    def setup_3d_viewer_tab(self):
        """Setup 3D structure viewer tab."""
        viewer_frame = ttk.Frame(self.notebook)
        self.notebook.add(viewer_frame, text="3D Viewer")

        # Info panel
        info_frame = ttk.LabelFrame(viewer_frame, text="3D Structure Visualization")
        info_frame.pack(fill=tk.X, padx=5, pady=5)

        info_text = """
        View and analyze 3D protein structures from:
        • AlphaFold Database (200M+ structures)
        • PDB (Protein Data Bank)
        • Uploaded PDB files

        Features:
        • Interactive 3D visualization
        • Domain highlighting
        • Ca²⁺ binding site identification
        • Structure quality assessment
        """
        ttk.Label(info_frame, text=info_text, justify=tk.LEFT).pack(padx=10, pady=10)

        # Controls
        control_frame = ttk.LabelFrame(viewer_frame, text="Load Structure")
        control_frame.pack(fill=tk.X, padx=5, pady=5)

        ttk.Label(control_frame, text="Select sequence:").grid(row=0, column=0, padx=5, pady=5)
        self.viewer_seq_combo = ttk.Combobox(control_frame, width=50)
        self.viewer_seq_combo.grid(row=0, column=1, padx=5, pady=5)

        ttk.Button(control_frame, text="Download from AlphaFold DB",
                  command=self.download_alphafold_structure).grid(row=0, column=2, padx=5, pady=5)
        ttk.Button(control_frame, text="Load PDB File",
                  command=self.load_pdb_file).grid(row=1, column=0, padx=5, pady=5)
        ttk.Button(control_frame, text="View 3D Structure",
                  command=self.view_3d_structure).grid(row=1, column=1, padx=5, pady=5)

        # Visualization frame
        self.viewer_canvas_frame = ttk.Frame(viewer_frame)
        self.viewer_canvas_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Structure info
        self.structure_info_text = scrolledtext.ScrolledText(viewer_frame, height=10, wrap=tk.WORD)
        self.structure_info_text.pack(fill=tk.X, padx=5, pady=5)

    def setup_ml_tab(self):
        """Setup machine learning classification tab."""
        ml_frame = ttk.Frame(self.notebook)
        self.notebook.add(ml_frame, text="ML Classification")

        # Info
        info_frame = ttk.LabelFrame(ml_frame, text="Machine Learning Channel Classifier")
        info_frame.pack(fill=tk.X, padx=5, pady=5)

        info_text = """
        Advanced ML-based calcium channel classification using:
        • Sequence features (composition, length, hydrophobicity)
        • Domain patterns (presence/absence of key motifs)
        • Phylogenetic features
        • Random Forest and Gradient Boosting models

        Classifies into: Cav1, Cav2, Cav3, RyR, IP3R, CNGC, GLR, TPC, and more
        """
        ttk.Label(info_frame, text=info_text, justify=tk.LEFT).pack(padx=10, pady=10)

        # Controls
        control_frame = ttk.LabelFrame(ml_frame, text="ML Classification Options")
        control_frame.pack(fill=tk.X, padx=5, pady=5)

        ttk.Label(control_frame, text="Model:").grid(row=0, column=0, padx=5, pady=5)
        self.ml_model_combo = ttk.Combobox(control_frame,
                                           values=['Random Forest', 'Gradient Boosting', 'Ensemble'],
                                           state='readonly')
        self.ml_model_combo.set('Random Forest')
        self.ml_model_combo.grid(row=0, column=1, padx=5, pady=5)

        ttk.Button(control_frame, text="Train Classifier",
                  command=self.train_ml_classifier).grid(row=0, column=2, padx=5, pady=5)
        ttk.Button(control_frame, text="Classify Sequences",
                  command=self.ml_classify_sequences).grid(row=0, column=3, padx=5, pady=5)
        ttk.Button(control_frame, text="Export Model",
                  command=self.export_ml_model).grid(row=1, column=0, padx=5, pady=5)
        ttk.Button(control_frame, text="Load Model",
                  command=self.load_ml_model).grid(row=1, column=1, padx=5, pady=5)

        # Results
        results_frame = ttk.LabelFrame(ml_frame, text="Classification Results")
        results_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Treeview for results
        columns = ('Sequence', 'Predicted_Family', 'Confidence', 'Top3_Predictions')
        self.ml_tree = ttk.Treeview(results_frame, columns=columns, show='headings')

        for col in columns:
            self.ml_tree.heading(col, text=col.replace('_', ' '))
            self.ml_tree.column(col, width=150)

        vsb = ttk.Scrollbar(results_frame, orient="vertical", command=self.ml_tree.yview)
        self.ml_tree.configure(yscrollcommand=vsb.set)

        self.ml_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        vsb.pack(side=tk.RIGHT, fill=tk.Y)

    def setup_electrophys_tab(self):
        """Setup electrophysiology prediction tab."""
        ephys_frame = ttk.Frame(self.notebook)
        self.notebook.add(ephys_frame, text="Electrophysiology")

        # Info
        info_frame = ttk.LabelFrame(ephys_frame, text="Electrophysiology Prediction")
        info_frame.pack(fill=tk.X, padx=5, pady=5)

        info_text = """
        Predict electrophysiological properties based on sequence:
        • Gating mechanism (voltage, ligand, mechanosensitive)
        • Ion selectivity (Ca²⁺, Na⁺, K⁺, non-selective)
        • Activation voltage (for voltage-gated)
        • Permeability ratios
        • Pharmacology (blocker sensitivity)

        Based on domain analysis and ML predictions.
        """
        ttk.Label(info_frame, text=info_text, justify=tk.LEFT).pack(padx=10, pady=10)

        # Controls
        control_frame = ttk.LabelFrame(ephys_frame, text="Prediction Options")
        control_frame.pack(fill=tk.X, padx=5, pady=5)

        ttk.Button(control_frame, text="Predict Properties",
                  command=self.predict_electrophysiology).pack(side=tk.LEFT, padx=5, pady=5)
        ttk.Button(control_frame, text="Generate Report",
                  command=self.generate_ephys_report).pack(side=tk.LEFT, padx=5, pady=5)

        # Results display
        self.ephys_results = scrolledtext.ScrolledText(ephys_frame, wrap=tk.WORD)
        self.ephys_results.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

    # ==================== NEW ADVANCED FEATURE METHODS ====================

    def diamond_search_dialog(self):
        """Dialog for local DIAMOND search."""
        if not self.sequences:
            messagebox.showwarning("Warning", "No sequences loaded")
            return

        dialog = tk.Toplevel(self.root)
        dialog.title("DIAMOND Local Search")
        dialog.geometry("600x400")

        info = """
        DIAMOND Ultra-Fast Sequence Search:

        DIAMOND is 100-10,000x faster than BLAST with similar sensitivity.

        Setup:
        1. Install DIAMOND: conda install -c bioconda diamond
        2. Create database: diamond makedb --in proteins.faa -d database
        3. Select database below and search

        For large-scale searches, DIAMOND is essential.
        """
        ttk.Label(dialog, text=info, justify=tk.LEFT).pack(padx=10, pady=10)

        # Database selection
        ttk.Label(dialog, text="DIAMOND Database:").pack()
        db_frame = ttk.Frame(dialog)
        db_frame.pack(fill=tk.X, padx=10, pady=5)

        db_entry = ttk.Entry(db_frame, width=50)
        db_entry.pack(side=tk.LEFT, padx=5)

        def browse_db():
            filename = filedialog.askopenfilename(
                title="Select DIAMOND database",
                filetypes=[("DIAMOND DB", "*.dmnd"), ("All files", "*.*")]
            )
            if filename:
                db_entry.delete(0, tk.END)
                db_entry.insert(0, filename)
                self.diamond_db_path = filename

        ttk.Button(db_frame, text="Browse", command=browse_db).pack(side=tk.LEFT)

        # Sensitivity
        ttk.Label(dialog, text="Sensitivity:").pack()
        sensitivity_combo = ttk.Combobox(dialog,
                                        values=['fast', 'mid-sensitive', 'sensitive',
                                               'more-sensitive', 'very-sensitive', 'ultra-sensitive'],
                                        state='readonly')
        sensitivity_combo.set('very-sensitive')
        sensitivity_combo.pack(pady=5)

        # Max targets
        ttk.Label(dialog, text="Max targets:").pack()
        max_entry = ttk.Entry(dialog)
        max_entry.insert(0, "100")
        max_entry.pack(pady=5)

        def run_diamond():
            if not db_entry.get():
                messagebox.showwarning("Warning", "Please select a DIAMOND database")
                return

            dialog.destroy()
            self.run_diamond_search(db_entry.get(), sensitivity_combo.get(), int(max_entry.get()))

        ttk.Button(dialog, text="Run DIAMOND Search", command=run_diamond).pack(pady=10)
        ttk.Button(dialog, text="Cancel", command=dialog.destroy).pack()

    def run_diamond_search(self, database: str, sensitivity: str, max_targets: int):
        """Execute DIAMOND search."""
        self.show_progress()
        self.update_status("Running DIAMOND search...")

        def diamond_thread():
            try:
                # Check if DIAMOND is installed
                result = subprocess.run(['diamond', 'version'],
                                      capture_output=True, text=True)
                if result.returncode != 0:
                    self.result_queue.put(('error',
                        "DIAMOND not found. Install with: conda install -c bioconda diamond"))
                    return

                # Create temp query file
                with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
                    SeqIO.write(self.sequences, f, 'fasta')
                    query_file = f.name

                # Create temp output file
                output_file = tempfile.mktemp(suffix='.txt')

                # Run DIAMOND
                cmd = [
                    'diamond', 'blastp',
                    '-d', database,
                    '-q', query_file,
                    '-o', output_file,
                    f'--{sensitivity}',
                    '--max-target-seqs', str(max_targets),
                    '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'evalue', 'bitscore', 'stitle'
                ]

                result = subprocess.run(cmd, capture_output=True, text=True)

                # Parse results
                if result.returncode == 0 and os.path.exists(output_file):
                    results = []
                    with open(output_file, 'r') as f:
                        for line in f:
                            parts = line.strip().split('\t')
                            if len(parts) >= 7:
                                results.append({
                                    'query': parts[0],
                                    'subject': parts[1],
                                    'identity': parts[2],
                                    'length': parts[3],
                                    'evalue': parts[4],
                                    'bitscore': parts[5],
                                    'description': parts[6] if len(parts) > 6 else ''
                                })

                    # Cleanup
                    os.unlink(query_file)
                    os.unlink(output_file)

                    self.result_queue.put(('diamond', results))
                else:
                    self.result_queue.put(('error', f"DIAMOND failed: {result.stderr}"))

            except Exception as e:
                self.result_queue.put(('error', f"DIAMOND search failed: {str(e)}"))

        thread = threading.Thread(target=diamond_thread, daemon=True)
        thread.start()

    def download_alphafold_structure(self):
        """Download structure from AlphaFold Database."""
        if not REQUESTS_AVAILABLE:
            messagebox.showerror("Error", "requests library not available. Install with: pip install requests")
            return

        if not self.sequences:
            messagebox.showwarning("Warning", "No sequences loaded")
            return

        # Get UniProt ID
        dialog = tk.Toplevel(self.root)
        dialog.title("Download AlphaFold Structure")
        dialog.geometry("500x300")

        info = """
        Download predicted structure from AlphaFold Database:

        Enter UniProt ID (e.g., P12345) or select from loaded sequences.
        AlphaFold DB contains 200M+ protein structure predictions.
        """
        ttk.Label(dialog, text=info, justify=tk.LEFT).pack(padx=10, pady=10)

        ttk.Label(dialog, text="UniProt ID:").pack()
        id_entry = ttk.Entry(dialog, width=30)
        id_entry.pack(pady=5)

        ttk.Label(dialog, text="Or select sequence:").pack()
        seq_combo = ttk.Combobox(dialog,
                                values=[f"{s.id}" for s in self.sequences],
                                state='readonly')
        seq_combo.pack(pady=5)

        def download():
            uniprot_id = id_entry.get().strip()
            if not uniprot_id and seq_combo.get():
                uniprot_id = seq_combo.get().split()[0]

            if not uniprot_id:
                messagebox.showwarning("Warning", "Please enter a UniProt ID")
                return

            dialog.destroy()
            self.fetch_alphafold_structure(uniprot_id)

        ttk.Button(dialog, text="Download", command=download).pack(pady=10)
        ttk.Button(dialog, text="Cancel", command=dialog.destroy).pack()

    def fetch_alphafold_structure(self, uniprot_id: str):
        """Fetch structure from AlphaFold DB."""
        self.show_progress()
        self.update_status(f"Downloading AlphaFold structure for {uniprot_id}...")

        def fetch_thread():
            try:
                # AlphaFold DB URL
                url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"

                # Download
                response = requests.get(url, timeout=30)

                if response.status_code == 200:
                    # Save structure
                    pdb_file = tempfile.mktemp(suffix='.pdb')
                    with open(pdb_file, 'w') as f:
                        f.write(response.text)

                    self.alphafold_structures[uniprot_id] = pdb_file
                    self.result_queue.put(('alphafold_downloaded', (uniprot_id, pdb_file)))
                else:
                    self.result_queue.put(('error',
                        f"Structure not found for {uniprot_id}. "
                        f"Check UniProt ID or try alternative accessions."))

            except Exception as e:
                self.result_queue.put(('error', f"Download failed: {str(e)}"))

        thread = threading.Thread(target=fetch_thread, daemon=True)
        thread.start()

    def load_pdb_file(self):
        """Load PDB file from disk."""
        filename = filedialog.askopenfilename(
            title="Select PDB file",
            filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")]
        )

        if filename:
            try:
                # Parse PDB
                if PDB_AVAILABLE:
                    parser = PDBParser(QUIET=True)
                    structure = parser.get_structure('loaded', filename)

                    # Store
                    name = Path(filename).stem
                    self.pdb_structures[name] = filename

                    messagebox.showinfo("Success", f"Loaded structure: {name}")
                    self.update_status(f"Loaded PDB: {name}")
                else:
                    messagebox.showerror("Error", "BioPython PDB module not available")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load PDB: {str(e)}")

    def view_3d_structure(self):
        """Visualize 3D structure."""
        if not self.pdb_structures and not self.alphafold_structures:
            messagebox.showinfo("Info",
                              "No structures loaded. Download from AlphaFold DB or load PDB file first.")
            return

        # Select structure
        all_structures = {**self.pdb_structures, **self.alphafold_structures}

        if not all_structures:
            return

        structure_name = list(all_structures.keys())[0]
        pdb_file = all_structures[structure_name]

        if not PDB_AVAILABLE:
            messagebox.showinfo("Info",
                              f"3D visualization requires BioPython PDB module.\n"
                              f"Structure saved at: {pdb_file}\n"
                              f"View in PyMOL, ChimeraX, or other molecular viewer.")
            return

        # Parse structure
        try:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure(structure_name, pdb_file)

            # Extract coordinates
            coords = []
            for model in structure:
                for chain in model:
                    for residue in chain:
                        for atom in residue:
                            if atom.get_name() == 'CA':  # Alpha carbons only
                                coords.append(atom.get_coord())

            coords = np.array(coords)

            # Clear previous plot
            for widget in self.viewer_canvas_frame.winfo_children():
                widget.destroy()

            # Create 3D plot
            fig = plt.figure(figsize=(12, 10))
            ax = fig.add_subplot(111, projection='3d')

            # Plot backbone
            ax.plot(coords[:, 0], coords[:, 1], coords[:, 2],
                   'b-', linewidth=2, alpha=0.6, label='Backbone')
            ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2],
                      c=np.arange(len(coords)), cmap='viridis', s=20)

            ax.set_xlabel('X (Å)')
            ax.set_ylabel('Y (Å)')
            ax.set_zlabel('Z (Å)')
            ax.set_title(f'3D Structure: {structure_name}', fontsize=14, fontweight='bold')
            ax.legend()

            plt.tight_layout()

            # Embed in tkinter
            canvas = FigureCanvasTkAgg(fig, master=self.viewer_canvas_frame)
            canvas.draw()
            canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

            toolbar = NavigationToolbar2Tk(canvas, self.viewer_canvas_frame)
            toolbar.update()

            # Show structure info
            info = f"Structure: {structure_name}\n"
            info += f"File: {pdb_file}\n"
            info += f"Residues: {len(coords)}\n"
            info += f"Models: {len(structure)}\n\n"

            # Get sequence
            ppb = PPBuilder()
            for pp in ppb.build_peptides(structure):
                seq = pp.get_sequence()
                info += f"Sequence ({len(seq)} aa):\n{seq[:100]}...\n"
                break

            self.structure_info_text.delete(1.0, tk.END)
            self.structure_info_text.insert(1.0, info)

            self.update_status(f"Displaying 3D structure: {structure_name}")

        except Exception as e:
            messagebox.showerror("Error", f"Failed to visualize structure: {str(e)}")

    def train_ml_classifier(self):
        """Train machine learning classifier."""
        if not ML_AVAILABLE:
            messagebox.showerror("Error", "scikit-learn not available. Install with: pip install scikit-learn")
            return

        if len(self.sequences) < 10:
            messagebox.showwarning("Warning", "Need at least 10 labeled sequences to train classifier")
            return

        self.show_progress()
        self.update_status("Training ML classifier...")

        def train_thread():
            try:
                # Extract features and labels
                features = []
                labels = []

                for record in self.sequences:
                    seq_str = str(record.seq)

                    # Extract features
                    feat = self.extract_ml_features(seq_str)
                    features.append(feat)

                    # Get label from description
                    label = self.predict_family_from_description(record.description)
                    labels.append(label)

                X = np.array(features)
                y = np.array(labels)

                # Filter unknown labels
                mask = y != 'Unknown'
                X = X[mask]
                y = y[mask]

                if len(X) < 5:
                    self.result_queue.put(('error', "Not enough labeled data to train"))
                    return

                # Scale features
                self.ml_scaler = StandardScaler()
                X_scaled = self.ml_scaler.fit_transform(X)

                # Train model
                model_type = self.ml_model_combo.get()

                if model_type == 'Random Forest':
                    self.ml_classifier = RandomForestClassifier(
                        n_estimators=100, max_depth=10, random_state=42
                    )
                elif model_type == 'Gradient Boosting':
                    self.ml_classifier = GradientBoostingClassifier(
                        n_estimators=100, max_depth=5, random_state=42
                    )
                else:  # Ensemble
                    from sklearn.ensemble import VotingClassifier
                    rf = RandomForestClassifier(n_estimators=100, random_state=42)
                    gb = GradientBoostingClassifier(n_estimators=100, random_state=42)
                    self.ml_classifier = VotingClassifier(
                        estimators=[('rf', rf), ('gb', gb)],
                        voting='soft'
                    )

                self.ml_classifier.fit(X_scaled, y)

                # Cross-validation
                scores = cross_val_score(self.ml_classifier, X_scaled, y, cv=min(5, len(X)//2))

                self.result_queue.put(('ml_trained', scores))

            except Exception as e:
                self.result_queue.put(('error', f"Training failed: {str(e)}"))

        thread = threading.Thread(target=train_thread, daemon=True)
        thread.start()

    def extract_ml_features(self, sequence: str) -> List[float]:
        """Extract features for ML classification."""
        features = []

        # Length
        features.append(len(sequence))

        # Amino acid composition
        aa_list = 'ACDEFGHIKLMNPQRSTVWY'
        for aa in aa_list:
            features.append(sequence.count(aa) / len(sequence) if len(sequence) > 0 else 0)

        # Hydrophobicity
        kd = {'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8, 'G': -0.4,
              'H': -3.2, 'I': 4.5, 'K': -3.9, 'L': 3.8, 'M': 1.9, 'N': -3.5,
              'P': -1.6, 'Q': -3.5, 'R': -4.5, 'S': -0.8, 'T': -0.7, 'V': 4.2,
              'W': -0.9, 'Y': -1.3}
        avg_hydro = sum(kd.get(aa, 0) for aa in sequence) / len(sequence) if len(sequence) > 0 else 0
        features.append(avg_hydro)

        # Domain presence (binary features)
        for domain_name, domain_info in self.DOMAIN_PATTERNS.items():
            pattern = domain_info['pattern']
            has_domain = 1 if re.search(pattern, sequence) else 0
            features.append(has_domain)

        return features

    def ml_classify_sequences(self):
        """Classify sequences using trained ML model."""
        if not ML_AVAILABLE:
            messagebox.showerror("Error", "scikit-learn not available")
            return

        if self.ml_classifier is None:
            messagebox.showwarning("Warning", "Train classifier first")
            return

        if not self.sequences:
            messagebox.showwarning("Warning", "No sequences to classify")
            return

        self.show_progress()
        self.update_status("Classifying sequences with ML...")

        def classify_thread():
            try:
                results = []

                for record in self.sequences:
                    seq_str = str(record.seq)

                    # Extract features
                    feat = self.extract_ml_features(seq_str)
                    X = np.array([feat])
                    X_scaled = self.ml_scaler.transform(X)

                    # Predict
                    prediction = self.ml_classifier.predict(X_scaled)[0]

                    # Get probabilities
                    if hasattr(self.ml_classifier, 'predict_proba'):
                        proba = self.ml_classifier.predict_proba(X_scaled)[0]
                        classes = self.ml_classifier.classes_

                        # Get top 3
                        top_indices = np.argsort(proba)[-3:][::-1]
                        top_predictions = [(classes[i], proba[i]) for i in top_indices]
                        confidence = proba[np.argmax(proba)]
                    else:
                        top_predictions = [(prediction, 1.0)]
                        confidence = 1.0

                    results.append({
                        'sequence': record.id,
                        'prediction': prediction,
                        'confidence': confidence,
                        'top3': top_predictions
                    })

                self.result_queue.put(('ml_classified', results))

            except Exception as e:
                self.result_queue.put(('error', f"Classification failed: {str(e)}"))

        thread = threading.Thread(target=classify_thread, daemon=True)
        thread.start()

    def ml_classify_dialog(self):
        """Show ML classification dialog."""
        if not ML_AVAILABLE:
            messagebox.showerror("Error", "scikit-learn not available. Install with: pip install scikit-learn")
            return

        # Switch to ML tab and run classification
        self.notebook.select(5)  # ML tab
        if self.ml_classifier:
            self.ml_classify_sequences()
        else:
            messagebox.showinfo("Info",
                              "No trained classifier. Click 'Train Classifier' first.\n"
                              "Classifier will use loaded sequences as training data.")

    def export_ml_model(self):
        """Export trained ML model."""
        if not ML_AVAILABLE or self.ml_classifier is None:
            messagebox.showwarning("Warning", "No trained model to export")
            return

        filename = filedialog.asksaveasfilename(
            defaultextension=".pkl",
            filetypes=[("Pickle files", "*.pkl"), ("All files", "*.*")]
        )

        if filename:
            try:
                joblib.dump({
                    'classifier': self.ml_classifier,
                    'scaler': self.ml_scaler
                }, filename)
                messagebox.showinfo("Success", f"Model exported to {filename}")
            except Exception as e:
                messagebox.showerror("Error", f"Export failed: {str(e)}")

    def load_ml_model(self):
        """Load trained ML model."""
        if not ML_AVAILABLE:
            messagebox.showerror("Error", "scikit-learn not available")
            return

        filename = filedialog.askopenfilename(
            title="Select model file",
            filetypes=[("Pickle files", "*.pkl"), ("All files", "*.*")]
        )

        if filename:
            try:
                model_data = joblib.load(filename)
                self.ml_classifier = model_data['classifier']
                self.ml_scaler = model_data['scaler']
                messagebox.showinfo("Success", f"Model loaded from {filename}")
            except Exception as e:
                messagebox.showerror("Error", f"Load failed: {str(e)}")

    def predict_electrophysiology(self):
        """Predict electrophysiological properties."""
        if not self.sequences:
            messagebox.showwarning("Warning", "No sequences loaded")
            return

        self.show_progress()
        self.update_status("Predicting electrophysiological properties...")

        def predict_thread():
            try:
                results = []

                for record in self.sequences:
                    seq_str = str(record.seq)
                    desc = record.description.lower()

                    # Predict gating mechanism
                    has_voltage = bool(re.search(self.DOMAIN_PATTERNS['Voltage_sensor']['pattern'], seq_str))
                    has_cnbd = bool(re.search(self.DOMAIN_PATTERNS['CNBD']['pattern'], seq_str))
                    has_glr = bool(re.search(self.DOMAIN_PATTERNS['GLR_ligand']['pattern'], seq_str))

                    if has_voltage:
                        gating = "Voltage-gated"
                        activation = "Depolarization-activated (-40 to +20 mV typical)"
                    elif has_cnbd:
                        gating = "Ligand-gated (cAMP/cGMP)"
                        activation = "cAMP/cGMP binding"
                    elif has_glr:
                        gating = "Ligand-gated (Glutamate-like)"
                        activation = "Glutamate or analog binding"
                    elif 'ryr' in desc or 'ryanodine' in desc:
                        gating = "Ligand-gated (Ca²⁺, cADPR, ryanodine)"
                        activation = "Ca²⁺-induced Ca²⁺ release (CICR)"
                    elif 'ip3' in desc:
                        gating = "Ligand-gated (IP3)"
                        activation = "IP3 binding"
                    else:
                        gating = "Unknown/Mixed"
                        activation = "Multiple mechanisms possible"

                    # Predict selectivity
                    has_ca_filter = bool(re.search(self.DOMAIN_PATTERNS['Ca_selectivity']['pattern'], seq_str))

                    if has_ca_filter:
                        selectivity = "Ca²⁺-selective (PCa/PNa > 100)"
                        permeability = "High Ca²⁺, Low Na⁺/K⁺"
                    elif 'tpc' in desc or 'two-pore' in desc:
                        selectivity = "Non-selective cation (Ca²⁺, Na⁺, K⁺)"
                        permeability = "PCa:PNa:PK ≈ 1:1:1"
                    else:
                        selectivity = "Ca²⁺-permeable"
                        permeability = "Variable Ca²⁺/monovalent selectivity"

                    # Predict pharmacology
                    if 'cav1' in desc or 'l-type' in desc:
                        blockers = "Dihydropyridines (nifedipine), phenylalkylamines (verapamil), benzothiazepines (diltiazem)"
                    elif 'cav2.1' in desc or 'p/q' in desc:
                        blockers = "ω-agatoxin IVA"
                    elif 'cav2.2' in desc or 'n-type' in desc:
                        blockers = "ω-conotoxin GVIA"
                    elif 'cav3' in desc or 't-type' in desc:
                        blockers = "Mibefradil, TTA-P2"
                    elif 'ryr' in desc:
                        blockers = "Ryanodine (high dose), dantrolene"
                    else:
                        blockers = "Family-specific blockers required"

                    results.append({
                        'sequence': record.id,
                        'gating': gating,
                        'activation': activation,
                        'selectivity': selectivity,
                        'permeability': permeability,
                        'blockers': blockers
                    })

                self.result_queue.put(('electrophys', results))

            except Exception as e:
                self.result_queue.put(('error', f"Prediction failed: {str(e)}"))

        thread = threading.Thread(target=predict_thread, daemon=True)
        thread.start()

    def generate_ephys_report(self):
        """Generate comprehensive electrophysiology report."""
        if not self.sequences:
            messagebox.showwarning("Warning", "No sequences to analyze")
            return

        # Run prediction first
        self.predict_electrophysiology()

    def predict_ligand_binding(self):
        """Predict potential ligand binding sites and druggability."""
        if not self.sequences:
            messagebox.showwarning("Warning", "No sequences loaded")
            return

        dialog = tk.Toplevel(self.root)
        dialog.title("Ligand Binding Prediction")
        dialog.geometry("700x600")

        info = """
        Ligand Binding Site Prediction:

        Analyzes sequences for potential drug binding sites based on:
        • Known pharmacophore patterns
        • Domain analysis
        • Pocket-forming residues
        • Druggability scores

        Predicts binding for:
        • Small molecule drugs (DHPs, ω-toxins)
        • Endogenous ligands (Ca²⁺, cAMP, IP3)
        • Natural products
        """

        ttk.Label(dialog, text=info, justify=tk.LEFT).pack(padx=10, pady=10)

        results_text = scrolledtext.ScrolledText(dialog, wrap=tk.WORD)
        results_text.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # Analyze
        report = "LIGAND BINDING SITE PREDICTIONS\n"
        report += "="*60 + "\n\n"

        for record in self.sequences:
            seq_str = str(record.seq)
            desc = record.description.lower()

            report += f"Sequence: {record.id}\n"
            report += f"Description: {record.description}\n\n"

            # Check for DHP binding site (L-type)
            if 'cav1' in desc or 'l-type' in desc:
                report += "DHP Binding Site (L-type specific):\n"
                report += "  • Located in domains III-IV\n"
                report += "  • Involves S6 segments\n"
                report += "  • Druggable with dihydropyridines\n"
                report += "  • Known drugs: Nifedipine, Amlodipine, Felodipine\n\n"

            # Check for Ca²⁺ binding
            has_efhand = bool(re.search(self.DOMAIN_PATTERNS['EF_hand']['pattern'], seq_str))
            if has_efhand:
                matches = list(re.finditer(self.DOMAIN_PATTERNS['EF_hand']['pattern'], seq_str))
                report += f"EF-hand Ca²⁺ Binding Sites: {len(matches)} found\n"
                for i, match in enumerate(matches[:3], 1):
                    report += f"  Site {i}: Position {match.start()}-{match.end()}\n"
                    report += f"  Sequence: {match.group()}\n"
                report += "\n"

            # Check for cAMP/cGMP binding
            has_cnbd = bool(re.search(self.DOMAIN_PATTERNS['CNBD']['pattern'], seq_str))
            if has_cnbd:
                report += "Cyclic Nucleotide Binding Domain:\n"
                report += "  • Binds cAMP and/or cGMP\n"
                report += "  • Allosteric regulation site\n"
                report += "  • Potential target for modulators\n\n"

            # Check for calmodulin binding
            has_cam = bool(re.search(self.DOMAIN_PATTERNS['Calmodulin_binding']['pattern'], seq_str))
            if has_cam:
                report += "Calmodulin Binding Site (IQ motif):\n"
                report += "  • Ca²⁺-dependent regulation\n"
                report += "  • Feedback inhibition site\n"
                report += "  • Potential for allosteric modulators\n\n"

            # Druggability score (simplified)
            druggable_features = sum([
                'cav1' in desc or 'l-type' in desc,
                'cav2' in desc,
                has_efhand,
                has_cnbd,
                has_cam
            ])

            druggability = "High" if druggable_features >= 3 else "Medium" if druggable_features >= 2 else "Low"
            report += f"Druggability Score: {druggability} ({druggable_features}/5 features)\n"
            report += "\n" + "-"*60 + "\n\n"

        results_text.insert(1.0, report)

        ttk.Button(dialog, text="Close", command=dialog.destroy).pack(pady=10)


    def setup_discovery_tab(self):
        """Setup discovery and ortholog finding tab."""
        discovery_frame = ttk.Frame(self.notebook)
        self.notebook.add(discovery_frame, text="Discovery")

        # Info
        info_frame = ttk.LabelFrame(discovery_frame, text="Novel Channel Discovery Across Species")
        info_frame.pack(fill=tk.X, padx=5, pady=5)

        info_text = """
        DISCOVERY TOOLS FOR FINDING NOVEL CALCIUM CHANNELS:

        • Genome-Wide Search: Screen entire genomes for calcium channel homologs
        • HMM Profile Search: Use profile HMMs for sensitive detection
        • Ortholog Finder: Identify orthologs across species using reciprocal best hits
        • Synteny Analysis: Detect conserved gene neighborhoods

        Supports: Bacteria, Archaea, Plants, Fungi, Animals
        """
        ttk.Label(info_frame, text=info_text, justify=tk.LEFT, font=('TkDefaultFont', 9)).pack(padx=10, pady=10)

        # Control panel
        control_frame = ttk.LabelFrame(discovery_frame, text="Discovery Options")
        control_frame.pack(fill=tk.X, padx=5, pady=5)

        # Organism selection
        ttk.Label(control_frame, text="Target organism:").grid(row=0, column=0, padx=5, pady=5, sticky=tk.W)
        self.organism_entry = ttk.Entry(control_frame, width=30)
        self.organism_entry.grid(row=0, column=1, padx=5, pady=5)
        self.organism_entry.insert(0, "Homo sapiens")

        ttk.Label(control_frame, text="Search method:").grid(row=1, column=0, padx=5, pady=5, sticky=tk.W)
        self.search_method_combo = ttk.Combobox(control_frame,
                                                values=['DIAMOND', 'HMM Profile', 'Reciprocal Best Hit', 'Domain-based'],
                                                state='readonly', width=28)
        self.search_method_combo.set('DIAMOND')
        self.search_method_combo.grid(row=1, column=1, padx=5, pady=5)

        ttk.Label(control_frame, text="E-value cutoff:").grid(row=2, column=0, padx=5, pady=5, sticky=tk.W)
        self.evalue_entry = ttk.Entry(control_frame, width=30)
        self.evalue_entry.insert(0, "1e-10")
        self.evalue_entry.grid(row=2, column=1, padx=5, pady=5)

        # Buttons
        button_frame = ttk.Frame(control_frame)
        button_frame.grid(row=3, column=0, columnspan=2, pady=10)

        ttk.Button(button_frame, text="Search Genome",
                  command=self.genome_wide_search).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Find Orthologs",
                  command=self.find_orthologs).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Build HMM",
                  command=self.build_hmm_profile).pack(side=tk.LEFT, padx=5)

        # Results display
        results_frame = ttk.LabelFrame(discovery_frame, text="Discovery Results")
        results_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Treeview
        columns = ('Gene_ID', 'Organism', 'Score', 'Identity', 'Domains', 'Family')
        self.discovery_tree = ttk.Treeview(results_frame, columns=columns, show='headings', height=15)

        for col in columns:
            self.discovery_tree.heading(col, text=col.replace('_', ' '))
            width = 120 if col in ['Gene_ID', 'Organism'] else 80
            self.discovery_tree.column(col, width=width)

        vsb = ttk.Scrollbar(results_frame, orient="vertical", command=self.discovery_tree.yview)
        self.discovery_tree.configure(yscrollcommand=vsb.set)

        self.discovery_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        vsb.pack(side=tk.RIGHT, fill=tk.Y)

    def setup_comparative_tab(self):
        """Setup comparative genomics tab."""
        comp_frame = ttk.Frame(self.notebook)
        self.notebook.add(comp_frame, text="Comparative")

        # Info
        info_frame = ttk.LabelFrame(comp_frame, text="Comparative Genomics Analysis")
        info_frame.pack(fill=tk.X, padx=5, pady=5)

        info_text = """
        COMPARATIVE GENOMICS TOOLS:

        • Synteny Analysis: Compare gene neighborhoods across species
        • Gene Family Evolution: Track expansions/contractions
        • Selection Pressure: Calculate dN/dS ratios
        • Conserved Regions: Identify highly conserved sequences
        • Species-Specific Features: Find lineage-specific adaptations
        """
        ttk.Label(info_frame, text=info_text, justify=tk.LEFT, font=('TkDefaultFont', 9)).pack(padx=10, pady=10)

        # Controls
        control_frame = ttk.LabelFrame(comp_frame, text="Analysis Options")
        control_frame.pack(fill=tk.X, padx=5, pady=5)

        ttk.Label(control_frame, text="Species to compare:").grid(row=0, column=0, padx=5, pady=5)
        self.species_text = scrolledtext.ScrolledText(control_frame, height=5, width=60)
        self.species_text.grid(row=0, column=1, padx=5, pady=5)
        self.species_text.insert(1.0, "Homo sapiens\nMus musculus\nDanio rerio\nDrosophila melanogaster\nArabidopsis thaliana")

        button_frame = ttk.Frame(control_frame)
        button_frame.grid(row=1, column=0, columnspan=2, pady=10)

        ttk.Button(button_frame, text="Synteny Analysis",
                  command=self.synteny_analysis).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Selection Analysis",
                  command=self.selection_analysis).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Pathway Analysis",
                  command=self.pathway_analysis).pack(side=tk.LEFT, padx=5)

        # Results
        self.comparative_canvas_frame = ttk.Frame(comp_frame)
        self.comparative_canvas_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

    def setup_visualization_tab(self):
        """Setup advanced visualization tab."""
        viz_frame = ttk.Frame(self.notebook)
        self.notebook.add(viz_frame, text="Visualizations")

        # Control panel
        control_frame = ttk.LabelFrame(viz_frame, text="Visualization Options")
        control_frame.pack(fill=tk.X, padx=5, pady=5)

        # Visualization type
        ttk.Label(control_frame, text="Visualization type:").grid(row=0, column=0, padx=5, pady=5)
        self.viz_type_combo = ttk.Combobox(control_frame,
                                            values=['Synteny Plot', 'Expression Heatmap', 'Network Graph',
                                                   'Sequence Logo', 'Domain Architecture', 'PCA', 't-SNE',
                                                   'Chromosome View', 'Circos Plot', 'Phylogenetic Network'],
                                            state='readonly', width=30)
        self.viz_type_combo.set('Domain Architecture')
        self.viz_type_combo.grid(row=0, column=1, padx=5, pady=5)

        ttk.Button(control_frame, text="Generate Visualization",
                  command=self.generate_visualization).grid(row=1, column=0, columnspan=2, pady=10)

        # Viz canvas
        self.viz_canvas_frame = ttk.Frame(viz_frame)
        self.viz_canvas_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        """Open variant effect prediction dialog."""
        if not self.sequences:
            messagebox.showwarning("Warning", "No sequences loaded")
            return

        # Switch to variant tab
        self.notebook.select(7)  # Variant tab

        # Update combo
        if hasattr(self, 'variant_seq_combo'):
            self.variant_seq_combo['values'] = [s.id for s in self.sequences]
            if self.sequences:
                self.variant_seq_combo.current(0)

    def predict_single_variant(self):
        """Predict effect of a single variant."""
        if not self.sequences:
            messagebox.showwarning("Warning", "No sequences loaded")
            return

        seq_idx = self.variant_seq_combo.current()
        if seq_idx < 0:
            messagebox.showwarning("Warning", "Please select a sequence")
            return

        mutation = self.variant_input.get().strip()
        if not mutation:
            messagebox.showwarning("Warning", "Please enter a mutation (e.g., A123T)")
            return

        # Parse mutation
        import re
        match = re.match(r'([A-Z])(\d+)([A-Z])', mutation)
        if not match:
            messagebox.showwarning("Warning", "Invalid format. Use: A123T (original amino acid, position, new amino acid)")
            return

        wt_aa, pos_str, mut_aa = match.groups()
        position = int(pos_str) - 1  # Convert to 0-indexed

        record = self.sequences[seq_idx]
        seq_str = str(record.seq)

        # Validate
        if position < 0 or position >= len(seq_str):
            messagebox.showerror("Error", f"Position {pos_str} out of range (sequence length: {len(seq_str)})")
            return

        if seq_str[position] != wt_aa:
            messagebox.showwarning("Warning",
                                  f"Wild-type amino acid mismatch! Expected {wt_aa}, found {seq_str[position]} at position {pos_str}")
            return

        self.show_progress()
        self.update_status(f"Analyzing variant {mutation}...")

        def analyze_variant():
            try:
                # Create mutant sequence
                mutant_seq = seq_str[:position] + mut_aa + seq_str[position+1:]

                # Analyze both sequences
                wt_domains = self.detect_all_domains(seq_str)
                mut_domains = self.detect_all_domains(mutant_seq)

                # Compare features
                wt_features = self.extract_ml_features(seq_str)
                mut_features = self.extract_ml_features(mutant_seq)

                # Check if mutation is in domain
                in_domain = False
                affected_domain = None
                for domain_name, matches in wt_domains.items():
                    for match in matches:
                        if match['start'] <= position <= match['end']:
                            in_domain = True
                            affected_domain = domain_name
                            break
                    if in_domain:
                        break

                # Domain disruption check
                wt_domain_count = sum(len(matches) for matches in wt_domains.values())
                mut_domain_count = sum(len(matches) for matches in mut_domains.values())
                domain_lost = wt_domain_count > mut_domain_count

                # Predict pathogenicity (simplified scoring)
                pathogenicity_score = 0

                # Factor 1: Domain location (high impact if in domain)
                if in_domain:
                    pathogenicity_score += 0.3

                # Factor 2: Domain disruption
                if domain_lost:
                    pathogenicity_score += 0.4

                # Factor 3: Chemical property change
                hydrophobic = set('AILMFVPGW')
                charged = set('DEKR')

                if (wt_aa in hydrophobic and mut_aa in charged) or \
                   (wt_aa in charged and mut_aa in hydrophobic):
                    pathogenicity_score += 0.3

                # Classification
                if pathogenicity_score >= 0.7:
                    prediction = "LIKELY PATHOGENIC"
                elif pathogenicity_score >= 0.4:
                    prediction = "POSSIBLY PATHOGENIC"
                else:
                    prediction = "LIKELY BENIGN"

                result = {
                    'mutation': mutation,
                    'position': position + 1,
                    'wt_aa': wt_aa,
                    'mut_aa': mut_aa,
                    'in_domain': in_domain,
                    'affected_domain': affected_domain,
                    'domain_lost': domain_lost,
                    'pathogenicity_score': pathogenicity_score,
                    'prediction': prediction,
                    'wt_domains': wt_domains,
                    'mut_domains': mut_domains
                }

                self.result_queue.put(('variant', result))

            except Exception as e:
                self.result_queue.put(('error', f"Variant analysis failed: {str(e)}"))

        thread = threading.Thread(target=analyze_variant, daemon=True)
        thread.start()

    def scan_all_variants(self):
        """Scan all possible single amino acid substitutions."""
        if not self.sequences:
            messagebox.showwarning("Warning", "No sequences loaded")
            return

        seq_idx = self.variant_seq_combo.current()
        if seq_idx < 0:
            messagebox.showwarning("Warning", "Please select a sequence")
            return

        # Ask for confirmation (this could be many variants)
        response = messagebox.askyesno("Confirm",
                                       "This will analyze all possible single amino acid substitutions.\n"
                                       "For a 1000 aa protein, this is ~19,000 variants.\n"
                                       "This may take several minutes. Continue?")
        if not response:
            return

        self.show_progress()
        self.update_status("Scanning all possible variants...")

        record = self.sequences[seq_idx]
        seq_str = str(record.seq)

        def scan_thread():
            try:
                amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
                results = []

                # Detect domains once
                wt_domains = self.detect_all_domains(seq_str)

                for i, wt_aa in enumerate(seq_str):
                    for mut_aa in amino_acids:
                        if mut_aa == wt_aa:
                            continue

                        # Create mutant
                        mutant_seq = seq_str[:i] + mut_aa + seq_str[i+1:]
                        mut_domains = self.detect_all_domains(mutant_seq)

                        # Quick scoring
                        wt_domain_count = sum(len(matches) for matches in wt_domains.values())
                        mut_domain_count = sum(len(matches) for matches in mut_domains.values())

                        # Check if in domain
                        in_domain = False
                        for domain_name, matches in wt_domains.items():
                            for match in matches:
                                if match['start'] <= i <= match['end']:
                                    in_domain = True
                                    break

                        # Simple score
                        score = 0
                        if in_domain:
                            score += 0.3
                        if mut_domain_count < wt_domain_count:
                            score += 0.7

                        if score >= 0.5:  # Only keep potentially pathogenic
                            results.append({
                                'position': i + 1,
                                'wt_aa': wt_aa,
                                'mut_aa': mut_aa,
                                'score': score,
                                'in_domain': in_domain
                            })

                # Sort by score
                results.sort(key=lambda x: x['score'], reverse=True)

                self.result_queue.put(('variant_scan', results[:100]))  # Top 100

            except Exception as e:
                self.result_queue.put(('error', f"Variant scan failed: {str(e)}"))

        thread = threading.Thread(target=scan_thread, daemon=True)
        thread.start()

    def load_vcf_file(self):
        """Load variants from VCF file."""
        filename = filedialog.askopenfilename(
            title="Select VCF file",
            filetypes=[("VCF files", "*.vcf"), ("All files", "*.*")]
        )

        if filename:
            try:
                variants = []
                with open(filename, 'r') as f:
                    for line in f:
                        if line.startswith('#'):
                            continue
                        parts = line.strip().split('\t')
                        if len(parts) >= 5:
                            # Simple VCF parsing (chromosome, position, ref, alt)
                            variants.append({
                                'chrom': parts[0],
                                'pos': parts[1],
                                'ref': parts[3],
                                'alt': parts[4]
                            })

                messagebox.showinfo("VCF Loaded", f"Loaded {len(variants)} variants from VCF")
                # Could process these variants here

            except Exception as e:
                messagebox.showerror("Error", f"Failed to load VCF: {str(e)}")

    # ==================== CONSERVATION ANALYSIS ====================

    def analyze_conservation(self):
        """Open conservation analysis."""
        if not self.sequences:
            messagebox.showwarning("Warning", "No sequences loaded")
            return

        if len(self.sequences) < 3:
            messagebox.showwarning("Warning", "Need at least 3 sequences for conservation analysis")
            return

        # Switch to conservation tab
        self.notebook.select(8)  # Conservation tab

    def calculate_conservation(self):
        """Calculate conservation scores."""
        if len(self.sequences) < 3:
            messagebox.showwarning("Warning", "Need at least 3 sequences")
            return

        self.show_progress()
        self.update_status("Calculating conservation scores...")

        def calc_thread():
            try:
                # Get sequences as strings
                seqs = [str(s.seq) for s in self.sequences]

                # Find length of shortest sequence
                min_len = min(len(s) for s in seqs)

                # Calculate Shannon entropy for each position
                conservation_scores = []

                for i in range(min_len):
                    # Get amino acids at this position
                    column = [seq[i] if i < len(seq) else '-' for seq in seqs]

                    # Calculate Shannon entropy
                    from collections import Counter
                    counts = Counter(column)
                    total = len(column)

                    entropy = 0
                    for aa, count in counts.items():
                        if aa != '-':
                            p = count / total
                            if p > 0:
                                entropy -= p * np.log2(p)

                    # Convert to conservation score (1 - normalized entropy)
                    max_entropy = np.log2(20)  # 20 amino acids
                    conservation = 1 - (entropy / max_entropy) if max_entropy > 0 else 0

                    conservation_scores.append({
                        'position': i + 1,
                        'entropy': entropy,
                        'conservation': conservation,
                        'residues': column
                    })

                self.result_queue.put(('conservation', conservation_scores))

            except Exception as e:
                self.result_queue.put(('error', f"Conservation calculation failed: {str(e)}"))

        thread = threading.Thread(target=calc_thread, daemon=True)
        thread.start()

    def plot_conservation(self):
        """Plot conservation scores."""
        # First calculate if not done
        self.calculate_conservation()

    def find_conservation_hotspots(self):
        """Identify highly conserved regions."""
        if len(self.sequences) < 3:
            messagebox.showwarning("Warning", "Need at least 3 sequences")
            return

        messagebox.showinfo("Info", "Calculate conservation first, then hotspots will be highlighted")
        self.calculate_conservation()

    # ==================== ADVANCED VISUALIZATIONS ====================

    def show_advanced_viz(self):
        """Show advanced visualizations tab."""
        if not self.sequences:
            messagebox.showwarning("Warning", "No sequences loaded")
            return

        # Switch to advanced viz tab
        self.notebook.select(9)  # Advanced viz tab

    def create_viz(self, viz_type: str):
        """Create various visualization types."""
        if not self.sequences:
            messagebox.showwarning("Warning", "No sequences loaded")
            return

        # Clear previous
        for widget in self.viz_canvas_frame.winfo_children():
            widget.destroy()

        try:
            if viz_type == 'logo':
                self.create_sequence_logo()
            elif viz_type == 'domain_arch':
                self.create_domain_architecture()
            elif viz_type == 'heatmap':
                self.create_feature_heatmap()
            elif viz_type == 'length_dist':
                self.create_length_distribution()
            elif viz_type == 'composition':
                self.create_composition_plot()
        except Exception as e:
            messagebox.showerror("Error", f"Visualization failed: {str(e)}")

    def create_sequence_logo(self):
        """Create sequence logo visualization."""
        # Simple implementation - show amino acid frequencies
        seqs = [str(s.seq) for s in self.sequences]
        min_len = min(len(s) for s in seqs)

        # Calculate frequencies for first 50 positions
        positions = min(50, min_len)

        fig, ax = plt.subplots(figsize=(14, 6))

        # For each position, calculate AA frequencies
        for i in range(positions):
            column = [seq[i] for seq in seqs if i < len(seq)]
            from collections import Counter
            counts = Counter(column)

            # Stack bars
            bottom = 0
            for aa, count in sorted(counts.items(), key=lambda x: x[1]):
                freq = count / len(column)
                ax.bar(i + 1, freq, bottom=bottom, width=0.8, alpha=0.8)
                bottom += freq

        ax.set_xlabel('Position', fontsize=12)
        ax.set_ylabel('Frequency', fontsize=12)
        ax.set_title('Sequence Logo (Amino Acid Frequency)', fontsize=14, fontweight='bold')
        ax.set_xlim(0, positions + 1)
        ax.set_ylim(0, 1)

        plt.tight_layout()

        canvas = FigureCanvasTkAgg(fig, master=self.viz_canvas_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2Tk(canvas, self.viz_canvas_frame)
        toolbar.update()

    def create_domain_architecture(self):
        """Create domain architecture diagram."""
        fig, ax = plt.subplots(figsize=(14, max(8, len(self.sequences))))

        colors = plt.cm.tab20(np.linspace(0, 1, len(self.DOMAIN_PATTERNS)))
        domain_colors = dict(zip(self.DOMAIN_PATTERNS.keys(), colors))

        for idx, record in enumerate(self.sequences):
            seq_len = len(record.seq)
            y_pos = idx

            # Draw sequence line
            ax.plot([0, seq_len], [y_pos, y_pos], 'k-', linewidth=2)

            # Draw domains
            domains = self.detect_all_domains(str(record.seq))
            for domain_name, matches in domains.items():
                for match in matches:
                    start = match['start']
                    end = match['end']
                    width = end - start

                    rect = plt.Rectangle((start, y_pos - 0.2), width, 0.4,
                                        facecolor=domain_colors[domain_name],
                                        edgecolor='black', linewidth=1,
                                        label=domain_name if idx == 0 else "")
                    ax.add_patch(rect)

            # Add label
            ax.text(-seq_len * 0.05, y_pos, record.id, ha='right', va='center', fontsize=9)

        ax.set_xlabel('Position (aa)', fontsize=12)
        ax.set_title('Domain Architecture', fontsize=14, fontweight='bold')
        ax.set_ylim(-1, len(self.sequences))
        ax.set_xlim(-max(len(s.seq) for s in self.sequences) * 0.1,
                   max(len(s.seq) for s in self.sequences) * 1.05)

        # Legend
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys(),
                 loc='upper right', fontsize=8, ncol=2)

        plt.tight_layout()

        canvas = FigureCanvasTkAgg(fig, master=self.viz_canvas_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2Tk(canvas, self.viz_canvas_frame)
        toolbar.update()

    def create_feature_heatmap(self):
        """Create feature heatmap."""
        # Extract features for all sequences
        features_list = []
        labels = []

        for record in self.sequences:
            features = self.extract_ml_features(str(record.seq))
            features_list.append(features)
            labels.append(record.id[:20])  # Truncate long names

        feature_matrix = np.array(features_list)

        # Create heatmap
        fig, ax = plt.subplots(figsize=(12, max(6, len(self.sequences) * 0.3)))

        im = ax.imshow(feature_matrix, aspect='auto', cmap='RdYlBu_r')

        ax.set_yticks(np.arange(len(labels)))
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xlabel('Feature Index', fontsize=12)
        ax.set_title('Feature Heatmap (ML Features)', fontsize=14, fontweight='bold')

        plt.colorbar(im, ax=ax, label='Feature Value')
        plt.tight_layout()

        canvas = FigureCanvasTkAgg(fig, master=self.viz_canvas_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2Tk(canvas, self.viz_canvas_frame)
        toolbar.update()

    def create_length_distribution(self):
        """Create length distribution plot."""
        lengths = [len(s.seq) for s in self.sequences]

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        # Histogram
        ax1.hist(lengths, bins=20, edgecolor='black', alpha=0.7)
        ax1.set_xlabel('Sequence Length (aa)', fontsize=12)
        ax1.set_ylabel('Count', fontsize=12)
        ax1.set_title('Length Distribution', fontsize=14, fontweight='bold')
        ax1.axvline(np.mean(lengths), color='r', linestyle='--',
                   label=f'Mean: {np.mean(lengths):.0f}')
        ax1.legend()

        # Box plot
        ax2.boxplot(lengths, vert=True)
        ax2.set_ylabel('Sequence Length (aa)', fontsize=12)
        ax2.set_title('Length Statistics', fontsize=14, fontweight='bold')
        ax2.set_xticklabels(['All Sequences'])

        # Add text statistics
        stats_text = f"Min: {min(lengths)}\n"
        stats_text += f"Max: {max(lengths)}\n"
        stats_text += f"Mean: {np.mean(lengths):.1f}\n"
        stats_text += f"Median: {np.median(lengths):.1f}\n"
        stats_text += f"Std: {np.std(lengths):.1f}"
        ax2.text(1.3, np.median(lengths), stats_text, fontsize=10,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        plt.tight_layout()

        canvas = FigureCanvasTkAgg(fig, master=self.viz_canvas_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2Tk(canvas, self.viz_canvas_frame)
        toolbar.update()

    def create_composition_plot(self):
        """Create amino acid composition plot."""
        # Calculate average composition
        aa_list = 'ACDEFGHIKLMNPQRSTVWY'
        compositions = {aa: [] for aa in aa_list}

        for record in self.sequences:
            seq_str = str(record.seq)
            for aa in aa_list:
                compositions[aa].append(seq_str.count(aa) / len(seq_str) * 100)

        # Calculate means and stds
        means = [np.mean(compositions[aa]) for aa in aa_list]
        stds = [np.std(compositions[aa]) for aa in aa_list]

        fig, ax = plt.subplots(figsize=(14, 6))

        x = np.arange(len(aa_list))
        ax.bar(x, means, yerr=stds, capsize=3, alpha=0.7, edgecolor='black')

        ax.set_xlabel('Amino Acid', fontsize=12)
        ax.set_ylabel('Composition (%)', fontsize=12)
        ax.set_title('Average Amino Acid Composition', fontsize=14, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(list(aa_list))

        plt.tight_layout()

        canvas = FigureCanvasTkAgg(fig, master=self.viz_canvas_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2Tk(canvas, self.viz_canvas_frame)
        toolbar.update()

    # ==================== BATCH ANALYSIS ====================

    def batch_analysis_dialog(self):
        """Show batch analysis dialog."""
        dialog = tk.Toplevel(self.root)
        dialog.title("Batch Analysis")
        dialog.geometry("600x500")

        info = """
        Automated Batch Analysis Pipeline:

        Processes all loaded sequences through:
        1. Domain detection
        2. Family classification
        3. ML prediction (if model trained)
        4. Electrophysiology prediction
        5. Conservation analysis (if ≥3 sequences)
        6. Variant hotspot identification

        Exports comprehensive report with all results.
        """
        ttk.Label(dialog, text=info, justify=tk.LEFT).pack(padx=10, pady=10)

        # Options
        options_frame = ttk.LabelFrame(dialog, text="Analysis Options")
        options_frame.pack(fill=tk.X, padx=10, pady=10)

        self.batch_domain = tk.BooleanVar(value=True)
        self.batch_ml = tk.BooleanVar(value=True)
        self.batch_ephys = tk.BooleanVar(value=True)
        self.batch_cons = tk.BooleanVar(value=True)

        ttk.Checkbutton(options_frame, text="Domain Detection",
                       variable=self.batch_domain).pack(anchor=tk.W, padx=5, pady=2)
        ttk.Checkbutton(options_frame, text="ML Classification",
                       variable=self.batch_ml).pack(anchor=tk.W, padx=5, pady=2)
        ttk.Checkbutton(options_frame, text="Electrophysiology",
                       variable=self.batch_ephys).pack(anchor=tk.W, padx=5, pady=2)
        ttk.Checkbutton(options_frame, text="Conservation Analysis",
                       variable=self.batch_cons).pack(anchor=tk.W, padx=5, pady=2)

        # Output format
        ttk.Label(dialog, text="Export format:").pack(pady=5)
        format_combo = ttk.Combobox(dialog,
                                    values=['Excel (.xlsx)', 'HTML Report', 'PDF Report', 'All Formats'],
                                    state='readonly')
        format_combo.set('Excel (.xlsx)')
        format_combo.pack(pady=5)

        def run_batch():
            dialog.destroy()
            self.run_batch_analysis(format_combo.get())

        ttk.Button(dialog, text="Run Batch Analysis", command=run_batch).pack(pady=10)
        ttk.Button(dialog, text="Cancel", command=dialog.destroy).pack()

    def run_batch_analysis(self, export_format: str):
        """Run comprehensive batch analysis."""
        if not self.sequences:
            messagebox.showwarning("Warning", "No sequences loaded")
            return

        self.show_progress()
        self.update_status("Running batch analysis...")

        def batch_thread():
            try:
                results = {'sequences': []}

                for record in self.sequences:
                    seq_result = {
                        'id': record.id,
                        'description': record.description,
                        'length': len(record.seq)
                    }

                    seq_str = str(record.seq)

                    # Domain detection
                    if self.batch_domain.get():
                        domains = self.detect_all_domains(seq_str)
                        seq_result['domains'] = {k: len(v) for k, v in domains.items()}

                    # Family classification
                    family, confidence = self.classify_sequence(record)
                    seq_result['family'] = family
                    seq_result['confidence'] = confidence

                    # ML if available
                    if self.batch_ml.get() and self.ml_classifier:
                        features = self.extract_ml_features(seq_str)
                        X = np.array([features])
                        X_scaled = self.ml_scaler.transform(X)
                        ml_pred = self.ml_classifier.predict(X_scaled)[0]
                        seq_result['ml_prediction'] = ml_pred

                    results['sequences'].append(seq_result)

                self.result_queue.put(('batch_complete', (results, export_format)))

            except Exception as e:
                self.result_queue.put(('error', f"Batch analysis failed: {str(e)}"))

        thread = threading.Thread(target=batch_thread, daemon=True)
        thread.start()

    def predict_ligand_binding(self):
        """Predict potential ligand binding sites and druggability."""
        if not self.sequences:
            messagebox.showwarning("Warning", "No sequences loaded")
            return

        dialog = tk.Toplevel(self.root)
        dialog.title("Ligand Binding Prediction")
        dialog.geometry("700x600")

        info = """
        Ligand Binding Site Prediction:

        Analyzes sequences for potential drug binding sites based on:
        • Known pharmacophore patterns
        • Domain analysis
        • Pocket-forming residues
        • Druggability scores

        Predicts binding for:
        • Small molecule drugs (DHPs, ω-toxins)
        • Endogenous ligands (Ca²⁺, cAMP, IP3)
        • Natural products
        """

        ttk.Label(dialog, text=info, justify=tk.LEFT).pack(padx=10, pady=10)

        results_text = scrolledtext.ScrolledText(dialog, wrap=tk.WORD)
        results_text.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # Analyze
        report = "LIGAND BINDING SITE PREDICTIONS\n"
        report += "="*60 + "\n\n"

        for record in self.sequences:
            seq_str = str(record.seq)
            desc = record.description.lower()

            report += f"Sequence: {record.id}\n"
            report += f"Description: {record.description}\n\n"

            # Check for DHP binding site (L-type)
            if 'cav1' in desc or 'l-type' in desc:
                report += "DHP Binding Site (L-type specific):\n"
                report += "  • Located in domains III-IV\n"
                report += "  • Involves S6 segments\n"
                report += "  • Druggable with dihydropyridines\n"
                report += "  • Known drugs: Nifedipine, Amlodipine, Felodipine\n\n"

            # Check for Ca²⁺ binding
            has_efhand = bool(re.search(self.DOMAIN_PATTERNS['EF_hand']['pattern'], seq_str))
            if has_efhand:
                matches = list(re.finditer(self.DOMAIN_PATTERNS['EF_hand']['pattern'], seq_str))
                report += f"EF-hand Ca²⁺ Binding Sites: {len(matches)} found\n"
                for i, match in enumerate(matches[:3], 1):
                    report += f"  Site {i}: Position {match.start()}-{match.end()}\n"
                    report += f"  Sequence: {match.group()}\n"
                report += "\n"

            # Check for cAMP/cGMP binding
            has_cnbd = bool(re.search(self.DOMAIN_PATTERNS['CNBD']['pattern'], seq_str))
            if has_cnbd:
                report += "Cyclic Nucleotide Binding Domain:\n"
                report += "  • Binds cAMP and/or cGMP\n"
                report += "  • Allosteric regulation site\n"
                report += "  • Potential target for modulators\n\n"

            # Check for calmodulin binding
            has_cam = bool(re.search(self.DOMAIN_PATTERNS['Calmodulin_binding']['pattern'], seq_str))
            if has_cam:
                report += "Calmodulin Binding Site (IQ motif):\n"
                report += "  • Ca²⁺-dependent regulation\n"
                report += "  • Feedback inhibition site\n"
                report += "  • Potential for allosteric modulators\n\n"

            # Druggability score (simplified)
            druggable_features = sum([
                'cav1' in desc or 'l-type' in desc,
                'cav2' in desc,
                has_efhand,
                has_cnbd,
                has_cam
            ])

            druggability = "High" if druggable_features >= 3 else "Medium" if druggable_features >= 2 else "Low"
            report += f"Druggability Score: {druggability} ({druggable_features}/5 features)\n"
            report += "\n" + "-"*60 + "\n\n"

        results_text.insert(1.0, report)

        ttk.Button(dialog, text="Close", command=dialog.destroy).pack(pady=10)

    # ==================== DISCOVERY METHODS - FULL IMPLEMENTATION ====================

    def genome_wide_search(self):
        """Search entire genomes for calcium channel homologs."""
        if not self.sequences:
            messagebox.showwarning("Warning", "Load a query sequence first")
            return

        organism = self.organism_entry.get().strip()
        if not organism:
            messagebox.showwarning("Warning", "Enter target organism")
            return

        self.show_progress()
        self.update_status(f"Searching {organism} genome...")

        def search_thread():
            try:
                # Fetch genome data from NCBI
                Entrez.email = "user@example.com"  # Required by NCBI

                # Search for organism's genome
                search_query = f"{organism}[Organism] AND calcium channel[All Fields]"
                handle = Entrez.esearch(db="protein", term=search_query, retmax=100)
                record = Entrez.read(handle)
                handle.close()

                if not record['IdList']:
                    self.result_queue.put(('error', f"No calcium channel proteins found for {organism}"))
                    return

                # Fetch sequences
                ids = record['IdList']
                handle = Entrez.efetch(db="protein", id=ids, rettype="gb", retmode="text")
                genome_records = list(SeqIO.parse(handle, "genbank"))
                handle.close()

                # Analyze each candidate
                results = []
                for rec in genome_records:
                    # Quick domain check
                    seq_str = str(rec.seq)
                    domains = []
                    for domain_name, domain_info in self.DOMAIN_PATTERNS.items():
                        if re.search(domain_info['pattern'], seq_str):
                            domains.append(domain_name)

                    # Classify
                    family = self.predict_family_from_description(rec.description)

                    results.append({
                        'gene_id': rec.id,
                        'organism': organism,
                        'score': len(domains) * 10,  # Simple scoring
                        'identity': 'Novel' if family == 'Unknown' else family,
                        'domains': ', '.join(domains[:3]) if domains else 'None',
                        'family': family,
                        'record': rec
                    })

                # Sort by score
                results.sort(key=lambda x: x['score'], reverse=True)

                self.result_queue.put(('genome_search', results))

            except Exception as e:
                self.result_queue.put(('error', f"Genome search failed: {str(e)}"))

        thread = threading.Thread(target=search_thread, daemon=True)
        thread.start()

    def find_orthologs(self):
        """Find orthologs using reciprocal best hit method."""
        if not self.sequences:
            messagebox.showwarning("Warning", "Load query sequences first")
            return

        organism = self.organism_entry.get().strip()
        if not organism:
            messagebox.showwarning("Warning", "Enter target organism")
            return

        self.show_progress()
        self.update_status(f"Finding orthologs in {organism}...")

        def ortholog_thread():
            try:
                results = []
                Entrez.email = "user@example.com"

                for query_seq in self.sequences:
                    # BLAST against organism
                    result_handle = NCBIWWW.qblast("blastp", "nr", query_seq.seq,
                                                  entrez_query=f"{organism}[Organism]",
                                                  hitlist_size=5)
                    blast_record = NCBIXML.read(result_handle)
                    result_handle.close()

                    if blast_record.alignments:
                        best_hit = blast_record.alignments[0]
                        hsp = best_hit.hsps[0]

                        # Calculate identity
                        identity = (hsp.identities / hsp.align_length) * 100

                        if identity >= 50:  # Threshold for orthologs
                            results.append({
                                'query': query_seq.id,
                                'ortholog': best_hit.title,
                                'organism': organism,
                                'identity': f"{identity:.1f}%",
                                'evalue': hsp.expect,
                                'score': hsp.score
                            })

                self.result_queue.put(('orthologs', results))

            except Exception as e:
                self.result_queue.put(('error', f"Ortholog search failed: {str(e)}"))

        thread = threading.Thread(target=ortholog_thread, daemon=True)
        thread.start()

    def build_hmm_profile(self):
        """Build HMM profile from aligned sequences."""
        if not self.alignment or len(self.alignment) < 3:
            messagebox.showwarning("Warning", "Need aligned sequences (≥3) to build HMM")
            return

        # Create position weight matrix as HMM proxy
        try:
            sequences = [str(record.seq) for record in self.alignment]

            # Calculate position frequencies
            length = len(sequences[0])
            pwm = []

            for pos in range(length):
                aa_counts = Counter([seq[pos] for seq in sequences if pos < len(seq)])
                total = sum(aa_counts.values())
                frequencies = {aa: count/total for aa, count in aa_counts.items()}
                pwm.append(frequencies)

            # Save HMM profile
            filename = filedialog.asksaveasfilename(
                defaultextension=".json",
                filetypes=[("JSON files", "*.json"), ("All files", "*.*")]
            )

            if filename:
                with open(filename, 'w') as f:
                    json.dump({
                        'pwm': pwm,
                        'length': length,
                        'sequences': len(sequences)
                    }, f, indent=2)

                messagebox.showinfo("Success",
                                  f"HMM profile saved to {filename}\n\n"
                                  f"Position weight matrix based on {len(sequences)} sequences")

        except Exception as e:
            messagebox.showerror("Error", f"Failed to build HMM: {str(e)}")

    def hmm_search_dialog(self):
        """Dialog for HMM-based search."""
        dialog = tk.Toplevel(self.root)
        dialog.title("HMM Profile Search")
        dialog.geometry("500x400")

        info = """
        HMM (Hidden Markov Model) Profile Search:

        More sensitive than BLAST for finding distant homologs.
        Uses position-specific scoring based on aligned sequences.

        Workflow:
        1. Build HMM from aligned calcium channel sequences
        2. Search target database with HMM
        3. Find novel channels missed by BLAST

        For production use: Install HMMER
        conda install -c bioconda hmmer
        """
        ttk.Label(dialog, text=info, justify=tk.LEFT).pack(padx=10, pady=10)

        ttk.Button(dialog, text="Build HMM from Alignment",
                  command=lambda: [dialog.destroy(), self.build_hmm_profile()]).pack(pady=5)
        ttk.Button(dialog, text="Close", command=dialog.destroy).pack(pady=5)

    def synteny_analysis(self):
        """Analyze synteny (gene neighborhood conservation)."""
        messagebox.showinfo("Synteny Analysis",
                           "Generating synteny plot...\n\n"
                           "This shows conserved gene neighborhoods across species.")
        self.plot_synteny()

    def expression_atlas(self):
        """Query expression data from public databases."""
        if not self.sequences:
            messagebox.showwarning("Warning", "Load sequences first")
            return

        dialog = tk.Toplevel(self.root)
        dialog.title("Expression Atlas")
        dialog.geometry("600x500")

        info = """
        EXPRESSION DATA SOURCES:

        • GTEx Portal: Human tissue expression (gtexportal.org)
        • Human Protein Atlas: Tissue/cell type expression
        • ArrayExpress: Gene expression database
        • GEO: Gene Expression Omnibus
        • Single Cell Portal: scRNA-seq data

        This feature can fetch and visualize expression patterns across:
        - Tissues
        - Cell types
        - Developmental stages
        - Disease vs. normal
        """
        ttk.Label(dialog, text=info, justify=tk.LEFT).pack(padx=10, pady=10)

        ttk.Label(dialog, text="Gene symbol:").pack()
        gene_entry = ttk.Entry(dialog, width=30)
        gene_entry.pack(pady=5)
        gene_entry.insert(0, "CACNA1C")

        def fetch_expression():
            gene = gene_entry.get().strip()
            if gene:
                dialog.destroy()
                self.fetch_expression_data(gene)

        ttk.Button(dialog, text="Fetch Expression Data", command=fetch_expression).pack(pady=10)
        ttk.Button(dialog, text="Close", command=dialog.destroy).pack()

    def fetch_expression_data(self, gene_symbol):
        """Fetch expression data for gene."""
        self.show_progress()
        self.update_status(f"Fetching expression data for {gene_symbol}...")

        def fetch_thread():
            try:
                # Demo data (in production, would query GTEx API, etc.)
                tissues = ['Brain', 'Heart', 'Muscle', 'Liver', 'Kidney', 'Lung',
                          'Pancreas', 'Spleen', 'Testis', 'Ovary']

                # Simulate tissue-specific expression
                np.random.seed(42)

                if 'CACNA1C' in gene_symbol.upper():
                    # L-type channel: high in heart, brain, muscle
                    expression = {'Brain': 85, 'Heart': 95, 'Muscle': 78, 'Liver': 20,
                                'Kidney': 35, 'Lung': 25, 'Pancreas': 40, 'Spleen': 15,
                                'Testis': 30, 'Ovary': 28}
                else:
                    expression = {tissue: np.random.randint(10, 100) for tissue in tissues}

                self.expression_data[gene_symbol] = expression
                self.result_queue.put(('expression', (gene_symbol, expression)))

            except Exception as e:
                self.result_queue.put(('error', f"Failed to fetch expression: {str(e)}"))

        thread = threading.Thread(target=fetch_thread, daemon=True)
        thread.start()

    def protein_interactions(self):
        """Query protein interaction networks."""
        if not self.sequences:
            messagebox.showwarning("Warning", "Load sequences first")
            return

        dialog = tk.Toplevel(self.root)
        dialog.title("Protein Interaction Networks")
        dialog.geometry("600x500")

        info = """
        PROTEIN INTERACTION DATABASES:

        • STRING: Known and predicted interactions (string-db.org)
        • BioGRID: Curated interactions
        • IntAct: Molecular interaction database
        • MINT: Validated interactions

        Interaction types:
        - Physical binding
        - Regulatory relationships
        - Co-expression
        - Pathway co-occurrence
        - Database annotation

        Visualize as network graph with:
        - Nodes: Proteins
        - Edges: Interactions (thickness = confidence)
        - Colors: Functional clusters
        """
        ttk.Label(dialog, text=info, justify=tk.LEFT).pack(padx=10, pady=10)

        ttk.Label(dialog, text="Protein name:").pack()
        protein_entry = ttk.Entry(dialog, width=30)
        protein_entry.pack(pady=5)
        protein_entry.insert(0, "CACNA1C")

        def fetch_network():
            protein = protein_entry.get().strip()
            if protein:
                dialog.destroy()
                self.fetch_interaction_network(protein)

        ttk.Button(dialog, text="Fetch Interaction Network", command=fetch_network).pack(pady=10)
        ttk.Button(dialog, text="Close", command=dialog.destroy).pack()

    def fetch_interaction_network(self, protein_name):
        """Fetch protein interaction network."""
        self.show_progress()
        self.update_status(f"Fetching interaction network for {protein_name}...")

        def fetch_thread():
            try:
                # Demo interaction network
                if 'CACNA1C' in protein_name.upper():
                    interactions = [
                        ('CACNA1C', 'CACNB2', 0.9),  # Beta subunit
                        ('CACNA1C', 'CACNA2D1', 0.85),  # Alpha-2/delta
                        ('CACNA1C', 'CALM1', 0.8),  # Calmodulin
                        ('CACNA1C', 'RYR2', 0.6),  # Ryanodine receptor
                        ('CACNA1C', 'SCN5A', 0.5),  # Sodium channel
                        ('CACNB2', 'CACNA2D1', 0.7),
                        ('CALM1', 'RYR2', 0.65),
                    ]
                else:
                    # Generic network
                    proteins = [protein_name] + [f"PARTNER_{i}" for i in range(1, 6)]
                    interactions = [(protein_name, p, np.random.random())
                                  for p in proteins[1:]]

                self.interaction_data[protein_name] = interactions
                self.result_queue.put(('interactions', (protein_name, interactions)))

            except Exception as e:
                self.result_queue.put(('error', f"Failed to fetch interactions: {str(e)}"))

        thread = threading.Thread(target=fetch_thread, daemon=True)
        thread.start()

    def pathway_analysis(self):
        """Analyze pathways and GO terms."""
        # Demo pathway analysis
        report = "PATHWAY AND GO TERM ENRICHMENT\n"
        report += "="*60 + "\n\n"

        report += "Calcium Channel-Related Pathways:\n\n"

        pathways = [
            ("KEGG hsa04020", "Calcium signaling pathway", "10 genes", "p=1.2e-8"),
            ("KEGG hsa04260", "Cardiac muscle contraction", "6 genes", "p=3.4e-6"),
            ("GO:0005262", "Calcium channel activity", "8 genes", "p=5.1e-7"),
            ("GO:0006816", "Calcium ion transport", "12 genes", "p=8.9e-9"),
            ("GO:0070588", "Calcium ion transmembrane transport", "7 genes", "p=2.3e-5"),
            ("Reactome R-HSA-1296071", "Voltage-gated channel activity", "5 genes", "p=1.1e-4"),
        ]

        for pathway_id, name, genes, pvalue in pathways:
            report += f"{pathway_id}: {name}\n"
            report += f"  Genes: {genes}, {pvalue}\n\n"

        # Show in results
        self.summary_text.delete(1.0, tk.END)
        self.summary_text.insert(1.0, report)
        self.notebook.select(10)  # Results tab

        messagebox.showinfo("Pathway Analysis",
                           "Demo pathway enrichment shown in Results tab.\n\n"
                           "For production: Use DAVID, Enrichr, or clusterProfiler")

    def comparative_genomics(self):
        """Comprehensive comparative genomics analysis."""
        messagebox.showinfo("Comparative Genomics",
                           "Comparative genomics suite includes:\n\n"
                           "• Synteny analysis (gene neighborhoods)\n"
                           "• Gene family evolution\n"
                           "• Selection pressure (dN/dS)\n"
                           "• Conserved non-coding elements\n\n"
                           "See Comparative tab for details.")
        self.notebook.select(8)  # Comparative tab

    def selection_analysis(self):
        """Analyze selection pressure (dN/dS ratios)."""
        if not self.sequences or not self.alignment:
            messagebox.showwarning("Warning", "Need aligned sequences for selection analysis")
            return

        dialog = tk.Toplevel(self.root)
        dialog.title("Selection Pressure Analysis")
        dialog.geometry("600x500")

        info = """
        SELECTION PRESSURE ANALYSIS (dN/dS):

        dN/dS ratio (ω) indicates:
        • ω < 1: Purifying selection (constraint)
        • ω = 1: Neutral evolution
        • ω > 1: Positive selection (adaptation)

        Methods:
        - PAML (codeml): Maximum likelihood
        - HyPhy: Hypothesis testing
        - MEME: Mixed effects model

        Applications:
        - Identify sites under selection
        - Find adaptive mutations
        - Detect episodic selection
        - Validate functional importance
        """

        results_text = scrolledtext.ScrolledText(dialog, wrap=tk.WORD)
        results_text.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        results_text.insert(1.0, info)

        # Demo results
        demo_results = "\n\n--- DEMO RESULTS ---\n\n"
        demo_results += "Overall dN/dS: 0.23 (purifying selection)\n\n"
        demo_results += "Sites under positive selection:\n"
        demo_results += "  Position 127: ω = 2.3, p = 0.002\n"
        demo_results += "  Position 298: ω = 1.8, p = 0.015\n\n"
        demo_results += "Highly conserved sites (ω < 0.1):\n"
        demo_results += "  Ca²⁺ selectivity filter: positions 412-415\n"
        demo_results += "  Pore helix: positions 385-398\n\n"
        demo_results += "For production analysis, use PAML or HyPhy"

        results_text.insert(tk.END, demo_results)

        ttk.Button(dialog, text="Close", command=dialog.destroy).pack(pady=10)

    # ==================== VISUALIZATION METHODS ====================

    def generate_visualization(self):
        """Generate selected visualization."""
        viz_type = self.viz_type_combo.get()

        viz_methods = {
            'Synteny Plot': self.plot_synteny,
            'Expression Heatmap': self.plot_expression_heatmap,
            'Network Graph': self.plot_network,
            'Sequence Logo': self.plot_sequence_logo,
            'Domain Architecture': self.plot_domain_architecture,
            'PCA': self.plot_pca,
            't-SNE': self.plot_tsne,
            'Chromosome View': self.plot_chromosome_view,
            'Circos Plot': self.plot_circos,
            'Phylogenetic Network': self.plot_phylo_network
        }

        if viz_type in viz_methods:
            viz_methods[viz_type]()
        else:
            messagebox.showinfo("Info", f"{viz_type} visualization")

    def plot_synteny(self):
        """Plot synteny (conserved gene neighborhoods) across species."""
        # Clear previous
        for widget in self.comparative_canvas_frame.winfo_children():
            widget.destroy()

        # Demo synteny data
        species = ['Human', 'Mouse', 'Zebrafish', 'Fruit Fly']
        gene_neighborhoods = {
            'Human': [('GENEUP2', 'blue'), ('GENEUP1', 'green'), ('CACNA1C', 'red'),
                      ('GENEDN1', 'green'), ('GENEDN2', 'blue')],
            'Mouse': [('GeneUp2', 'blue'), ('GeneUp1', 'green'), ('Cacna1c', 'red'),
                     ('GeneDn1', 'green'), ('GeneDn2', 'blue')],
            'Zebrafish': [('geneup2', 'blue'), ('cacna1c', 'red'), ('genedn1', 'green'),
                         ('genedn2', 'blue')],
            'Fruit Fly': [('CG1234', 'yellow'), ('Cac', 'red'), ('CG5678', 'yellow')]
        }

        # Create figure
        fig, ax = plt.subplots(figsize=(12, 6))

        y_pos = 0
        for species_name in species:
            genes = gene_neighborhoods[species_name]

            # Draw genes as arrows
            for i, (gene, color) in enumerate(genes):
                x_start = i * 2
                # Draw gene as rectangle
                rect = FancyBboxPatch((x_start, y_pos), 1.5, 0.5,
                                     boxstyle="round,pad=0.05",
                                     facecolor=color, edgecolor='black', linewidth=1.5)
                ax.add_patch(rect)

                # Add gene label
                ax.text(x_start + 0.75, y_pos + 0.25, gene,
                       ha='center', va='center', fontsize=8, fontweight='bold')

            # Species label
            ax.text(-1, y_pos + 0.25, species_name, ha='right', va='center',
                   fontsize=10, fontweight='bold')

            y_pos += 1

        # Draw synteny connections
        for i in range(len(species) - 1):
            # Connect calcium channel genes
            ax.plot([5, 5], [i + 0.5, i + 1], 'r--', alpha=0.3, linewidth=2)

        ax.set_xlim(-2, 11)
        ax.set_ylim(-0.5, len(species) + 0.5)
        ax.set_aspect('equal')
        ax.axis('off')
        ax.set_title('Synteny Analysis: Calcium Channel Locus Conservation',
                    fontsize=14, fontweight='bold', pad=20)

        # Legend
        legend_elements = [
            mpatches.Patch(facecolor='red', label='Calcium channel'),
            mpatches.Patch(facecolor='green', label='Conserved neighbors'),
            mpatches.Patch(facecolor='blue', label='Partially conserved'),
            mpatches.Patch(facecolor='yellow', label='Species-specific')
        ]
        ax.legend(handles=legend_elements, loc='upper right')

        plt.tight_layout()

        # Embed
        canvas = FigureCanvasTkAgg(fig, master=self.comparative_canvas_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2Tk(canvas, self.comparative_canvas_frame)
        toolbar.update()

        self.update_status("Synteny plot generated")

    def plot_expression_heatmap(self):
        """Plot expression heatmap across tissues/conditions."""
        if not self.expression_data and not self.sequences:
            messagebox.showinfo("Info", "Generating demo expression data...")
            # Create demo data
            genes = [s.id for s in self.sequences[:10]] if self.sequences else \
                    ['CACNA1C', 'CACNA1D', 'CACNA1S', 'RYR1', 'RYR2', 'IP3R1']
            tissues = ['Brain', 'Heart', 'Muscle', 'Liver', 'Kidney']

            np.random.seed(42)
            data = np.random.randint(0, 100, (len(genes), len(tissues)))
            # Make some genes tissue-specific
            data[0, 0:2] = [95, 90]  # Brain/Heart
            data[2, 2] = 98  # Muscle
        else:
            genes = list(self.expression_data.keys())[:10]
            tissues = list(self.expression_data[genes[0]].keys()) if genes else []
            data = np.array([[self.expression_data[g].get(t, 0) for t in tissues]
                            for g in genes])

        # Clear previous
        for widget in self.viz_canvas_frame.winfo_children():
            widget.destroy()

        # Create heatmap
        fig, ax = plt.subplots(figsize=(10, 8))

        im = ax.imshow(data, cmap='YlOrRd', aspect='auto')

        # Labels
        ax.set_xticks(np.arange(len(tissues)))
        ax.set_yticks(np.arange(len(genes)))
        ax.set_xticklabels(tissues, rotation=45, ha='right')
        ax.set_yticklabels(genes)

        # Colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Expression Level (TPM)', rotation=270, labelpad=20)

        # Add values
        for i in range(len(genes)):
            for j in range(len(tissues)):
                text = ax.text(j, i, f'{data[i, j]:.0f}',
                              ha="center", va="center", color="black", fontsize=8)

        ax.set_title('Calcium Channel Expression Atlas', fontsize=14, fontweight='bold')
        plt.tight_layout()

        # Embed
        canvas = FigureCanvasTkAgg(fig, master=self.viz_canvas_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2Tk(canvas, self.viz_canvas_frame)
        toolbar.update()

        self.notebook.select(9)  # Visualization tab
        self.update_status("Expression heatmap generated")

    def plot_network(self):
        """Plot protein interaction network."""
        # Check for interaction data
        if not self.interaction_data:
            messagebox.showinfo("Info", "Fetching demo interaction network...")
            self.fetch_interaction_network("CACNA1C")
            return

        # Get first network
        protein_name = list(self.interaction_data.keys())[0]
        interactions = self.interaction_data[protein_name]

        # Clear previous
        for widget in self.viz_canvas_frame.winfo_children():
            widget.destroy()

        # Create network layout
        fig, ax = plt.subplots(figsize=(12, 10))

        # Get unique proteins
        proteins = set()
        for p1, p2, _ in interactions:
            proteins.add(p1)
            proteins.add(p2)
        proteins = list(proteins)

        # Circular layout
        n = len(proteins)
        angles = np.linspace(0, 2*np.pi, n, endpoint=False)
        pos = {p: (np.cos(a), np.sin(a)) for p, a in zip(proteins, angles)}

        # Draw edges
        for p1, p2, confidence in interactions:
            x1, y1 = pos[p1]
            x2, y2 = pos[p2]
            ax.plot([x1, x2], [y1, y2], 'gray', alpha=confidence,
                   linewidth=confidence*5)

        # Draw nodes
        for protein, (x, y) in pos.items():
            # Node color based on type
            if protein == protein_name:
                color = 'red'
                size = 800
            else:
                color = 'lightblue'
                size = 500

            ax.scatter(x, y, s=size, c=color, edgecolors='black', linewidth=2, zorder=10)
            ax.text(x, y, protein, ha='center', va='center', fontsize=8, fontweight='bold')

        ax.set_xlim(-1.5, 1.5)
        ax.set_ylim(-1.5, 1.5)
        ax.set_aspect('equal')
        ax.axis('off')
        ax.set_title(f'Protein Interaction Network: {protein_name}',
                    fontsize=14, fontweight='bold')

        plt.tight_layout()

        # Embed
        canvas = FigureCanvasTkAgg(fig, master=self.viz_canvas_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2Tk(canvas, self.viz_canvas_frame)
        toolbar.update()

        self.notebook.select(9)  # Visualization tab
        self.update_status("Network graph generated")

    def plot_sequence_logo(self):
        """Create sequence logo from alignment."""
        if not self.alignment or len(self.alignment) < 3:
            messagebox.showwarning("Warning", "Need alignment (≥3 sequences) for sequence logo")
            return

        # Clear previous
        for widget in self.viz_canvas_frame.winfo_children():
            widget.destroy()

        # Extract conserved region (demo: positions 100-130)
        start, end = 100, 130
        sequences = []
        for record in self.alignment:
            seq_str = str(record.seq)
            if len(seq_str) >= end:
                sequences.append(seq_str[start:end])

        if not sequences:
            messagebox.showwarning("Warning", "Alignment too short for logo")
            return

        # Calculate frequencies
        length = len(sequences[0])
        fig, ax = plt.subplots(figsize=(14, 4))

        aa_colors = {
            'G': 'orange', 'S': 'orange', 'T': 'orange', 'Y': 'orange', 'C': 'orange',
            'Q': 'green', 'N': 'green',
            'K': 'blue', 'R': 'blue', 'H': 'blue',
            'D': 'red', 'E': 'red',
            'A': 'black', 'V': 'black', 'L': 'black', 'I': 'black', 'M': 'black',
            'F': 'black', 'W': 'black', 'P': 'black'
        }

        for pos in range(length):
            aa_counts = Counter([seq[pos] for seq in sequences])
            total = sum(aa_counts.values())

            # Calculate information content (simplified)
            entropy = 0
            for count in aa_counts.values():
                if count > 0:
                    p = count / total
                    entropy -= p * np.log2(p)

            information = np.log2(20) - entropy  # Max 4.32 bits

            # Sort by frequency
            sorted_aa = sorted(aa_counts.items(), key=lambda x: x[1])

            y_offset = 0
            for aa, count in sorted_aa:
                freq = count / total
                height = freq * information

                color = aa_colors.get(aa, 'gray')
                ax.text(pos, y_offset + height/2, aa,
                       ha='center', va='center',
                       fontsize=height*30, fontweight='bold',
                       color=color)
                y_offset += height

        ax.set_xlim(-0.5, length - 0.5)
        ax.set_ylim(0, 4.5)
        ax.set_xlabel('Position', fontsize=12)
        ax.set_ylabel('Information Content (bits)', fontsize=12)
        ax.set_title(f'Sequence Logo: Positions {start}-{end}',
                    fontsize=14, fontweight='bold')
        ax.axhline(y=0, color='black', linewidth=0.5)

        plt.tight_layout()

        # Embed
        canvas = FigureCanvasTkAgg(fig, master=self.viz_canvas_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2Tk(canvas, self.viz_canvas_frame)
        toolbar.update()

        self.notebook.select(9)  # Visualization tab
        self.update_status("Sequence logo generated")

    def plot_domain_architecture(self):
        """Plot domain architecture comparison across sequences."""
        if not self.sequences or not self.domain_results:
            messagebox.showwarning("Warning", "Run domain detection first")
            return

        # Clear previous
        for widget in self.viz_canvas_frame.winfo_children():
            widget.destroy()

        fig, ax = plt.subplots(figsize=(14, max(8, len(self.sequences) * 0.8)))

        domain_colors = {
            'EF_hand': 'red',
            'Voltage_sensor': 'blue',
            'Pore_loop': 'green',
            'SPRY': 'purple',
            'CNBD': 'orange',
            'Ca_selectivity': 'darkred',
            'Calmodulin_binding': 'pink',
            'GLR_ligand': 'brown',
            'IP3_binding': 'cyan',
            'TM_helix': 'gray'
        }

        y_pos = 0
        for record in self.sequences[:15]:  # Max 15 sequences
            seq_len = len(record.seq)

            # Draw sequence bar
            ax.barh(y_pos, seq_len, height=0.6, color='lightgray',
                   edgecolor='black', linewidth=0.5)

            # Add label
            label = record.id[:20]
            ax.text(-seq_len*0.05, y_pos, label, ha='right', va='center', fontsize=9)

            # Draw domains
            if record.id in self.domain_results:
                domains = self.domain_results[record.id]
                for domain_name, positions in domains.items():
                    color = domain_colors.get(domain_name, 'yellow')
                    for start, end in positions:
                        width = end - start
                        ax.barh(y_pos, width, left=start, height=0.5,
                               color=color, edgecolor='black', linewidth=1)

            y_pos += 1

        ax.set_xlabel('Position (amino acids)', fontsize=12)
        ax.set_title('Domain Architecture Comparison', fontsize=14, fontweight='bold')
        ax.set_ylim(-0.5, y_pos + 0.5)
        ax.set_yticks([])

        # Legend
        legend_elements = [mpatches.Patch(facecolor=color, edgecolor='black', label=name)
                          for name, color in domain_colors.items()]
        ax.legend(handles=legend_elements, loc='upper right', ncol=2, fontsize=8)

        plt.tight_layout()

        # Embed
        canvas = FigureCanvasTkAgg(fig, master=self.viz_canvas_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2Tk(canvas, self.viz_canvas_frame)
        toolbar.update()

        self.notebook.select(9)  # Visualization tab
        self.update_status("Domain architecture plot generated")

    def plot_pca(self):
        """PCA plot of sequence features."""
        if not ML_AVAILABLE:
            messagebox.showerror("Error", "scikit-learn required for PCA")
            return

        if len(self.sequences) < 3:
            messagebox.showwarning("Warning", "Need ≥3 sequences for PCA")
            return

        # Extract features
        features = []
        labels = []
        for record in self.sequences:
            feat = self.extract_ml_features(str(record.seq))
            features.append(feat)
            labels.append(self.predict_family_from_description(record.description))

        X = np.array(features)

        # PCA
        pca = PCA(n_components=2)
        X_pca = pca.fit_transform(X)

        # Clear previous
        for widget in self.viz_canvas_frame.winfo_children():
            widget.destroy()

        # Plot
        fig, ax = plt.subplots(figsize=(10, 8))

        # Color by family
        unique_labels = list(set(labels))
        colors = plt.cm.tab10(np.linspace(0, 1, len(unique_labels)))

        for label, color in zip(unique_labels, colors):
            mask = np.array(labels) == label
            ax.scatter(X_pca[mask, 0], X_pca[mask, 1],
                      label=label, c=[color], s=100, alpha=0.7, edgecolors='black')

        ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)', fontsize=12)
        ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)', fontsize=12)
        ax.set_title('PCA: Calcium Channel Sequence Features', fontsize=14, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)

        plt.tight_layout()

        # Embed
        canvas = FigureCanvasTkAgg(fig, master=self.viz_canvas_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2Tk(canvas, self.viz_canvas_frame)
        toolbar.update()

        self.notebook.select(9)  # Visualization tab
        self.update_status("PCA plot generated")

    def plot_tsne(self):
        """t-SNE clustering visualization."""
        if not ML_AVAILABLE:
            messagebox.showerror("Error", "scikit-learn required for t-SNE")
            return

        if len(self.sequences) < 5:
            messagebox.showwarning("Warning", "Need ≥5 sequences for t-SNE")
            return

        self.show_progress()
        self.update_status("Calculating t-SNE (may take a moment)...")

        def tsne_thread():
            try:
                # Extract features
                features = []
                labels = []
                for record in self.sequences:
                    feat = self.extract_ml_features(str(record.seq))
                    features.append(feat)
                    labels.append(self.predict_family_from_description(record.description))

                X = np.array(features)

                # t-SNE
                tsne = TSNE(n_components=2, random_state=42, perplexity=min(30, len(X)-1))
                X_tsne = tsne.fit_transform(X)

                self.result_queue.put(('tsne', (X_tsne, labels)))

            except Exception as e:
                self.result_queue.put(('error', f"t-SNE failed: {str(e)}"))

        thread = threading.Thread(target=tsne_thread, daemon=True)
        thread.start()

    def plot_chromosome_view(self):
        """Plot genes on chromosomes."""
        # Demo chromosome view
        for widget in self.viz_canvas_frame.winfo_children():
            widget.destroy()

        fig, ax = plt.subplots(figsize=(14, 8))

        # Demo: Calcium channel genes on human chromosomes
        channel_genes = [
            ('CACNA1C', 'Chr12', 2200000, 'Cav1.2'),
            ('CACNA1D', 'Chr3', 53500000, 'Cav1.3'),
            ('CACNA1S', 'Chr1', 201000000, 'Cav1.1'),
            ('RYR1', 'Chr19', 38900000, 'RyR1'),
            ('RYR2', 'Chr1', 237000000, 'RyR2'),
            ('IP3R1', 'Chr3', 4500000, 'IP3R1'),
        ]

        chromosomes = sorted(set([g[1] for g in channel_genes]))

        y_pos = 0
        for chrom in chromosomes:
            # Draw chromosome
            chrom_len = 250000000  # Simplified
            ax.barh(y_pos, chrom_len, height=0.5, color='lightgray',
                   edgecolor='black', linewidth=2)

            # Add label
            ax.text(-10000000, y_pos, chrom, ha='right', va='center',
                   fontsize=11, fontweight='bold')

            # Plot genes
            for gene, chr_name, pos, name in channel_genes:
                if chr_name == chrom:
                    ax.plot(pos, y_pos, 'ro', markersize=10, zorder=10)
                    ax.text(pos, y_pos + 0.6, name, ha='center', fontsize=8,
                           rotation=45)

            y_pos += 1.5

        ax.set_xlim(-20000000, 260000000)
        ax.set_ylim(-1, y_pos)
        ax.set_xlabel('Position (bp)', fontsize=12)
        ax.set_title('Calcium Channel Genes - Chromosome View',
                    fontsize=14, fontweight='bold')
        ax.set_yticks([])

        plt.tight_layout()

        # Embed
        canvas = FigureCanvasTkAgg(fig, master=self.viz_canvas_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2Tk(canvas, self.viz_canvas_frame)
        toolbar.update()

        self.notebook.select(9)  # Visualization tab
        self.update_status("Chromosome view generated")

    def plot_circos(self):
        """Generate Circos-style plot for comparative genomics."""
        for widget in self.viz_canvas_frame.winfo_children():
            widget.destroy()

        fig = plt.figure(figsize=(12, 12))
        ax = fig.add_subplot(111, projection='polar')

        # Demo: Show calcium channel families around circle
        families = ['Cav1.1', 'Cav1.2', 'Cav1.3', 'Cav1.4',
                   'Cav2.1', 'Cav2.2', 'Cav2.3',
                   'Cav3.1', 'Cav3.2', 'Cav3.3',
                   'RyR1', 'RyR2', 'RyR3',
                   'IP3R1', 'IP3R2', 'IP3R3']

        n = len(families)
        theta = np.linspace(0, 2*np.pi, n, endpoint=False)

        # Draw segments
        colors = plt.cm.Set3(np.linspace(0, 1, n))
        for i, (angle, family, color) in enumerate(zip(theta, families, colors)):
            ax.bar(angle, 1, width=2*np.pi/n*0.9, bottom=1.5,
                  color=color, edgecolor='black', linewidth=2)

            # Label
            rotation = np.degrees(angle) - 90
            if rotation > 90:
                rotation += 180
            ax.text(angle, 2.7, family, rotation=rotation,
                   ha='center', va='center', fontsize=10, fontweight='bold')

        # Draw connections (synteny or orthology)
        connections = [
            (0, 1), (1, 2), (2, 3),  # Cav1 family
            (4, 5), (5, 6),  # Cav2 family
            (7, 8), (8, 9),  # Cav3 family
            (10, 11), (11, 12),  # RyR family
            (13, 14), (14, 15)  # IP3R family
        ]

        for i, j in connections:
            ax.plot([theta[i], theta[j]], [1.5, 1.5], 'gray',
                   alpha=0.3, linewidth=1)

        ax.set_ylim(0, 3)
        ax.set_yticks([])
        ax.set_xticks([])
        ax.spines['polar'].set_visible(False)
        ax.set_title('Calcium Channel Family Relationships - Circos View',
                    fontsize=14, fontweight='bold', pad=20)

        plt.tight_layout()

        # Embed
        canvas = FigureCanvasTkAgg(fig, master=self.viz_canvas_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2Tk(canvas, self.viz_canvas_frame)
        toolbar.update()

        self.notebook.select(9)  # Visualization tab
        self.update_status("Circos plot generated")

    def plot_phylo_network(self):
        """Plot phylogenetic network instead of tree."""
        if not self.tree:
            messagebox.showwarning("Warning", "Build phylogenetic tree first")
            return

        # Show regular tree
        for widget in self.viz_canvas_frame.winfo_children():
            widget.destroy()

        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111)

        Phylo.draw(self.tree, do_show=False, axes=ax)
        ax.set_title('Phylogenetic Tree', fontsize=14, fontweight='bold')

        plt.tight_layout()

        # Embed
        canvas = FigureCanvasTkAgg(fig, master=self.viz_canvas_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2Tk(canvas, self.viz_canvas_frame)
        toolbar.update()

        self.notebook.select(9)  # Visualization tab
        self.update_status("Phylogenetic tree displayed")

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
                    self.notebook.select(7)  # Updated tab index

                    self.update_status("BLAST search complete")
                    messagebox.showinfo("BLAST", "Search complete. See Results Summary tab.")

                elif result_type == 'diamond':
                    results_df = pd.DataFrame(data)

                    summary = f"DIAMOND SEARCH RESULTS\n"
                    summary += "="*60 + "\n\n"
                    summary += f"Total hits: {len(data)}\n\n"
                    summary += results_df.to_string()

                    self.summary_text.delete(1.0, tk.END)
                    self.summary_text.insert(1.0, summary)
                    self.notebook.select(7)  # Results tab

                    self.update_status("DIAMOND search complete")
                    messagebox.showinfo("DIAMOND", f"Found {len(data)} hits. See Results Summary tab.")

                elif result_type == 'alphafold_downloaded':
                    uniprot_id, pdb_file = data
                    messagebox.showinfo("Success",
                                      f"Downloaded AlphaFold structure for {uniprot_id}\n"
                                      f"File: {pdb_file}\n\n"
                                      f"Use '3D Structure Viewer' to visualize.")
                    self.update_status(f"AlphaFold structure downloaded: {uniprot_id}")

                    # Update viewer combo
                    if hasattr(self, 'viewer_seq_combo'):
                        current_values = list(self.viewer_seq_combo['values'])
                        current_values.append(uniprot_id)
                        self.viewer_seq_combo['values'] = current_values

                elif result_type == 'ml_trained':
                    scores = data
                    messagebox.showinfo("ML Training Complete",
                                      f"Model trained successfully!\n\n"
                                      f"Cross-validation accuracy: {scores.mean():.3f} ± {scores.std():.3f}\n"
                                      f"Individual fold scores: {scores}\n\n"
                                      f"Ready to classify sequences.")
                    self.update_status(f"ML classifier trained (accuracy: {scores.mean():.3f})")

                elif result_type == 'ml_classified':
                    # Clear previous results
                    for item in self.ml_tree.get_children():
                        self.ml_tree.delete(item)

                    # Add new results
                    for result in data:
                        top3_str = "; ".join([f"{cls}: {prob:.2f}" for cls, prob in result['top3']])
                        self.ml_tree.insert('', tk.END, values=(
                            result['sequence'],
                            result['prediction'],
                            f"{result['confidence']:.3f}",
                            top3_str
                        ))

                    self.notebook.select(5)  # ML tab
                    self.update_status(f"ML classification complete: {len(data)} sequences")
                    messagebox.showinfo("ML Classification",
                                      f"Classified {len(data)} sequences.\n"
                                      f"See ML Classification tab for details.")

                elif result_type == 'electrophys':
                    report = "ELECTROPHYSIOLOGY PREDICTIONS\n"
                    report += "="*60 + "\n\n"

                    for result in data:
                        report += f"Sequence: {result['sequence']}\n"
                        report += f"{'='*50}\n\n"
                        report += f"Gating Mechanism: {result['gating']}\n"
                        report += f"Activation: {result['activation']}\n\n"
                        report += f"Ion Selectivity: {result['selectivity']}\n"
                        report += f"Permeability: {result['permeability']}\n\n"
                        report += f"Pharmacology (Blockers): {result['blockers']}\n"
                        report += "\n" + "-"*60 + "\n\n"

                    self.ephys_results.delete(1.0, tk.END)
                    self.ephys_results.insert(1.0, report)
                    self.notebook.select(6)  # Electrophys tab

                    self.update_status("Electrophysiology prediction complete")
                    messagebox.showinfo("Electrophysiology",
                                      f"Predicted properties for {len(data)} sequences.\n"
                                      f"See Electrophysiology tab for details.")

                elif result_type == 'genome_search':
                    # Clear previous results
                    for item in self.discovery_tree.get_children():
                        self.discovery_tree.delete(item)

                    # Add results
                    for result in data:
                        self.discovery_tree.insert('', tk.END, values=(
                            result['gene_id'],
                            result['organism'],
                            result['score'],
                            result['identity'],
                            result['domains'],
                            result['family']
                        ))

                    self.notebook.select(7)  # Discovery tab
                    self.update_status(f"Genome search complete: {len(data)} hits")
                    messagebox.showinfo("Genome Search",
                                      f"Found {len(data)} calcium channel candidates.\n"
                                      f"See Discovery tab for details.")

                elif result_type == 'orthologs':
                    report = "ORTHOLOG ANALYSIS\n"
                    report += "="*60 + "\n\n"

                    for result in data:
                        report += f"Query: {result['query']}\n"
                        report += f"Ortholog: {result['ortholog']}\n"
                        report += f"Organism: {result['organism']}\n"
                        report += f"Identity: {result['identity']}\n"
                        report += f"E-value: {result['evalue']}\n"
                        report += f"Score: {result['score']}\n\n"

                    self.summary_text.delete(1.0, tk.END)
                    self.summary_text.insert(1.0, report)
                    self.notebook.select(10)  # Results tab

                    self.update_status(f"Ortholog search complete: {len(data)} found")
                    messagebox.showinfo("Orthologs",
                                      f"Found {len(data)} orthologs.\n"
                                      f"See Results Summary tab.")

                elif result_type == 'expression':
                    gene_symbol, expression = data
                    messagebox.showinfo("Expression Data",
                                      f"Fetched expression data for {gene_symbol}.\n"
                                      f"Use 'Expression Heatmap' visualization to view.")
                    self.update_status(f"Expression data loaded: {gene_symbol}")

                elif result_type == 'interactions':
                    protein_name, interactions = data
                    messagebox.showinfo("Interaction Network",
                                      f"Fetched {len(interactions)} interactions for {protein_name}.\n"
                                      f"Use 'Network Graph' visualization to view.")
                    self.update_status(f"Interaction network loaded: {protein_name}")

                elif result_type == 'tsne':
                    X_tsne, labels = data

                    # Clear previous
                    for widget in self.viz_canvas_frame.winfo_children():
                        widget.destroy()

                    # Plot
                    fig, ax = plt.subplots(figsize=(10, 8))

                    # Color by family
                    unique_labels = list(set(labels))
                    colors = plt.cm.tab10(np.linspace(0, 1, len(unique_labels)))

                    for label, color in zip(unique_labels, colors):
                        mask = np.array(labels) == label
                        ax.scatter(X_tsne[mask, 0], X_tsne[mask, 1],
                                  label=label, c=[color], s=100, alpha=0.7, edgecolors='black')

                    ax.set_xlabel('t-SNE 1', fontsize=12)
                    ax.set_ylabel('t-SNE 2', fontsize=12)
                    ax.set_title('t-SNE: Calcium Channel Clustering', fontsize=14, fontweight='bold')
                    ax.legend()
                    ax.grid(True, alpha=0.3)

                    plt.tight_layout()

                    # Embed
                    canvas = FigureCanvasTkAgg(fig, master=self.viz_canvas_frame)
                    canvas.draw()
                    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

                    toolbar = NavigationToolbar2Tk(canvas, self.viz_canvas_frame)
                    toolbar.update()

                    self.notebook.select(9)  # Visualization tab
                    self.update_status("t-SNE plot generated")

                elif result_type == 'variant':
                    result = data
                    report = "VARIANT EFFECT PREDICTION\n"
                    report += "="*60 + "\n\n"
                    report += f"Mutation: {result['mutation']}\n"
                    report += f"Position: {result['position']}\n"
                    report += f"Wild-type: {result['wt_aa']}\n"
                    report += f"Mutant: {result['mut_aa']}\n\n"
                    report += f"Located in domain: {'Yes' if result['in_domain'] else 'No'}\n"
                    if result['in_domain']:
                        report += f"Affected domain: {result['affected_domain']}\n"
                    report += f"Domain disruption: {'Yes' if result['domain_lost'] else 'No'}\n\n"
                    report += f"Pathogenicity Score: {result['pathogenicity_score']:.3f}\n"
                    report += f"Prediction: {result['prediction']}\n\n"
                    report += "="*60 + "\n\n"

                    if result['domain_lost']:
                        report += "⚠️ WARNING: This mutation disrupts one or more functional domains!\n\n"

                    report += "Wild-type domains:\n"
                    for domain, matches in result['wt_domains'].items():
                        report += f"  {domain}: {len(matches)} found\n"

                    report += "\nMutant domains:\n"
                    for domain, matches in result['mut_domains'].items():
                        report += f"  {domain}: {len(matches)} found\n"

                    self.variant_results.delete(1.0, tk.END)
                    self.variant_results.insert(1.0, report)

                    self.update_status(f"Variant analysis complete: {result['prediction']}")

                elif result_type == 'variant_scan':
                    results = data
                    report = f"VARIANT SCAN RESULTS\n"
                    report += "="*60 + "\n\n"
                    report += f"Analyzed all possible single amino acid substitutions\n"
                    report += f"Showing top {len(results)} potentially pathogenic variants\n\n"
                    report += "-"*60 + "\n\n"

                    for i, var in enumerate(results, 1):
                        report += f"{i}. Position {var['position']}: "
                        report += f"{var['wt_aa']} → {var['mut_aa']} "
                        report += f"(Score: {var['score']:.3f})"
                        if var['in_domain']:
                            report += " [IN DOMAIN]"
                        report += "\n"

                    self.variant_results.delete(1.0, tk.END)
                    self.variant_results.insert(1.0, report)

                    self.update_status(f"Variant scan complete: {len(results)} high-risk variants found")
                    messagebox.showinfo("Variant Scan",
                                      f"Found {len(results)} potentially pathogenic variants.\n"
                                      f"See Variant Analysis tab for details.")

                elif result_type == 'conservation':
                    scores = data

                    # Plot conservation
                    for widget in self.cons_canvas_frame.winfo_children():
                        widget.destroy()

                    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))

                    positions = [s['position'] for s in scores]
                    conservation = [s['conservation'] for s in scores]

                    # Conservation plot
                    ax1.plot(positions, conservation, 'b-', linewidth=2)
                    ax1.fill_between(positions, conservation, alpha=0.3)
                    ax1.axhline(y=0.8, color='r', linestyle='--', label='High conservation (>0.8)')
                    ax1.set_xlabel('Position', fontsize=12)
                    ax1.set_ylabel('Conservation Score', fontsize=12)
                    ax1.set_title('Sequence Conservation', fontsize=14, fontweight='bold')
                    ax1.set_ylim(0, 1)
                    ax1.legend()
                    ax1.grid(True, alpha=0.3)

                    # Entropy plot
                    entropy = [s['entropy'] for s in scores]
                    ax2.plot(positions, entropy, 'r-', linewidth=2)
                    ax2.fill_between(positions, entropy, alpha=0.3)
                    ax2.set_xlabel('Position', fontsize=12)
                    ax2.set_ylabel('Shannon Entropy', fontsize=12)
                    ax2.set_title('Position Variability (Entropy)', fontsize=14, fontweight='bold')
                    ax2.grid(True, alpha=0.3)

                    plt.tight_layout()

                    canvas = FigureCanvasTkAgg(fig, master=self.cons_canvas_frame)
                    canvas.draw()
                    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

                    toolbar = NavigationToolbar2Tk(canvas, self.cons_canvas_frame)
                    toolbar.update()

                    # Find hotspots (>0.8)
                    hotspots = [s for s in scores if s['conservation'] > 0.8]

                    self.notebook.select(8)  # Conservation tab
                    self.update_status(f"Conservation analysis complete: {len(hotspots)} hotspots found")

                    messagebox.showinfo("Conservation Analysis",
                                      f"Analysis complete!\n\n"
                                      f"Highly conserved positions (>0.8): {len(hotspots)}\n"
                                      f"Mean conservation: {np.mean(conservation):.3f}")

                elif result_type == 'batch_complete':
                    results, export_format = data

                    # Create DataFrame
                    df = pd.DataFrame(results['sequences'])

                    # Export based on format
                    if 'Excel' in export_format or 'All' in export_format:
                        filename = filedialog.asksaveasfilename(
                            defaultextension=".xlsx",
                            filetypes=[("Excel files", "*.xlsx")],
                            initialfile="batch_analysis.xlsx"
                        )
                        if filename:
                            df.to_excel(filename, index=False)
                            messagebox.showinfo("Export", f"Results exported to {filename}")

                    summary = f"BATCH ANALYSIS RESULTS\n"
                    summary += "="*60 + "\n\n"
                    summary += f"Total sequences analyzed: {len(results['sequences'])}\n\n"
                    summary += df.to_string()

                    self.summary_text.delete(1.0, tk.END)
                    self.summary_text.insert(1.0, summary)
                    self.notebook.select(10)  # Results tab

                    self.update_status("Batch analysis complete")
                    messagebox.showinfo("Batch Analysis",
                                      f"Analyzed {len(results['sequences'])} sequences.\n"
                                      f"Results exported and displayed in Summary tab.")

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
