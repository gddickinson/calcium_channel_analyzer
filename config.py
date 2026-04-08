"""
Calcium Channel Analyzer -- Configuration and Constants
========================================================

Domain patterns, channel family definitions, and shared constants
extracted from the monolithic analyzer for reuse and testability.
"""

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

# Kyte-Doolittle hydrophobicity scale
KD_HYDROPHOBICITY = {
    'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8, 'G': -0.4,
    'H': -3.2, 'I': 4.5, 'K': -3.9, 'L': 3.8, 'M': 1.9, 'N': -3.5,
    'P': -1.6, 'Q': -3.5, 'R': -4.5, 'S': -0.8, 'T': -0.7, 'V': 4.2,
    'W': -0.9, 'Y': -1.3
}

# Standard amino acid list
AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'
