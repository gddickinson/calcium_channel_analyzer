"""
Calcium Channel Analyzer -- Domain Detection
==============================================

Functions for detecting conserved domains in calcium channel sequences.
Uses regex-based pattern matching against known domain signatures.
"""

import re
from typing import List, Dict, Optional

from config import DOMAIN_PATTERNS


def detect_domains_in_sequence(
    sequence: str,
    seq_id: str = "unknown",
    domain_patterns: Optional[Dict] = None,
    min_length: Optional[int] = None,
    use_pattern_min_length: bool = False,
) -> List[Dict]:
    """
    Detect conserved domains in a protein sequence.

    Every non-overlapping regex match is kept (as the desktop app's Detect
    Domains does) unless a minimum length is asked for.

    Args:
        sequence: Amino acid sequence string.
        seq_id: Identifier for the sequence.
        domain_patterns: Domain pattern dict (defaults to DOMAIN_PATTERNS).
        min_length: If given, drop matches shorter than this many residues.
        use_pattern_min_length: If True, drop matches shorter than each
            pattern's own 'min_length'.

    Returns:
        List of dicts with keys: sequence, domain, position, start, end,
        match, score. Positions are 1-based and inclusive: 'position' is
        "start-end" and sequence[start - 1:end] == match.
    """
    if domain_patterns is None:
        domain_patterns = DOMAIN_PATTERNS

    results = []
    for domain_name, domain_info in domain_patterns.items():
        pattern = domain_info['pattern']
        for match in re.finditer(pattern, sequence):
            matched_seq = match.group()
            if min_length is not None and len(matched_seq) < min_length:
                continue
            if use_pattern_min_length and len(matched_seq) < domain_info.get('min_length', 1):
                continue
            score = calculate_domain_score(matched_seq, domain_name, domain_patterns)
            start, end = match.start() + 1, match.end()
            results.append({
                'sequence': seq_id,
                'domain': domain_name,
                'position': f"{start}-{end}",
                'start': start,
                'end': end,
                'match': matched_seq,
                'score': f"{score:.1f}",
            })
    return results


def calculate_domain_score(
    sequence: str,
    domain_name: str,
    domain_patterns: Optional[Dict] = None,
) -> float:
    """
    Calculate a simple conservation score for a domain match.

    Args:
        sequence: The matched subsequence.
        domain_name: Name of the domain.
        domain_patterns: Domain pattern dict.

    Returns:
        Numeric score (0-100).
    """
    if domain_patterns is None:
        domain_patterns = DOMAIN_PATTERNS

    domain_info = domain_patterns[domain_name]
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


def predict_family_from_description(description: str) -> str:
    """
    Predict channel family from sequence description string.

    Args:
        description: Sequence description/header text.

    Returns:
        Predicted family name string.
    """
    description_lower = description.lower()

    family_keywords = {
        'Cav1 (L-type)': ['cav1', 'l-type', 'cacna1s', 'cacna1c',
                           'cacna1d', 'cacna1f', 'dihydropyridine'],
        'Cav2 (P/Q/N/R-type)': ['cav2', 'p/q-type', 'n-type', 'r-type',
                                  'cacna1a', 'cacna1b', 'cacna1e'],
        'Cav3 (T-type)': ['cav3', 't-type', 'cacna1g', 'cacna1h', 'cacna1i'],
        'RyR': ['ryanodine', 'ryr'],
        'IP3R': ['ip3', 'inositol', 'itpr'],
        'CNGC': ['cyclic nucleotide', 'cngc', 'cng channel'],
        'GLR': ['glutamate receptor-like', 'glr', 'iglur'],
        'TPC': ['two-pore', 'tpc', 'tpcn'],
        'TRP': ['transient receptor', 'trp', 'trpc', 'trpv', 'trpm'],
        'Orai': ['orai', 'crac', 'stim'],
        'P2X': ['p2x', 'purinergic'],
        'MSL': ['mechanosensitive', 'msl', 'mscs'],
        'MCA': ['mid1-complementing', 'mca'],
    }

    for family, keywords in family_keywords.items():
        for keyword in keywords:
            if keyword in description_lower:
                return family

    return 'Unknown'
