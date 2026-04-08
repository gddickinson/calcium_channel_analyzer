"""
Calcium Channel Analyzer -- ML Feature Extraction
===================================================

Functions for extracting machine learning features from protein sequences.
"""

import re
from typing import List, Dict, Optional

from config import DOMAIN_PATTERNS, KD_HYDROPHOBICITY, AMINO_ACIDS


def extract_ml_features(
    sequence: str,
    domain_patterns: Optional[Dict] = None,
) -> List[float]:
    """
    Extract numeric features from a protein sequence for ML classification.

    Features include:
      - Sequence length
      - Amino acid composition (20 features)
      - Average hydrophobicity (Kyte-Doolittle)
      - Binary domain presence flags (one per domain pattern)

    Args:
        sequence: Amino acid sequence string.
        domain_patterns: Domain patterns dict (defaults to DOMAIN_PATTERNS).

    Returns:
        List of float feature values.
    """
    if domain_patterns is None:
        domain_patterns = DOMAIN_PATTERNS

    features = []

    # Length
    seq_len = len(sequence)
    features.append(float(seq_len))

    # Amino acid composition
    for aa in AMINO_ACIDS:
        features.append(sequence.count(aa) / seq_len if seq_len > 0 else 0.0)

    # Average hydrophobicity
    avg_hydro = (
        sum(KD_HYDROPHOBICITY.get(aa, 0) for aa in sequence) / seq_len
        if seq_len > 0
        else 0.0
    )
    features.append(avg_hydro)

    # Domain presence (binary features)
    for domain_name, domain_info in domain_patterns.items():
        pattern = domain_info['pattern']
        has_domain = 1.0 if re.search(pattern, sequence) else 0.0
        features.append(has_domain)

    return features
