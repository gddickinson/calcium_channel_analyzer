"""
Tests for ML feature extraction.
"""

import sys
import os
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from ml_features import extract_ml_features
from config import DOMAIN_PATTERNS, AMINO_ACIDS


class TestMLFeatures:
    """Test ML feature extraction from sequences."""

    def test_feature_count(self):
        """Features = 1 (length) + 20 (AA comp) + 1 (hydro) + N (domains)."""
        seq = "ACDEFGHIKLMNPQRSTVWY"
        features = extract_ml_features(seq)
        expected = 1 + 20 + 1 + len(DOMAIN_PATTERNS)
        assert len(features) == expected

    def test_length_feature(self):
        seq = "ACDEFGHIKLMNPQRSTVWY"
        features = extract_ml_features(seq)
        assert features[0] == float(len(seq))

    def test_composition_sums_to_one(self):
        """AA composition features (indices 1-20) should sum to ~1.0."""
        seq = "ACDEFGHIKLMNPQRSTVWY"
        features = extract_ml_features(seq)
        aa_comp = features[1:21]
        assert abs(sum(aa_comp) - 1.0) < 0.01

    def test_polyalanine_composition(self):
        """Poly-A should have A=1.0 and all others 0."""
        seq = "A" * 50
        features = extract_ml_features(seq)
        # A is the first amino acid in AMINO_ACIDS
        assert features[1] == 1.0
        # All other AAs should be 0
        for i in range(2, 21):
            assert features[i] == 0.0

    def test_hydrophobicity_positive_for_hydrophobic(self):
        """Poly-I (isoleucine, very hydrophobic) should have positive hydro."""
        seq = "I" * 50
        features = extract_ml_features(seq)
        hydro_idx = 21
        assert features[hydro_idx] > 0

    def test_hydrophobicity_negative_for_charged(self):
        """Poly-D (aspartate, very hydrophilic) should have negative hydro."""
        seq = "D" * 50
        features = extract_ml_features(seq)
        hydro_idx = 21
        assert features[hydro_idx] < 0

    def test_empty_sequence(self):
        """Empty sequence should return all zeros except length."""
        features = extract_ml_features("")
        assert features[0] == 0.0
        assert all(f == 0.0 for f in features[1:])
