"""
Tests for domain detection on known calcium channel sequences.
"""

import sys
import os
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from domain_detection import (
    detect_domains_in_sequence,
    calculate_domain_score,
    predict_family_from_description,
)
from config import DOMAIN_PATTERNS


class TestDomainDetection:
    """Test domain detection against known sequences."""

    def test_ef_hand_detected(self):
        """EF-hand motif: D.{1,3}D.{1,3}[DNS].{1,3}[EDNST] with min_length 12."""
        # Use 3-char gaps to ensure match length >= 12
        # D + 3 + D + 3 + D + 3 + E = 13 chars
        seq = "AAADAAADAAADAAAEAAAA"
        results = detect_domains_in_sequence(seq, seq_id="test_ef")
        domain_names = [r['domain'] for r in results]
        assert 'EF_hand' in domain_names

    def test_pore_loop_detected(self):
        """Ion channel selectivity filter: T[LI]LG[LI]G."""
        seq = "AAATLLGLGAAAA"
        results = detect_domains_in_sequence(seq, seq_id="test_pore")
        domain_names = [r['domain'] for r in results]
        assert 'Pore_loop' in domain_names

    def test_voltage_sensor_detected(self):
        """Voltage-sensing domain S4: R/K repeats."""
        seq = "AAARKRKVVVVVRKRKAAAA"
        results = detect_domains_in_sequence(seq, seq_id="test_vs")
        domain_names = [r['domain'] for r in results]
        assert 'Voltage_sensor' in domain_names

    def test_empty_sequence_returns_empty(self):
        results = detect_domains_in_sequence("", seq_id="empty")
        assert results == []

    def test_no_domains_in_polyalanine(self):
        """Poly-alanine should not match most domains."""
        seq = "A" * 100
        results = detect_domains_in_sequence(seq, seq_id="polyA")
        # Should only potentially match TM_helix (hydrophobic), filter rest
        non_tm = [r for r in results if r['domain'] != 'TM_helix']
        # Most domains require charged/specific residues
        assert len(non_tm) == 0


class TestDomainScore:
    """Test domain scoring."""

    def test_score_is_positive(self):
        score = calculate_domain_score("DKDLDNLE", "EF_hand")
        assert score > 0

    def test_tm_helix_score_uses_hydrophobicity(self):
        # All hydrophobic residues
        score_hydro = calculate_domain_score("AILVMFGWAILVMFGWAILVMFGW", "TM_helix")
        # Mixed residues
        score_mixed = calculate_domain_score("AILVDDDKAILVDDDKAILVDDDK", "TM_helix")
        assert score_hydro > score_mixed


class TestFamilyPrediction:
    """Test channel family prediction from descriptions."""

    def test_ryr_detected(self):
        assert predict_family_from_description("Ryanodine receptor 1") == "RyR"

    def test_ltype_detected(self):
        result = predict_family_from_description("CACNA1C L-type calcium channel")
        assert "L-type" in result

    def test_unknown_for_garbage(self):
        assert predict_family_from_description("random protein xyz") == "Unknown"

    def test_cngc_detected(self):
        result = predict_family_from_description("Cyclic nucleotide gated channel 4")
        assert result == "CNGC"

    def test_tpc_detected(self):
        result = predict_family_from_description("Two-pore channel 1")
        assert result == "TPC"
