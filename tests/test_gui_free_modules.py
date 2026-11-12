"""Tests for the GUI-free modules factored out for the web app: sequence_analysis, channel_rules, alignment_tree.

They copy logic that is still inline in calcium_channel_analyzer.py; test_desktop_constants_unchanged guards against
the two drifting apart (it reads the desktop source as text, so Tkinter is not needed).
"""
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from alignment_tree import build_tree, identity_distances, make_aligner, pairwise_alignments  # noqa: E402
from channel_rules import classify_channel, predict_electrophysiology  # noqa: E402
from sequence_analysis import FastaError, hydropathy_profile, parse_fasta, predict_tm_regions, tm_prediction  # noqa: E402


def test_parse_fasta_ids_cleanup_and_warnings():
    recs, warns = parse_fasta(">sp|X|Y desc here\nmk v1*\n>empty\n>z\nAC-D\n")
    assert [r["id"] for r in recs] == ["sp|X|Y", "z"]
    assert recs[0]["sequence"] == "MKV" and recs[0]["description"] == "sp|X|Y desc here"
    assert recs[1]["sequence"] == "ACD"
    assert any("empty sequence" in w for w in warns) and any("removed characters -" in w for w in warns)


@pytest.mark.parametrize("text", ["", "MKV\n>a\nMKV", ">a\n\n"])
def test_parse_fasta_rejects(text):
    with pytest.raises(FastaError):
        parse_fasta(text)


def test_hydropathy_and_tm_regions():
    seq = "D" * 20 + "I" * 30 + "D" * 20
    pos, sc = hydropathy_profile(seq, 19)
    assert len(sc) == len(seq) - 18 and pos[0] == 9
    assert sc[0] == pytest.approx(-3.5) and max(sc) == pytest.approx(4.5)
    regions = predict_tm_regions(pos, sc, 1.6)
    assert len(regions) == 1
    s, e = regions[0]
    assert 15 < s < 25 and 45 < e < 55
    assert tm_prediction("ACD")["regions"] == []  # shorter than the window


def test_tm_region_open_at_end():
    pos, sc = hydropathy_profile("D" * 20 + "I" * 30, 19)
    assert predict_tm_regions(pos, sc, 1.6)[-1][1] == pos[-1]


def test_classification_rules_order():
    assert classify_channel("CACNA1C L-type", "AAAA")["classification"] == "Voltage-gated (Cav1/L-type)"
    assert classify_channel("ryanodine receptor", "AAAA")["confidence"] == "High"
    r = classify_channel("unknown", "KRRAAAAKRR" + "DAAEAAD")
    assert r["classification"] == "Voltage-gated Ca2+ channel" and r["confidence"] == "Medium"
    assert classify_channel("nothing", "AAAA")["classification"] == "Unknown"


def test_electrophysiology_rules():
    e = predict_electrophysiology("two-pore channel", "AAAA")
    assert e["selectivity"].startswith("Non-selective") and e["gating"] == "Unknown/Mixed"
    assert predict_electrophysiology("x", "KRRAAAAKRR")["gating"] == "Voltage-gated"


def test_alignment_and_tree():
    assert make_aligner().score("ACDEFGHIK", "ACDGHIKKK") == pytest.approx(5.8)
    pairs = pairwise_alignments(["a", "b"], ["ACDEFGHIK", "ACDGHIKKK"])
    assert pairs[0]["aligned1"].replace("-", "") == "ACDEFGHIK"
    d = identity_distances(["a", "b"], ["AAAA", "AA"])
    assert d[0][1] == pytest.approx(0.5)  # 'A'/'A' x2 match, 'A'/'-' x2 differ
    assert d[0][0] == 0
    with pytest.raises(ValueError):
        build_tree(["a", "b"], ["AA", "AC"])


def test_desktop_constants_unchanged():
    src = (ROOT / "calcium_channel_analyzer.py").read_text()
    tm = src[src.index("def predict_tm_domains"):src.index("def plot_hydropathy")]
    assert "window = 19" in tm and "tm_threshold = 1.6" in tm
    al = src[src.index("def run_alignment"):src.index("def build_tree")]
    for s in ("match_score = 1", "mismatch_score = 0", "open_gap_score = -0.5", "extend_gap_score = -0.1", "[:500]"):
        assert s in al
    tr = src[src.index("def build_tree"):src.index("def plot_tree")]
    assert "DistanceCalculator('identity')" in tr and re.search(r"'-' \* \(max_len - len\(s\)\)", tr)
