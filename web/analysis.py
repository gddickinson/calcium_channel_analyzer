"""Analysis bundles for the web API: FASTA text in, JSON-ready dicts out.

Every function calls the repo's own modules (config, domain_detection, ml_features, sequence_analysis, channel_rules,
alignment_tree); nothing here re-implements an analysis. The same functions generate the JS-engine fixture
(scripts/make_engine_fixture.py), so the browser engines are checked against exactly what the server returns.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from alignment_tree import (ALIGN_PREFIX, GAP_EXTEND, GAP_OPEN, MATCH, MISMATCH,  # noqa: E402
                            build_tree, identity_distances, pairwise_alignments, tree_to_dict, tree_to_newick)
from channel_rules import classify_channel, predict_electrophysiology  # noqa: E402
from config import AMINO_ACIDS, CHANNEL_FAMILIES, DOMAIN_PATTERNS, KD_HYDROPHOBICITY  # noqa: E402
from domain_detection import detect_domains_in_sequence, predict_family_from_description  # noqa: E402
from ml_features import extract_ml_features  # noqa: E402
from sequence_analysis import TM_THRESHOLD, TM_WINDOW, parse_fasta, tm_prediction  # noqa: E402

EXAMPLE_FASTA = ROOT / "example_ca_channels.fasta"
MAX_ALIGN_SEQS = 12
FEATURE_NAMES = (["length"] + [f"frac_{aa}" for aa in AMINO_ACIDS] + ["mean_kd_hydropathy"]
                 + [f"has_{d}" for d in DOMAIN_PATTERNS])


def config_payload() -> dict:
    """The constants the browser engines need (single source of truth: config.py)."""
    return {
        "domain_patterns": DOMAIN_PATTERNS,
        "channel_families": CHANNEL_FAMILIES,
        "kd_hydrophobicity": KD_HYDROPHOBICITY,
        "amino_acids": AMINO_ACIDS,
        "feature_names": FEATURE_NAMES,
        "tm": {"window": TM_WINDOW, "threshold": TM_THRESHOLD},
        "align": {"prefix": ALIGN_PREFIX, "match": MATCH, "mismatch": MISMATCH,
                  "gap_open": GAP_OPEN, "gap_extend": GAP_EXTEND, "max_sequences": MAX_ALIGN_SEQS},
    }


def example_fasta() -> str:
    return EXAMPLE_FASTA.read_text(encoding="utf-8")


def parse(fasta: str) -> dict:
    records, warnings = parse_fasta(fasta)
    for r in records:
        r["length"] = len(r["sequence"])
        r["family"] = predict_family_from_description(r["description"])
    return {"records": records, "warnings": warnings}


def domains(fasta: str, use_min_length: bool = False) -> dict:
    """domain_detection.detect_domains_in_sequence on every record (positions 1-based, inclusive). Every match is kept
    unless use_min_length, which applies each pattern's min_length."""
    records, warnings = parse_fasta(fasta)
    hits = []
    for r in records:
        hits += detect_domains_in_sequence(r["sequence"], r["id"], use_pattern_min_length=use_min_length)
    return {"use_min_length": use_min_length, "lengths": {r["id"]: len(r["sequence"]) for r in records}, "order": [r["id"] for r in records],
            "hits": hits, "warnings": warnings}


def hydropathy(fasta: str, index: int = 0, window: int = TM_WINDOW, threshold: float = TM_THRESHOLD) -> dict:
    records, warnings = parse_fasta(fasta)
    if not 0 <= index < len(records):
        raise ValueError(f"index must be between 0 and {len(records) - 1}")
    if not 3 <= window <= 51:
        raise ValueError("window must be between 3 and 51")
    r = records[index]
    return {"id": r["id"], "description": r["description"], **tm_prediction(r["sequence"], window, threshold),
            "warnings": warnings}


def classify(fasta: str) -> dict:
    records, warnings = parse_fasta(fasta)
    rows = []
    for r in records:
        rows.append({
            "id": r["id"], "description": r["description"],
            "family_from_description": predict_family_from_description(r["description"]),
            **classify_channel(r["description"], r["sequence"]),
            "electrophysiology": predict_electrophysiology(r["description"], r["sequence"]),
            "features": extract_ml_features(r["sequence"]),
        })
    return {"rows": rows, "feature_names": FEATURE_NAMES, "warnings": warnings}


def align(fasta: str) -> dict:
    records, warnings = parse_fasta(fasta)
    if len(records) < 2:
        raise ValueError("Need at least 2 sequences for alignment")
    if len(records) > MAX_ALIGN_SEQS:
        raise ValueError(f"At most {MAX_ALIGN_SEQS} sequences for all-pairs alignment")
    ids = [r["id"] for r in records]
    return {"ids": ids, "prefix": ALIGN_PREFIX,
            "pairs": pairwise_alignments(ids, [r["sequence"] for r in records]), "warnings": warnings}


def tree(fasta: str, method: str = "UPGMA") -> dict:
    records, warnings = parse_fasta(fasta)
    method = method.upper()
    if method not in ("UPGMA", "NJ"):
        raise ValueError("method must be UPGMA or NJ")
    ids, seqs = [r["id"] for r in records], [r["sequence"] for r in records]
    t = build_tree(ids, seqs, method)
    return {"ids": ids, "method": method, "distances": identity_distances(ids, seqs),
            "tree": tree_to_dict(t.root), "newick": tree_to_newick(t), "warnings": warnings}
