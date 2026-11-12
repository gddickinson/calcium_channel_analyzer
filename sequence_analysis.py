"""
Calcium Channel Analyzer -- Sequence analysis (GUI-free)
=========================================================

Pure functions behind the desktop app's sequence analyses, usable without
Tkinter (the web app and tests import them):

  - parse_fasta()          FASTA text -> records (id = first header word, like Bio.SeqIO)
  - hydropathy_profile()   Kyte-Doolittle sliding-window hydropathy
  - predict_tm_regions()   TM regions where the profile exceeds a threshold
  - tm_prediction()        both of the above for one sequence

The hydropathy/TM logic is the same as CalciumChannelAnalyzer.predict_tm_domains
(window 19, threshold 1.6, region = window centres while the score is above
the threshold).
"""

import re
from typing import Dict, List, Tuple

from config import KD_HYDROPHOBICITY

TM_WINDOW = 19
TM_THRESHOLD = 1.6
MAX_SEQUENCES = 200
MAX_RESIDUES = 500_000


class FastaError(ValueError):
    """Raised for FASTA input that cannot be used at all."""


def parse_fasta(text: str) -> Tuple[List[Dict], List[str]]:
    """
    Parse FASTA text into records.

    Each record is {'id', 'description', 'sequence'}: the id is the first word
    of the header and the description the whole header line (as Bio.SeqIO
    gives them). Sequences are upper-cased; whitespace, digits and a trailing
    '*' are removed; other non-letter characters are dropped with a warning.
    Records with no residues are skipped with a warning.

    Returns:
        (records, warnings)

    Raises:
        FastaError: no header line, no usable record, or input too large.
    """
    records: List[Dict] = []
    warnings: List[str] = []
    header = None
    chunks: List[str] = []
    lines = (text or "").replace("\r", "").split("\n")
    if not any(line.strip() for line in lines):
        raise FastaError("No sequences: paste FASTA text starting with a '>' header line.")

    def finish():
        if header is None:
            return
        raw = re.sub(r"[\s\d]", "", "".join(chunks)).upper().rstrip("*")
        bad = sorted(set(re.sub(r"[A-Z]", "", raw)))
        seq = re.sub(r"[^A-Z]", "", raw)
        ident = header.split()[0] if header.split() else f"seq{len(records) + 1}"
        if bad:
            warnings.append(f"{ident}: removed characters {' '.join(bad)}")
        if not seq:
            warnings.append(f"{ident}: empty sequence skipped")
            return
        records.append({"id": ident, "description": header, "sequence": seq})

    for n, line in enumerate(lines, 1):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            finish()
            header = line[1:].strip()
            chunks = []
        elif header is None:
            raise FastaError(f"Line {n}: sequence data before the first '>' header line.")
        else:
            chunks.append(line)
    finish()

    if not records:
        raise FastaError("No usable sequences found.")
    if len(records) > MAX_SEQUENCES:
        raise FastaError(f"Too many sequences ({len(records)}); the limit is {MAX_SEQUENCES}.")
    total = sum(len(r["sequence"]) for r in records)
    if total > MAX_RESIDUES:
        raise FastaError(f"Too many residues ({total}); the limit is {MAX_RESIDUES}.")
    seen: Dict[str, int] = {}
    for r in records:
        seen[r["id"]] = seen.get(r["id"], 0) + 1
    dups = [k for k, v in seen.items() if v > 1]
    if dups:
        warnings.append("Duplicate ids: " + ", ".join(dups))
    return records, warnings


def hydropathy_profile(sequence: str, window: int = TM_WINDOW) -> Tuple[List[int], List[float]]:
    """
    Kyte-Doolittle hydropathy averaged over a sliding window.

    Returns:
        (positions, scores): positions are 0-based window centres
        (i + window // 2); unknown residues score 0.
    """
    positions: List[int] = []
    scores: List[float] = []
    for i in range(len(sequence) - window + 1):
        window_seq = sequence[i:i + window]
        scores.append(sum(KD_HYDROPHOBICITY.get(aa, 0) for aa in window_seq) / window)
        positions.append(i + window // 2)
    return positions, scores


def predict_tm_regions(
    positions: List[int], scores: List[float], threshold: float = TM_THRESHOLD
) -> List[Tuple[int, int]]:
    """
    TM regions from a hydropathy profile: a region opens at the first window
    centre scoring above the threshold and closes at the first centre at or
    below it (or at the last centre). Returns (start, end) pairs, 0-based;
    the region's residues are sequence[start:end].
    """
    regions: List[Tuple[int, int]] = []
    in_tm = False
    start = 0
    for i, score in enumerate(scores):
        if score > threshold and not in_tm:
            start = positions[i]
            in_tm = True
        elif score <= threshold and in_tm:
            regions.append((start, positions[i]))
            in_tm = False
    if in_tm:
        regions.append((start, positions[-1]))
    return regions


def tm_prediction(sequence: str, window: int = TM_WINDOW, threshold: float = TM_THRESHOLD) -> Dict:
    """Hydropathy profile and predicted TM regions for one sequence."""
    positions, scores = hydropathy_profile(sequence, window)
    regions = predict_tm_regions(positions, scores, threshold)
    return {
        "length": len(sequence),
        "window": window,
        "threshold": threshold,
        "positions": positions,
        "scores": scores,
        "regions": [
            {"start": s, "end": e, "length": e - s, "sequence": sequence[s:e]} for s, e in regions
        ],
    }
