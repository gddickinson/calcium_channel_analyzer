"""
Calcium Channel Analyzer -- Rule-based channel classification (GUI-free)
=========================================================================

The keyword-and-motif rules used by the desktop app, as pure functions:

  - classify_channel()          rules of CalciumChannelAnalyzer.classify_channels
  - predict_electrophysiology() rules of CalciumChannelAnalyzer.predict_electrophysiology

Both look at the header description (lower-cased) first and the regex
markers from config.DOMAIN_PATTERNS second. They are heuristics: the strings
they return are the app's own rule outputs, not measurements.
"""

import re
from typing import Dict

from config import DOMAIN_PATTERNS


def _has(name: str, seq: str) -> bool:
    return bool(re.search(DOMAIN_PATTERNS[name]["pattern"], seq))


def channel_markers(sequence: str) -> Dict[str, bool]:
    """Presence of the five classification markers (regex search, anywhere)."""
    return {
        "voltage_sensor": _has("Voltage_sensor", sequence),
        "ca_selectivity": _has("Ca_selectivity", sequence),
        "cnbd": _has("CNBD", sequence),
        "spry": _has("SPRY", sequence),
        "ip3_binding": _has("IP3_binding", sequence),
    }


def classify_channel(description: str, sequence: str) -> Dict:
    """Classify one sequence: description keywords first, then domain markers."""
    desc = description.lower()
    m = channel_markers(sequence)
    classification, confidence = "Unknown", "Low"
    if any(term in desc for term in ["cav1", "l-type"]):
        classification, confidence = "Voltage-gated (Cav1/L-type)", "High"
    elif any(term in desc for term in ["cav2", "p/q", "n-type"]):
        classification, confidence = "Voltage-gated (Cav2/P-Q-N-type)", "High"
    elif any(term in desc for term in ["cav3", "t-type"]):
        classification, confidence = "Voltage-gated (Cav3/T-type)", "High"
    elif "ryanodine" in desc or "ryr" in desc:
        classification, confidence = "Ligand-gated (RyR)", "High"
    elif "ip3" in desc:
        classification, confidence = "Ligand-gated (IP3R)", "High"
    elif "cngc" in desc:
        classification, confidence = "Plant CNGC", "High"
    elif "tpc" in desc or "two-pore" in desc:
        classification, confidence = "Two-pore channel (TPC)", "High"
    elif m["voltage_sensor"] and m["ca_selectivity"]:
        classification, confidence = "Voltage-gated Ca2+ channel", "Medium"
    elif m["spry"]:
        classification, confidence = "Ligand-gated (likely RyR)", "Medium"
    elif m["ip3_binding"]:
        classification, confidence = "Ligand-gated (likely IP3R)", "Medium"
    elif m["cnbd"]:
        classification, confidence = "CNGC-like channel", "Medium"
    return {
        "classification": classification,
        "confidence": confidence,
        "length": len(sequence),
        "markers": m,
    }


def predict_electrophysiology(description: str, sequence: str) -> Dict[str, str]:
    """Gating, selectivity and pharmacology labels from the app's rules."""
    desc = description.lower()
    has_voltage = _has("Voltage_sensor", sequence)
    has_cnbd = _has("CNBD", sequence)
    has_glr = _has("GLR_ligand", sequence)

    if has_voltage:
        gating = "Voltage-gated"
        activation = "Depolarization-activated (-40 to +20 mV typical)"
    elif has_cnbd:
        gating = "Ligand-gated (cAMP/cGMP)"
        activation = "cAMP/cGMP binding"
    elif has_glr:
        gating = "Ligand-gated (Glutamate-like)"
        activation = "Glutamate or analog binding"
    elif "ryr" in desc or "ryanodine" in desc:
        gating = "Ligand-gated (Ca²⁺, cADPR, ryanodine)"
        activation = "Ca²⁺-induced Ca²⁺ release (CICR)"
    elif "ip3" in desc:
        gating = "Ligand-gated (IP3)"
        activation = "IP3 binding"
    else:
        gating = "Unknown/Mixed"
        activation = "Multiple mechanisms possible"

    if _has("Ca_selectivity", sequence):
        selectivity = "Ca²⁺-selective (PCa/PNa > 100)"
        permeability = "High Ca²⁺, Low Na⁺/K⁺"
    elif "tpc" in desc or "two-pore" in desc:
        selectivity = "Non-selective cation (Ca²⁺, Na⁺, K⁺)"
        permeability = "PCa:PNa:PK ≈ 1:1:1"
    else:
        selectivity = "Ca²⁺-permeable"
        permeability = "Variable Ca²⁺/monovalent selectivity"

    if "cav1" in desc or "l-type" in desc:
        blockers = "Dihydropyridines (nifedipine), phenylalkylamines (verapamil), benzothiazepines (diltiazem)"
    elif "cav2.1" in desc or "p/q" in desc:
        blockers = "ω-agatoxin IVA"
    elif "cav2.2" in desc or "n-type" in desc:
        blockers = "ω-conotoxin GVIA"
    elif "cav3" in desc or "t-type" in desc:
        blockers = "Mibefradil, TTA-P2"
    elif "ryr" in desc:
        blockers = "Ryanodine (high dose), dantrolene"
    else:
        blockers = "Family-specific blockers required"

    return {
        "gating": gating,
        "activation": activation,
        "selectivity": selectivity,
        "permeability": permeability,
        "blockers": blockers,
    }
