"""
Calcium Channel Analyzer -- Pairwise alignment and distance trees (GUI-free)
=============================================================================

The settings of the desktop app's "Multiple Alignment" and "Phylogenetic
Tree" actions, as pure functions (Biopython, no Tkinter):

  - make_aligner()          PairwiseAligner: global, match 1, mismatch 0,
                            gap open -0.5, gap extend -0.1
  - pairwise_alignments()   all pairs, each sequence cut to its first 500 aa
  - identity_distances()    sequences padded with '-' to equal length (no MSA),
                            Bio.Phylo DistanceCalculator('identity')
  - build_tree()            UPGMA or NJ (Bio.Phylo DistanceTreeConstructor)
  - tree_to_dict()          nested {name, branch_length, children} for JSON

As in the desktop app, the tree is built from unaligned, end-padded
sequences, so it reflects position-by-position identity from the N-terminus,
not homology from a multiple sequence alignment.
"""

from typing import Dict, List, Sequence

from Bio import Phylo
from Bio.Align import MultipleSeqAlignment, PairwiseAligner
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

ALIGN_PREFIX = 500
MATCH, MISMATCH, GAP_OPEN, GAP_EXTEND = 1.0, 0.0, -0.5, -0.1


def make_aligner() -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = MATCH
    aligner.mismatch_score = MISMATCH
    aligner.open_gap_score = GAP_OPEN
    aligner.extend_gap_score = GAP_EXTEND
    return aligner


def pairwise_alignments(ids: Sequence[str], seqs: Sequence[str], prefix: int = ALIGN_PREFIX) -> List[Dict]:
    """Global alignment of every pair (i < j) of the first `prefix` residues."""
    aligner = make_aligner()
    out = []
    for i in range(len(seqs)):
        for j in range(i + 1, len(seqs)):
            a, b = seqs[i][:prefix], seqs[j][:prefix]
            best = aligner.align(a, b)[0]
            top, bottom = str(best[0]), str(best[1])
            identical = sum(x == y and x != "-" for x, y in zip(top, bottom))
            out.append({
                "i": i, "j": j, "seq1": ids[i], "seq2": ids[j],
                "score": float(best.score),
                "aligned1": top, "aligned2": bottom,
                "identity": identical / len(top) if top else 0.0,
            })
    return out


def padded_alignment(ids: Sequence[str], seqs: Sequence[str]) -> MultipleSeqAlignment:
    max_len = max(len(s) for s in seqs)
    return MultipleSeqAlignment(
        [SeqRecord(Seq(s + "-" * (max_len - len(s))), id=i) for s, i in zip(seqs, ids)]
    )


def identity_distances(ids: Sequence[str], seqs: Sequence[str]) -> List[List[float]]:
    """Full symmetric identity-distance matrix of the end-padded sequences."""
    dm = DistanceCalculator("identity").get_distance(padded_alignment(ids, seqs))
    n = len(ids)
    return [[float(dm[i, j]) for j in range(n)] for i in range(n)]


def build_tree(ids: Sequence[str], seqs: Sequence[str], method: str = "UPGMA"):
    """Bio.Phylo tree from identity distances; method 'UPGMA' or 'NJ'."""
    if len(seqs) < 3:
        raise ValueError("Need at least 3 sequences for tree construction")
    if len(set(ids)) != len(ids):
        raise ValueError("Sequence ids must be unique to build a tree")
    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(padded_alignment(ids, seqs))
    constructor = DistanceTreeConstructor(calculator)
    return constructor.upgma(dm) if method.upper() == "UPGMA" else constructor.nj(dm)


def tree_to_dict(clade) -> Dict:
    """Nested dict of a Bio.Phylo clade (children in Biopython's order)."""
    return {
        "name": clade.name,
        "branch_length": float(clade.branch_length or 0.0),
        "children": [tree_to_dict(c) for c in clade.clades],
    }


def tree_to_newick(tree) -> str:
    from io import StringIO
    buf = StringIO()
    Phylo.write(tree, buf, "newick")
    return buf.getvalue().strip()
