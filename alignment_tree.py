"""
Calcium Channel Analyzer -- Pairwise alignment and distance trees (GUI-free)
=============================================================================

The settings of the desktop app's "Multiple Alignment" and "Phylogenetic
Tree" actions, as pure functions (Biopython, no Tkinter):

  - make_aligner()          PairwiseAligner: global, match 1, mismatch 0,
                            gap open -0.5, gap extend -0.1
  - pairwise_alignments()   all pairs, each sequence cut to its first 500 aa
  - alignment_distances()   1 - identity of each pair's Needleman-Wunsch (Gotoh)
                            global alignment: BLOSUM62, gap open -10, gap extend
                            -0.5, end gaps scored like internal gaps; identity =
                            identical columns / alignment length
  - build_tree()            UPGMA or NJ (Bio.Phylo DistanceTreeConstructor) on
                            those distances
  - tree_to_dict()          nested {name, branch_length, children} for JSON

The desktop app's own tree pads unaligned sequences with '-' and compares them
position by position; that loses the register after any insertion or deletion,
so the web app and these functions use pairwise alignments instead. The
alignment used for the identity is Biopython's first optimal alignment
(alignments[0]); the JS port reproduces it (tests/fixtures).
"""

from typing import Dict, List, Sequence

from Bio import Phylo
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

ALIGN_PREFIX = 500
MATCH, MISMATCH, GAP_OPEN, GAP_EXTEND = 1.0, 0.0, -0.5, -0.1
TREE_MATRIX, TREE_GAP_OPEN, TREE_GAP_EXTEND = "BLOSUM62", -10.0, -0.5
MAX_CELLS = 40_000_000  # longest pair product for one alignment (RyR2 x RyR2 is about 25M)


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


def tree_matrix():
    return substitution_matrices.load(TREE_MATRIX)


def make_tree_aligner() -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = tree_matrix()
    aligner.open_gap_score = TREE_GAP_OPEN
    aligner.extend_gap_score = TREE_GAP_EXTEND
    return aligner


def matrix_letters(seq: str, alphabet: str) -> str:
    """Residues missing from the substitution matrix (e.g. U, O, J) are scored as X."""
    return "".join(c if c in alphabet else "X" for c in seq)


def alignment_distances(ids: Sequence[str], seqs: Sequence[str]):
    """
    Distance matrix (1 - identity) from the global alignment of every pair.

    Returns:
        (distances, pairs): full symmetric matrix, and per pair
        {i, j, score, identity, length} (length = alignment columns).
    """
    aligner = make_tree_aligner()
    alphabet = aligner.substitution_matrix.alphabet
    clean = [matrix_letters(s, alphabet) for s in seqs]
    n = len(seqs)
    d = [[0.0] * n for _ in range(n)]
    pairs = []
    for i in range(n):
        for j in range(i + 1, n):
            if len(clean[i]) * len(clean[j]) > MAX_CELLS:
                raise ValueError(f"{ids[i]} x {ids[j]} is too long to align (limit {MAX_CELLS:,} cells)")
            best = aligner.align(clean[i], clean[j])[0]
            top, bottom = str(best[0]), str(best[1])
            same = sum(x == y and x != "-" for x, y in zip(top, bottom))
            identity = same / len(top)
            d[i][j] = d[j][i] = 1.0 - identity
            pairs.append({"i": i, "j": j, "score": float(best.score), "identity": identity, "length": len(top),
                          "aligned1": top, "aligned2": bottom})
    return d, pairs


def build_tree(ids: Sequence[str], seqs: Sequence[str], method: str = "UPGMA", distances=None):
    """Bio.Phylo tree from alignment distances; method 'UPGMA' or 'NJ'."""
    if len(seqs) < 3:
        raise ValueError("Need at least 3 sequences for tree construction")
    if len(set(ids)) != len(ids):
        raise ValueError("Sequence ids must be unique to build a tree")
    if distances is None:
        distances, _ = alignment_distances(ids, seqs)
    dm = DistanceMatrix(list(ids), [distances[i][:i + 1] for i in range(len(ids))])
    constructor = DistanceTreeConstructor()
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
