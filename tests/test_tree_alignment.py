"""Distance trees from pairwise global alignments (BLOSUM62), on a small set with known relationships.

Two families, each a random 'ancestor' with mutated copies; one copy also has a short deletion near its N-terminus.
Distances from real alignments must put every copy with its own ancestor; the previous method (end-padding unaligned
sequences) loses the register after an indel and mixes the families.
"""
import random
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from alignment_tree import alignment_distances, build_tree, tree_to_dict  # noqa: E402

AA = "ACDEFGHIKLMNPQRSTVWY"


def mutate(seq, n, rng):
    s = list(seq)
    for i in rng.sample(range(len(s)), n):
        s[i] = rng.choice([a for a in AA if a != s[i]])
    return "".join(s)


def family_set():
    rng = random.Random(42)
    a = "".join(rng.choice(AA) for _ in range(150))
    b = "".join(rng.choice(AA) for _ in range(150))
    a2 = mutate(a, 12, rng)
    a2 = a2[:10] + a2[14:]  # 4-residue deletion near the start
    ids = ["A", "A1", "A2", "B", "B1", "B2"]
    seqs = [a, mutate(a, 8, rng), a2, b, mutate(b, 8, rng), mutate(b, 15, rng)[3:]]
    return ids, seqs


def leaves(node):
    return [node["name"]] if not node["children"] else [x for c in node["children"] for x in leaves(c)]


def clades(node):
    out = [frozenset(leaves(node))]
    for c in node["children"]:
        out += clades(c)
    return out


@pytest.mark.parametrize("method", ["UPGMA", "NJ"])
def test_families_recovered(method):
    ids, seqs = family_set()
    t = tree_to_dict(build_tree(ids, seqs, method).root)
    cl = set(clades(t))
    # each family forms a clade (for the unrooted NJ tree: one family is a clade and the other its complement)
    fam_a, fam_b = frozenset({"A", "A1", "A2"}), frozenset({"B", "B1", "B2"})
    assert fam_a in cl or fam_b in cl
    if method == "UPGMA":
        assert fam_a in cl and fam_b in cl


def test_distances_follow_alignment_not_register():
    ids, seqs = family_set()
    d, pairs = alignment_distances(ids, seqs)
    i = {k: n for n, k in enumerate(ids)}
    # the deletion copy is still close to its ancestor (about 12 substitutions + 4 gap columns in 150)
    assert d[i["A"]][i["A2"]] < 0.2
    assert d[i["B"]][i["B2"]] < 0.2
    for x in ("A", "A1", "A2"):
        for y in ("B", "B1", "B2"):
            assert d[i[x]][i[y]] > 0.6
    assert all(d[k][k] == 0 for k in range(len(ids)))
    assert all(abs(d[p][q] - d[q][p]) < 1e-12 for p in range(len(ids)) for q in range(len(ids)))
    assert len(pairs) == 15 and {"score", "identity", "length"} <= set(pairs[0])
