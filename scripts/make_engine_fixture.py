"""Write tests/fixtures/engine_fixture.json: Python results that the JS engines must reproduce.

Inputs are the repo's example FASTA plus a small hand-written edge-case FASTA (lower case, digits, stop codon, an
empty record, descriptions that reach each classification rule). Outputs come from web/analysis.py, i.e. the repo's
own Python modules. tests/test_engine_fixture.py checks this file is current; static/*.test.mjs read it.

Run: python scripts/make_engine_fixture.py
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from web import analysis  # noqa: E402

OUT = ROOT / "tests" / "fixtures" / "engine_fixture.json"

EDGE_FASTA = """>t1 test cyclic nucleotide-gated cngc sequence
mkrlavgllvkkrrrlavrkrkllffmmiivvlllggggaaaavvvddeeldkdldnletgflgiaaa*
>t2 empty record
>t3 test ryanodine-like
AAADAAADAAADAAAEAAAALIVLIVLIVMMMMMMMMMMMMLIVLIVKRKRVVVVVRKRKAAATLLGLGAAAA
>t4 unknown protein 12
MDEEKRLLLIVFAAAAKRILVFRKAADEAAEDLLSAFGIVLLAVILLLLAVVAALLIIFFVVLLGGAPPKK
DEGHIKLMNPQRSTVWY
>t5 test two-pore channel
GGGGSSSSTTTT
>t6 test ip3 receptor
KAALAAIAAK
"""


def bundle(fasta: str) -> dict:
    n = len(analysis.parse(fasta)["records"])
    return {
        "fasta": fasta,
        "parse": analysis.parse(fasta),
        "domains": analysis.domains(fasta),
        "classify": analysis.classify(fasta),
        "hydropathy": [analysis.hydropathy(fasta, i) for i in range(n)],
        "hydropathy_w11_t1": [analysis.hydropathy(fasta, i, 11, 1.0) for i in range(n)],
        "align": analysis.align(fasta),
        "tree_upgma": analysis.tree(fasta, "UPGMA"),
        "tree_nj": analysis.tree(fasta, "NJ"),
    }


def make() -> dict:
    return {"config": analysis.config_payload(),
            "cases": {"example": bundle(analysis.example_fasta()), "edge": bundle(EDGE_FASTA)}}


if __name__ == "__main__":
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(make(), ensure_ascii=False, indent=0))
    print(f"wrote {OUT} ({OUT.stat().st_size / 1e3:.0f} kB)")
