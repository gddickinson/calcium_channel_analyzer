"""FastAPI web app for the Calcium Channel Analyzer: serves the front end in static/ and the JSON API under /api.

Run with `python run_web.py` (http://127.0.0.1:8135, local only). `scripts/build_static.py` exports every path in
STATIC_ENDPOINTS to JSON for the GitHub Pages site, so keep that list in step with the GET routes below. The POST
analyses run the repo's Python code here; on the static site static_api.js runs the JS ports (seq_engine.js,
align_engine.js), which node tests check against this code. /api/fetch (NCBI, UniProt) exists only in the local app.
"""
from pathlib import Path

from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from . import analysis
from .fetch import fetch as fetch_sequences
from sequence_analysis import FastaError, TM_THRESHOLD, TM_WINDOW  # noqa: E402  (repo root on sys.path via analysis)

STATIC = Path(__file__).resolve().parent / "static"
TITLE = "Calcium Channel Analyzer"

app = FastAPI(title=TITLE)
app.mount("/static", StaticFiles(directory=STATIC), name="static")


class Fasta(BaseModel):
    fasta: str
    index: int = 0
    window: int = TM_WINDOW
    threshold: float = TM_THRESHOLD
    method: str = "UPGMA"
    use_min_length: bool = False


class FetchRequest(BaseModel):
    source: str
    accessions: str


def run(fn, *args):
    try:
        return fn(*args)
    except (FastaError, ValueError) as e:
        raise HTTPException(400, str(e))


@app.get("/", include_in_schema=False)
def index():
    return FileResponse(STATIC / "index.html")


@app.get("/api/info")
def info():
    """Project summary for the home page: counts from config.py and the example FASTA."""
    ex = analysis.parse(analysis.example_fasta())["records"]
    fams = analysis.CHANNEL_FAMILIES
    return {
        "title": TITLE, "version": "1",
        "n_domain_patterns": len(analysis.DOMAIN_PATTERNS),
        "n_families": sum(len(v) for v in fams.values()),
        "n_feature_dims": len(analysis.FEATURE_NAMES),
        "n_examples": len(ex),
        "example_residues": sum(r["length"] for r in ex),
    }


@app.get("/api/config")
def config():
    """Domain patterns, channel families, Kyte-Doolittle scale and analysis settings (config.py)."""
    return analysis.config_payload()


@app.get("/api/examples")
def examples():
    """The example FASTA shipped with the repo, raw and parsed."""
    text = analysis.example_fasta()
    return {"file": analysis.EXAMPLE_FASTA.name, "fasta": text, **analysis.parse(text)}


@app.post("/api/parse")
def parse(req: Fasta):
    return run(analysis.parse, req.fasta)


@app.post("/api/domains")
def domains(req: Fasta):
    return run(analysis.domains, req.fasta, req.use_min_length)


@app.post("/api/hydropathy")
def hydropathy(req: Fasta):
    return run(analysis.hydropathy, req.fasta, req.index, req.window, req.threshold)


@app.post("/api/classify")
def classify(req: Fasta):
    return run(analysis.classify, req.fasta)


@app.post("/api/align")
def align(req: Fasta):
    return run(analysis.align, req.fasta)


@app.post("/api/tree")
def tree(req: Fasta):
    return run(analysis.tree, req.fasta, req.method)


@app.post("/api/fetch")
def fetch(req: FetchRequest):
    """Fetch FASTA from NCBI protein or UniProt (local app only; needs the network)."""
    return run(fetch_sequences, req.source, req.accessions)


# GET paths (without the /api prefix) written to the static site as api/<path>.json.
STATIC_ENDPOINTS = ["/info", "/config", "/examples"]
