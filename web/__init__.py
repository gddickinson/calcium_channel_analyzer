"""Web interface for the Calcium Channel Analyzer (FastAPI + plain ES-module front end).

Contents:
    app.py       FastAPI app: GET /api/info, /api/config, /api/examples; POST analyses; local-only /api/fetch
    analysis.py  FASTA text -> JSON bundles, calling the repo's own analysis modules
    fetch.py     NCBI / UniProt sequence fetching (local app only)
    static/      front end (pages, JS engines ported from the Python code, node tests)
"""
