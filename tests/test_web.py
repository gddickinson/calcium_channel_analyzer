"""Web interface tests: every exported endpoint answers, the static build is complete, and the JS engine tests pass.

Run: python -m pytest -q tests/test_web.py
"""
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from web.app import STATIC, STATIC_ENDPOINTS, app  # noqa: E402

client = TestClient(app)


def test_index_served():
    r = client.get("/")
    assert r.status_code == 200 and "static/main.js" in r.text


@pytest.mark.parametrize("path", STATIC_ENDPOINTS)
def test_endpoint(path):
    r = client.get("/api" + path)
    assert r.status_code == 200, r.text
    assert r.json() is not None


def test_static_build(tmp_path):
    sys.path.insert(0, str(ROOT / "scripts"))
    import build_static
    build_static.build(tmp_path / "site")
    assert (tmp_path / "site" / "index.html").read_text().count("APP_STATIC") == 1
    for path in STATIC_ENDPOINTS:
        assert (tmp_path / "site" / "api" / (path.strip("/") + ".json")).exists()


def test_static_assets_referenced_exist():
    html = (STATIC / "index.html").read_text()
    for name in ("style.css", "main.js"):
        assert f"static/{name}" in html and (STATIC / name).exists()


@pytest.mark.skipif(not shutil.which("node"), reason="node not installed")
def test_js_engine():
    tests = sorted(STATIC.glob("*.test.mjs"))
    for t in tests:
        r = subprocess.run(["node", str(t)], capture_output=True, text=True)
        assert r.returncode == 0, r.stdout + r.stderr


# --- project-specific API tests ------------------------------------------------------------------------------------
EXAMPLE = (ROOT / "example_ca_channels.fasta").read_text()


@pytest.mark.parametrize("path", ["parse", "domains", "hydropathy", "classify", "align", "tree"])
def test_analysis_endpoints_on_example(path):
    r = client.post("/api/" + path, json={"fasta": EXAMPLE})
    assert r.status_code == 200, r.text
    assert r.json()


def test_example_counts():
    info = client.get("/api/info").json()
    assert info["n_examples"] == 7 and info["n_domain_patterns"] == 10 and info["n_feature_dims"] == 32


def test_bad_fasta_is_400():
    for body in ({"fasta": ""}, {"fasta": "MKVL"}, {"fasta": ">only header\n"}):
        assert client.post("/api/domains", json=body).status_code == 400


def test_tree_and_align_limits():
    one = ">a\nMKV\n>b\nMKL\n"
    assert client.post("/api/tree", json={"fasta": one}).status_code == 400
    assert client.post("/api/tree", json={"fasta": EXAMPLE, "method": "ML"}).status_code == 400
    many = "".join(f">s{i}\nMKVLA\n" for i in range(13))
    assert client.post("/api/align", json={"fasta": many}).status_code == 400


def test_hydropathy_bad_index_and_window():
    assert client.post("/api/hydropathy", json={"fasta": EXAMPLE, "index": 99}).status_code == 400
    assert client.post("/api/hydropathy", json={"fasta": EXAMPLE, "window": 1}).status_code == 400


def test_fetch_is_mocked_never_networked(monkeypatch):
    from web import fetch as fetch_mod
    monkeypatch.setitem(fetch_mod.FETCHERS, "uniprot", lambda acc: f">sp|{acc}|X test\nMKV\n")
    r = client.post("/api/fetch", json={"source": "uniprot", "accessions": "Q13936, BAD!"})
    assert r.status_code == 400  # invalid accession rejected before any request
    r = client.post("/api/fetch", json={"source": "uniprot", "accessions": "Q13936"})
    assert r.status_code == 200 and r.json()["fasta"].startswith(">sp|Q13936|")
    assert client.post("/api/fetch", json={"source": "pdb", "accessions": "1ABC"}).status_code == 400


def test_static_site_has_no_fetch_endpoint():
    assert all("fetch" not in p for p in STATIC_ENDPOINTS)
