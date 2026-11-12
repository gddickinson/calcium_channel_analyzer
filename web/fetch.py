"""Local-only sequence fetching for the web app (never used by the static site, tests or CI).

- NCBI protein: Bio.Entrez.efetch(db="protein", rettype="fasta"), as the desktop app's "Load from GenBank".
- UniProt: the UniProt REST FASTA endpoint (https://rest.uniprot.org/uniprotkb/<accession>.fasta).

NCBI asks for a contact email: set NCBI_EMAIL, otherwise the address used by the desktop app is sent.
"""
import os
import re
import urllib.request

ACCESSION = re.compile(r"^[A-Za-z0-9_.\-]{1,30}$")
MAX_ACCESSIONS = 10
TIMEOUT = 30
NCBI_EMAIL = os.environ.get("NCBI_EMAIL", "george.dickinson@uci.edu")


def _accessions(text: str) -> list:
    accs = [a.strip() for a in re.split(r"[,\s]+", text or "") if a.strip()]
    if not accs:
        raise ValueError("Enter at least one accession")
    if len(accs) > MAX_ACCESSIONS:
        raise ValueError(f"At most {MAX_ACCESSIONS} accessions at a time")
    bad = [a for a in accs if not ACCESSION.match(a)]
    if bad:
        raise ValueError("Not an accession: " + ", ".join(bad))
    return accs


def fetch_ncbi(acc: str) -> str:
    from Bio import Entrez
    Entrez.email = NCBI_EMAIL
    with Entrez.efetch(db="protein", id=acc, rettype="fasta", retmode="text") as h:
        return h.read()


def fetch_uniprot(acc: str) -> str:
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    req = urllib.request.Request(url, headers={"User-Agent": "calcium-channel-analyzer"})
    with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
        return r.read().decode("utf-8")


FETCHERS = {"ncbi": fetch_ncbi, "uniprot": fetch_uniprot}


def fetch(source: str, accessions: str) -> dict:
    """FASTA text for each accession; failures are reported per accession, not raised."""
    if source not in FETCHERS:
        raise ValueError("source must be ncbi or uniprot")
    fasta, errors = [], []
    for acc in _accessions(accessions):
        try:
            text = FETCHERS[source](acc).strip()
            if not text.startswith(">"):
                raise ValueError("no FASTA record returned")
            fasta.append(text)
        except Exception as e:  # network or not-found: tell the user which accession failed
            errors.append(f"{acc}: {e}")
    return {"fasta": "\n".join(fasta) + ("\n" if fasta else ""), "errors": errors}
