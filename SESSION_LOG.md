# Session log

## 2026-11-12 -- Web interface and static site (branch `web-interface`)
- Added a FastAPI web app (`web/`, `run_web.py`, port 8135) and a static GitHub Pages site (`scripts/build_static.py`,
  `.github/workflows/pages.yml`) in the shared dark-neon web kit style. Pages: home, sequences, domain motifs,
  hydropathy & TM, classification & features, pairwise alignment, distance tree, about.
- Factored GUI-free modules out of the logic inline in `calcium_channel_analyzer.py`: `sequence_analysis.py`
  (FASTA parsing, KD hydropathy, TM regions), `channel_rules.py` (rule-based classification, electrophysiology labels),
  `alignment_tree.py` (PairwiseAligner settings, identity distances, UPGMA/NJ). The desktop file was NOT edited
  (George has uncommitted work on it); a test checks its constants still match the new modules.
- JS ports (`web/static/seq_engine.js`, `align_engine.js`) run every analysis on the static site; node tests compare
  them with `tests/fixtures/engine_fixture.json` from `scripts/make_engine_fixture.py` (rtol 1e-9). Alignment scores
  match Biopython; for co-optimal alignments the JS may show a different (equally scoring) alignment.
- NCBI / UniProt fetching is local-app only (`web/fetch.py`); tests mock it.
- Fixed README: install via `pip install -r requirements.txt` (the auto-install it described no longer exists) and
  removed the reference to the missing DISCOVERY_TOOLS_GUIDE.md.
- Found, not fixed (needs George): `detect_all_domains` and `classify_sequence` are called but never defined
  (variant prediction, variant scan, the create_viz domain architecture plot, batch analysis) -- that code has no menu
  entry, so it is unreachable today; `predict_ligand_binding` is defined twice (the second wins, same body);
  the GUI's Detect Domains does not apply `min_length` and reports 1-based positions, while
  `domain_detection.detect_domains_in_sequence` (used by the web) applies `min_length` with 0-based starts;
  `example_ca_channels.fasta` holds placeholder sequences (homopolymer/repeat runs), not the UniProt entries named.

## 2026-11-12 (later) -- Fixes requested by George (branch `web-interface`)
- Example FASTA replaced with the real UniProtKB/Swiss-Prot sequences of the seven proteins its headers named (four of
  the old accessions were wrong: P21457, Q9SRD3, Q12001, Q9LNQ0 -> Q92736, O65718, Q9M8W7, Q94KI8). Fetched once
  from the UniProt REST API (release 2026_04); source, checksum and CC BY 4.0 licence in THIRD_PARTY.md. No network in
  tests/CI.
- JS hydropathy now reproduces Python 3.12's compensated float sum() (real sequences exposed last-digit differences).
- domain_detection.detect_domains_in_sequence: 1-based inclusive positions (+ integer start/end) and every match kept,
  as in the GUI; min_length=... or use_pattern_min_length=True filter on request. Web page has a checkbox for it.
- Tree: distances are now 1 - identity of pairwise global alignments (BLOSUM62, gap -10/-0.5, end gaps penalised),
  not end-padded unaligned sequences. The JS Gotoh port uses Biopython's tie rule (epsilon 1e-6, M > X > Y) and
  reproduces Biopython's alignments character for character (fixture). Test: two synthetic families with an indel.
- Still differs in the desktop GUI (file not touched): its Detect Domains output format and its padded tree.
