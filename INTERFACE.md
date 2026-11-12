# Calcium Channel Analyzer Interface Map

## Project Structure

```
calcium_channel_analyzer/
  calcium_channel_analyzer.py  -- Main GUI application (CalciumChannelAnalyzer class, Tkinter)
  config.py                    -- Shared constants: DOMAIN_PATTERNS, CHANNEL_FAMILIES,
                                  KD_HYDROPHOBICITY, AMINO_ACIDS
  domain_detection.py          -- Domain detection logic: detect_domains_in_sequence(),
                                  calculate_domain_score(), predict_family_from_description()
  ml_features.py               -- ML feature extraction: extract_ml_features()
  sequence_analysis.py         -- GUI-free: parse_fasta(), hydropathy_profile(), predict_tm_regions(),
                                  tm_prediction() (the GUI's predict_tm_domains logic)
  channel_rules.py             -- GUI-free: classify_channel(), predict_electrophysiology(), channel_markers()
                                  (the GUI's classify_channels / predict_electrophysiology rules)
  alignment_tree.py            -- GUI-free (Biopython): make_aligner(), pairwise_alignments() (the GUI's alignment
                                  settings); make_tree_aligner(), alignment_distances() (BLOSUM62 global alignments,
                                  1 - identity), build_tree() (UPGMA/NJ), tree_to_dict(), tree_to_newick()
  requirements.txt             -- Python dependencies (desktop app)
  requirements-web.txt         -- Pinned dependencies for the web app and its tests
  run_web.py                   -- Starts the web app on http://127.0.0.1:8135
  example_ca_channels.fasta    -- Example FASTA: 7 UniProtKB/Swiss-Prot calcium channels (see THIRD_PARTY.md)
  THIRD_PARTY.md               -- Source, date, checksum and licence (CC BY 4.0) of the example sequences
  DOCUMENTATION.md             -- User documentation
  web/
    app.py                     -- FastAPI: GET /api/info|config|examples (exported to the static site),
                                  POST /api/parse|domains|hydropathy|classify|align|tree, local-only POST /api/fetch
    analysis.py                -- FASTA text -> JSON bundles, calling the modules above; config_payload()
    fetch.py                   -- NCBI (Entrez) / UniProt REST FASTA fetching (local app only)
    static/                    -- Front end (no build step): index.html, core.js (router, api), main.js, state.js
                                  (current sequence set), ui.js, static_api.js (static site: runs the JS engines),
                                  pages home/sequences/domains/hydropathy/classify/align/tree/about .js,
                                  seq_engine.js + align_engine.js (JS ports) and their *.test.mjs node tests
  scripts/
    build_static.py            -- Export the static GitHub Pages site to site/
    screenshots.py             -- Playwright smoke test + README screenshots (docs/img/)
    make_engine_fixture.py     -- Writes tests/fixtures/engine_fixture.json (Python results for the JS tests)
  tests/
    test_domain_detection.py   -- Tests for domain detection and family prediction
    test_ml_features.py        -- Tests for ML feature extraction
    test_gui_free_modules.py   -- Tests for sequence_analysis, channel_rules, alignment_tree (+ drift guard vs GUI)
    test_tree_alignment.py     -- Alignment-based tree distances on two synthetic families with an indel
    test_engine_fixture.py     -- The JS fixture is current (rtol 1e-9)
    test_web.py                -- API, static build and node engine tests
  .github/workflows/pages.yml  -- Builds and publishes the static site on push to main
```

## Key Entry Points

- `python calcium_channel_analyzer.py` -- Launch the GUI
- `python run_web.py` -- Launch the web app
- `config.DOMAIN_PATTERNS` -- All domain regex patterns and metadata
- `domain_detection.detect_domains_in_sequence(seq, seq_id)` -- Detect domains
- `ml_features.extract_ml_features(seq)` -- Extract ML features from a sequence
- `sequence_analysis.tm_prediction(seq)` -- Hydropathy profile and TM regions

## Module Dependencies

- `config.py` -- standalone, no project imports
- `domain_detection.py`, `ml_features.py`, `sequence_analysis.py`, `channel_rules.py` -- import from `config`
- `alignment_tree.py` -- Biopython only
- `calcium_channel_analyzer.py` -- imports config, domain_detection, ml_features (its TM, classification, alignment
  and tree code is still inline; the GUI-free modules copy it, and tests/test_gui_free_modules.py checks the settings;
  its own Detect Domains and tree differ from the web's, see SESSION_LOG.md)
- `web/analysis.py` -- imports all the analysis modules; `web/app.py` imports `web/analysis.py` and `web/fetch.py`
