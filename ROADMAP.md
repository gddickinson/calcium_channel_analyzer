# Calcium Channel Analyzer -- Roadmap

## Current State
A feature-rich bioinformatics platform for calcium channel analysis packed into a single monolithic file (`calcium_channel_analyzer.py`). Includes domain detection (10+ patterns), ML classification, AlphaFold3 integration, phylogenetic analysis, expression profiling, 11 visualization types, and a Tkinter GUI with 11 tabs. Auto-installs dependencies on first run. Has `DOCUMENTATION.md` and an example FASTA file. The single-file architecture is the main limitation -- this is likely thousands of lines.

## Short-term Improvements
- [x] Split `calcium_channel_analyzer.py` into modules: `config.py`, `domain_detection.py`, `ml_features.py` (extracted constants, domain logic, and ML feature extraction; main GUI file imports from these)
- [x] Add `requirements.txt` instead of auto-installing packages at runtime
- [x] Add unit tests for domain detection regex patterns against known channel sequences
- [ ] Add input validation for FASTA parsing -- handle malformed headers, empty sequences
- [x] Extract hardcoded domain patterns and channel family definitions into JSON/YAML config files
- [ ] Add proper error handling for network operations (NCBI, UniProt, AlphaFold API calls)
- [ ] Add logging instead of print statements for debugging API interactions

## Feature Enhancements
- [ ] Add protein multiple sequence alignment visualization (MSA viewer)
- [ ] Implement maximum likelihood phylogenetics (currently UPGMA only per benchmarks)
- [ ] Add InterPro/Pfam domain scanning for more accurate domain detection beyond regex
- [ ] Support protein structure superposition and RMSD calculation
- [ ] Add variant effect prediction for known channelopathy mutations
- [ ] Implement batch genome screening with progress tracking and result caching
- [ ] Add support for protein-protein docking predictions for channel complexes
- [ ] Create exportable analysis reports in PDF format

## Long-term Vision
- [ ] Refactor into a proper Python package with CLI and library API
- [ ] Build a web interface using Flask/Django for multi-user access
- [ ] Integrate with Galaxy bioinformatics platform as a tool
- [ ] Add deep learning-based function prediction (using ESM or ProtTrans embeddings)
- [ ] Create a curated calcium channel database with community curation
- [ ] Support custom HMM profiles trained on user-specific channel families

## Technical Debt
- [x] The entire application is in one file -- this is the single biggest issue to address
- [x] Auto-dependency installation at runtime is fragile and surprising -- replace with standard `pip install`
- [ ] Network API calls likely have no caching -- repeated queries hit external servers unnecessarily
- [x] No separation between business logic, data access, and GUI presentation
- [ ] ML classifier training likely happens in the GUI thread -- needs background threading
- [ ] `DOCUMENTATION.md` may not reflect the current codebase after changes -- needs audit
