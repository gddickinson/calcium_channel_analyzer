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
  requirements.txt             -- Pinned Python dependencies
  example_ca_channels.fasta    -- Example FASTA file for testing
  DOCUMENTATION.md             -- User documentation
  tests/
    test_domain_detection.py   -- Tests for domain detection and family prediction
    test_ml_features.py        -- Tests for ML feature extraction
```

## Key Entry Points

- `python calcium_channel_analyzer.py` -- Launch the GUI
- `config.DOMAIN_PATTERNS` -- All domain regex patterns and metadata
- `domain_detection.detect_domains_in_sequence(seq, seq_id)` -- Detect domains
- `ml_features.extract_ml_features(seq)` -- Extract ML features from a sequence

## Module Dependencies

- `config.py` -- standalone, no project imports
- `domain_detection.py` -- imports from `config`
- `ml_features.py` -- imports from `config`
- `calcium_channel_analyzer.py` -- imports from all above modules
