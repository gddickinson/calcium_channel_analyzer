# Advanced Calcium Channel Bioinformatics Platform

## 🧬 Overview

Bioinformatics tool for calcium channel research, incorporating the latest 2024-2025 methodologies including AlphaFold3 integration, modern sequence analysis, and domain detection.

**Developed by**: Dr. George Dickinson, UC Irvine, Department of Neurobiology & Behavior

## ✨ Key Features

### 🔬 Bioinformatics Integration
- **AlphaFold3 compatibility**: Prepare sequences for state-of-the-art structure prediction
- **Fast homology search**: Guidance toward DIAMOND/MMseqs2 (100x faster than BLAST)
- **Multi-database access**: NCBI, UniProt, AlphaFold Database

### 🧪 Comprehensive Analysis
- **Domain Detection**: EF-hand, CNBD, SPRY, voltage sensors, pore loops, and more
- **Family Classification**: Automatic identification of Cav1-3, RyR, IP3R, CNGC, GLR, TPC families
- **Phylogenetic Analysis**: UPGMA and Neighbor-Joining tree construction
- **Structural Prediction**: Transmembrane topology, hydropathy plots

### 📊 Visualizations
- Domain architecture maps
- Phylogenetic trees
- Hydropathy plots
- Interactive figures with matplotlib

### 📁 Export Options
- HTML reports
- Excel spreadsheets
- FASTA files
- Newick trees
- Publication-quality figures

## 🚀 Quick Start

### Installation
```bash
# The tool auto-installs dependencies on first run
# Simply execute:
python3 calcium_channel_analyzer.py
```

### Basic Workflow

#### 1. Load Sequences
```
File → Load FASTA
or
File → Load from GenBank
```

Try the example file: `example_ca_channels.fasta`

#### 2. Analyze Domains
```
Domain Analysis tab → Select domain types → Run Domain Detection
```

#### 3. Classify Channels
```
Analysis menu → Classify Channels
```

#### 4. Build Phylogenetic Tree
```
Phylogenetics tab → Select algorithm → Build Tree
```

#### 5. Predict Structure
```
Structure Prediction tab → Select sequence → Export for AlphaFold
```

## 📚 Example Analyses

### Example 1: Plant Calcium Channels
Load sequences:
- `NP_176967` (AtCNGC2 - Arabidopsis cyclic nucleotide-gated channel 2)
- `NP_568963` (AtCNGC4 - Arabidopsis cyclic nucleotide-gated channel 4)
- `Q9LNQ0` (AtTPC1 - Arabidopsis two-pore channel 1)

Detect domains and classify to identify:
- CNBD domains in CNGCs
- EF-hands in TPC1
- Calmodulin-binding sites
- Transmembrane topology

### Example 2: Comparative Evolution
Load orthologous Cav1.2 sequences:
- Human: `NP_000710`
- Mouse: `NP_031595`
- Zebrafish: `NP_571438`

Build phylogenetic tree and compare domain conservation.

### Example 3: Novel Channel Discovery
1. BLAST with known channel
2. Detect domains in hits
3. Classify predicted family
4. Export to AlphaFold3 for structure

## 🔬 Scientific Background

### Calcium Channel Families

#### Animal Channels
**Voltage-gated (Cav)**:
- **Cav1 (L-type)**: Long-lasting, cardiac/skeletal muscle
- **Cav2 (P/Q, N, R-type)**: Neurotransmitter release
- **Cav3 (T-type)**: Low-voltage activated

**Ligand-gated**:
- **RyR1-3**: Ryanodine receptors (ER Ca2+ release)
- **IP3R1-3**: IP3 receptors (intracellular signaling)

**Other**:
- **TRP**: Temperature, pain, sensation
- **P2X**: ATP-gated
- **Orai**: Store-operated Ca2+ entry

#### Plant Channels (Based on 2024-2025 Research)
- **CNGC**: Cyclic nucleotide-gated (pathogen defense, signaling)
- **GLR**: Glutamate receptor-like (stress responses)
- **TPC1**: Two-pore channel (vacuolar Ca2+ release)
- **MSL**: Mechanosensitive
- **MCA**: Mid1-complementing activity

### Recent Breakthroughs (2024-2025)

1. **AlphaFold3** (May 2024, Nobel Prize 2024):
   - Predicts all biomolecular interactions
   - 50%+ improvement over previous methods
   - Revolutionizing structural biology

2. **CNGC2-CNGC4 Complex** (June 2025):
   - Heteromeric channel formation
   - P2K1-mediated phosphorylation
   - Critical for plant immunity

3. **TPC Structure** (Ongoing research):
   - Atomic-resolution structures
   - Gating mechanisms revealed
   - Drug binding sites characterized

## 🛠️ Advanced Usage

### For Production Research

**Fast Sequence Searching**:
```bash
# Install DIAMOND (recommended)
conda install -c bioconda diamond

# Search (100x faster than BLAST)
diamond blastp -d nr -q query.faa -o results.txt --very-sensitive
```

**Multiple Sequence Alignment**:
```bash
# MAFFT (recommended for accuracy)
mafft --auto sequences.fasta > alignment.fasta

# MUSCLE (recommended for speed)
muscle -align sequences.fasta -output alignment.fasta
```

**AlphaFold3 Structure Prediction**:
1. Export sequence using tool
2. Visit: https://alphafoldserver.com/
3. Upload FASTA file
4. Download PDB structure
5. Analyze in PyMOL or ChimeraX

## 📖 Documentation

Complete documentation: `DOCUMENTATION.md`

Includes:
- Detailed feature descriptions
- Domain pattern explanations
- Workflow examples
- Troubleshooting guide
- Database resources
- Citation information

## 🎯 Key Capabilities

### Domain Detection
Identifies conserved motifs including:
- **EF-hand**: Ca2+-binding loops
- **CNBD**: Cyclic nucleotide binding
- **Voltage sensors**: S4 segments (R/K residues)
- **Pore loops**: Selectivity filters
- **SPRY domains**: RyR-specific
- **IP3-binding**: IP3R-specific
- **TM helices**: Hydrophobic segments
- **Ca2+ selectivity**: EEEE locus
- **CaM-binding**: IQ motifs
- **GLR ligand**: Plant-specific

### Phylogenetic Methods
- UPGMA: Ultrametric tree
- Neighbor-Joining: Unrooted tree
- Distance matrix calculation
- Bootstrap support (via external tools)

### Export Formats
- **HTML**: Comprehensive reports
- **Excel**: Tabular data
- **FASTA**: Sequence collections
- **Newick**: Phylogenetic trees
- **PNG/SVG**: High-quality figures

## 🌐 Resources

### Databases
- **NCBI**: https://www.ncbi.nlm.nih.gov/
- **UniProt**: https://www.uniprot.org/
- **AlphaFold DB**: https://alphafold.ebi.ac.uk/
- **Pfam**: https://pfam.xfam.org/
- **IUPHAR**: https://www.guidetopharmacology.org/

### Tools
- **AlphaFold Server**: https://alphafoldserver.com/
- **DIAMOND**: https://github.com/bbuchfink/diamond
- **MMseqs2**: https://github.com/soedinglab/MMseqs2
- **MAFFT**: https://mafft.cbrc.jp/
- **IQ-TREE**: http://www.iqtree.org/

## 🔍 Research Applications

### Use Cases
1. **Novel channel discovery**: Identify new calcium channels in genomes
2. **Evolutionary analysis**: Track channel evolution across species
3. **Structure-function**: Map domains to channel properties
4. **Drug discovery**: Identify binding sites for therapeutics
5. **Comparative genomics**: Analyze channel diversity
6. **Functional prediction**: Infer properties from sequence

### Example Research Questions
- What calcium channels exist in this plant species?
- How do plant and animal channels differ?
- What domains are conserved across evolution?
- Where do drugs bind to L-type channels?
- How did voltage sensing evolve?

## 🤝 Support

For questions or issues:
- Review `DOCUMENTATION.md`
- Check example workflows
- Contact: George Dickinson, UC Irvine

## 📄 License

Academic research use.
Developed at UC Irvine

## 🙏 Acknowledgments

Built with:
- **Biopython**: Sequence analysis framework
- **AlphaFold**: Structure prediction
- **BLAST**: Sequence similarity
- **Matplotlib/Seaborn**: Visualization
- **Pandas/NumPy**: Data analysis

Based on research from:
- UC Irvine (PIEZO1, calcium imaging)
- Plant calcium signaling (Demidchik, Dodd, Sanders)
- Voltage-gated channels (Catterall)
- AlphaFold (DeepMind)

---

**Version 2.0** | December 2024  
*Modern bioinformatics for calcium channel research*

## 🚦 System Requirements

- **Python**: 3.8 or higher
- **RAM**: 4GB minimum, 8GB recommended
- **Storage**: 1GB for software, more for databases
- **OS**: Linux, macOS, or Windows

## 📊 Performance Notes

- BLAST searches: Slow (use DIAMOND/MMseqs2 for production)
- Domain detection: Fast (< 1 second per sequence)
- Phylogenetic trees: Fast for small datasets (< 100 sequences)
- AlphaFold: Requires external service

## 🔮 Future Development

Planned features:
- [ ] Local DIAMOND integration
- [ ] AlphaFold structure retrieval API
- [ ] Machine learning classification
- [ ] 3D structure viewer
- [ ] Protein-ligand docking
- [ ] Electrophysiology prediction

---


