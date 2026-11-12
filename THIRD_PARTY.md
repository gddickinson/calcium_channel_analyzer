# Third-party data

## example_ca_channels.fasta -- UniProtKB/Swiss-Prot sequences

Seven reviewed protein sequences, copied unchanged (headers and residues) from UniProt:

| Accession | Entry | Protein | Organism | Length |
|---|---|---|---|---|
| Q13936 | CAC1C_HUMAN | Voltage-dependent L-type calcium channel subunit alpha-1C (Cav1.2, CACNA1C) | Homo sapiens | 2221 |
| O00555 | CAC1A_HUMAN | Voltage-dependent P/Q-type calcium channel subunit alpha-1A (Cav2.1, CACNA1A) | Homo sapiens | 2506 |
| Q92736 | RYR2_HUMAN | Ryanodine receptor 2 | Homo sapiens | 4967 |
| Q14643 | ITPR1_HUMAN | Inositol 1,4,5-trisphosphate-gated calcium channel ITPR1 | Homo sapiens | 2758 |
| O65718 | CNGC2_ARATH | Cyclic nucleotide-gated ion channel 2 | Arabidopsis thaliana | 726 |
| Q94KI8 | TPC1_ARATH | Two pore calcium channel protein 1 | Arabidopsis thaliana | 733 |
| Q9M8W7 | GLR11_ARATH | Glutamate receptor 1.1 | Arabidopsis thaliana | 808 |

- Source: UniProt REST API, `https://rest.uniprot.org/uniprotkb/<accession>.fasta` (canonical isoforms),
  UniProt release 2026_04, retrieved 2026-11-12. SHA-256 of the file:
  `a703d44a6327a0855b376896ba2f6ef4f87e087755e240e690270cd681abb2d1`.
- Licence: Creative Commons Attribution 4.0 International (CC BY 4.0), https://www.uniprot.org/help/license.
  Cite: The UniProt Consortium, "UniProt: the Universal Protein Knowledgebase", Nucleic Acids Research.
- Accessions verified against the entry names. The previous example file had headers naming these proteins but
  placeholder sequences (277-624 aa of homopolymer/repeat runs), and four of its accessions pointed at other proteins
  (P21457 = recoverin, Q9SRD3 = a vacuolar iron transporter, Q12001 = yeast ALG6; Q9LNQ0 does not exist). Those were
  replaced by Q92736, O65718, Q9M8W7 and Q94KI8.
- The file is fetched once by hand; tests and CI never contact UniProt.
