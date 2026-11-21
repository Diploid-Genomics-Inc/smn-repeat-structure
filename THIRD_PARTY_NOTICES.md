# Third-party material in this package

Everything not listed here was produced by the study's own code.

## NCBI records and sequence (public domain)

Files under `data/outgroup/`, `data/refs/` and `data/refs_seq/` are records and sequence retrieved from NCBI
(E-utilities and nuccore). The accession and the coordinates are in each file's name, and `data/SOURCES.md` gives
the query behind each.

- `data/outgroup/*_SMN_*.fa`: great-ape SMN gene bodies from the T2T primary assemblies.
- `data/outgroup/regions/`, `data/outgroup/regions_hap2/`: the stretches around each ape's SMN gene.
- `data/refs_seq/*.gb`: RefSeq records for SMN1/SMN2, NAIP and GTF2H2.
- `data/refs/*.json`: NCBI Gene summaries.

NCBI places these records in the public domain; the underlying assemblies are those of their submitters.

## Values transcribed from published articles

- `data/papers_tables/luo2014_table4.tsv`: values read from Table 4 of Luo et al. 2014.
- `data/refs/berzal2026_fig4_psv.tsv`: the paralogous-sequence-variant names read from Fig. 4 of
  Berzal-Serrano et al. 2026.

Each holds only the values the study compares against, transcribed by hand, and is cited in the manuscript.

## Not redistributed here

- Chen et al. 2023 supplementary Table S6 (`mmc2.xlsx`), which step `11` reads: download it from the article.
- The publishers' PDFs of the cited papers.
- The HPRC release-2 assemblies: public, from the Human Pangenome Reference Consortium. The contigs this study
  read, with their GenBank accessions, are listed in `data/hprc_r2_contigs_used.tsv`.
