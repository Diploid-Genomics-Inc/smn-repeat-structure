# Inputs and tools of the preprint pipeline (scripts/)

## Inputs
- `results/bundles/`, `results/copies/`, `tables/` (six tables), `data/derived/copy_seqs/` (not in git),
  `data/refs_seq/smn1_smn2_refseq.gb`: copied read-only from the clinical study `smn-1113-longread-clinical` at its
  commit 89cf028 (`CLINICAL_SUMMARY.md`, checksums in `CLINICAL_INPUTS.sha256`). They come from the HPRC release-2
  assemblies ($STRATA_DATA_HPRC: hprc_r2_full.agcrs/.mdbi/.mdbv/.midx) through that study's steps 01-37.
- `data/outgroup/*_SMN_*.fa`: great-ape SMN gene bodies, NCBI nuccore efetch of the T2T primary assemblies
  (accessions and coordinates in the file names; NCBI Gene 461829, 100984482, 101148291, 100171813;
  `data/outgroup/esearch_gene.json`, `esummary_gene.json`).
- `data/outgroup/regions/`, `data/outgroup/regions_hap2/`: 5-18 Mb around each ape's SMN gene on both haplotypes
  (NCBI efetch; assembly records in `data/outgroup/assembly/`). The 18-MB gorilla file is not in git (size);
  step 15 needs it: refetch with the efetch URL in its name.
- `papers/supp/chen2023/`: Chen et al. 2023 supplementary tables (haplogroup variants, Table S6; names, Table S3).

## Tools (versions read on this machine, 2026-11-14)
- pangea 0.27.0 (pangea-bin, main:56818da6, release build, macOS aarch64): `pangea cohort bundle dist`.
- MAFFT v7.526 (2024/Apr/26): copy alignment (FFT-NS-2), outgroup addition (--add --keeplength).
- IQ-TREE 2.4.0 (MacOS ARM 64-bit, built Mar 8 2025).
- Python 3.13.13 with NumPy 2.4.4, SciPy 1.17.1, pandas 3.0.2, Matplotlib 3.10.8 (the project's interpreter;
  requirements.lock; two packages differ from the lock, so the environment is not pinned).
