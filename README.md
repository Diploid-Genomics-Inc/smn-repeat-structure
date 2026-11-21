# The SMN region's large-scale repeat structure against SMN gene sequence (HPRC release 2)

This repository holds the code, the result tables and the figures behind the manuscript
*SMN1 lineages on a mirrored repeat* (`note/note.pdf`, built from `note/note.tex` and `note/supplement.tex`).

The study reads the SMN region whole in the resolved haplotype assemblies of HPRC release 2, decomposes each
haplotype's window into large repeat units, and compares that structure with the sequence of the SMN copies.

- Manuscript: `note/note.pdf` - main text, then the supplementary sections S1-S16.
- Every number the manuscript states is computed by `scripts/40_numbers.py` and inserted as a LaTeX macro
  (`note/numbers.tex`); `note/numbers.json` records, for each number, the line of code and the files it was read
  from, and `results/ledger.tsv` states the finding each result file supports.

[placeholder: repository URL]
[placeholder: licence - see LICENSE]
[placeholder: citation and preprint DOI]

## Layout

| Path | What it holds |
|---|---|
| `scripts/` | the pipeline, steps `00`-`43` in order, plus `smnlib.py` (shared code), `chainrule.py`, `trace_io.py` (records the files each step opens) and `run_all.sh` |
| `results/copies/`, `results/bundles/`, `tables/` | the inputs carried over from the source study (SMN copy calls, the finer bundle decomposition, per-haplotype tables); checksums in `CLINICAL_INPUTS.sha256` |
| `results/decomp/` | the bundle decomposition this study uses (`smn.bed`) |
| `results/sv/` | every result table and summary the manuscript reads |
| `results/ledger.tsv` | one row per result file: what it is, what it supports, where it is used |
| `figures/` | the manuscript's figures, PNG and PDF, each drawn by a numbered step |
| `note/` | the manuscript's sources, the generated numbers and tables, and the PDF |
| `data/refs/`, `data/refs_seq/`, `data/outgroup/` | reference records and great-ape sequences fetched from NCBI |
| `data/hprc_r2_contigs_used.tsv` | every contig the study read, with its GenBank accession and the window taken from it |
| `data/SOURCES.md` | where each input came from, and the tool versions used |
| `VERSION`, `SHA256SUMS` | the commit this package was built from, and the sha256 of every file in it |

## Reproducing

Run every command from the repository root. What you can rerun depends on what you have.

**1. The manuscript's numbers and tables, from the result tables.**

```
python3 scripts/40_numbers.py     # rewrites note/numbers.tex, note/numbers.json and the note/tab_*.tex tables
python3 scripts/41_provenance.py  # rewrites the provenance table and figure
```

`40_numbers.py` recomputes every number the manuscript states and asserts the claims it makes in words; it stops at
the first disagreement. Three of its inputs are not in this package, being large and derived: the copy alignment and
the alignment with the ape outgroups (`results/sv/copies.aln.fa`, `copies_with_outgroups.aln.fa`, steps `02` and
`03`) and the GRCh38 window (step `00`). Make those three first, or ask the authors for them.

Of the 59 steps, 35 have every input they read in this package and rerun as it stands; the other 24 read data that
is regenerated (the windows, the alignments, the structural distances) or re-fetched (the ape region sequences,
Chen et al.'s supplementary table). The PDF is built from `note/note.tex` with any LaTeX that has the packages it
loads, and does not depend on any of that.

**2. From the assemblies (needs the HPRC release-2 archive and the PANGEA binary).**

Set `STRATA_DATA_HPRC` to the directory holding `hprc_r2_full.agcrs` (and its `.mdbi`/`.mdbv`/`.midx`), then

```
sh scripts/run_all.sh
```

Step `00` cuts each haplotype's window out of the archive and decomposes it into bundles, writing about 790 MB into
`data/derived/windows/`; the later steps read it. The contigs and the window taken from each are listed in
`data/hprc_r2_contigs_used.tsv`, so the windows can be checked against the public assemblies. Steps `01` and `01b`
need PANGEA, steps `02` and `03` MAFFT, step `04a` IQ-TREE 2.

[placeholder: PANGEA macOS (Apple silicon) binary, version 0.27.0 - URL and SHA-256]

**The copy sequences the alignment reads.** Step `02` aligns one FASTA per haplotype of its SMN copies
(`data/derived/copy_seqs/`). Step `00b` states the rule that defines them - each copy is its window from
`win_start - 2000` to `win_end + 2000` (`results/copies/copy_calls.tsv`), reverse-complemented on the minus strand -
and rebuilds them from the windows (`--write`), or checks them against files already there (the default). Two kinds
of copy it cannot rebuild: 51 of the 872 aligned copies lie on 24 haplotypes that have no window of their own, and
HG00597 mat's SMN2 copy is joined from two pieces (step `02` leaves joined copies out). Rebuilding the alignment
exactly therefore needs the copy FASTAs of the source study; everything downstream of the alignment is in
`results/sv/`.

## What is not in this package

| Left out | Size | Why, and how to get it |
|---|---|---|
| `data/derived/windows/` | ~790 MB | written by step `00` from the HPRC archive |
| `data/derived/copy_seqs/` | ~27 MB | rebuilt by step `00b` from the windows, apart from the copies named above |
| `results/sv/copies.fa`, `copies.aln.fa`, `copies_with_outgroups.aln.fa`, `copies_body_og.fa` | ~106 MB | written by steps `02`-`04a` |
| `results/sv/smn_struct.dist`, `results/sv/fine/smn_struct.dist` | ~36 MB | written by steps `01` and `01b` (`pangea cohort bundle dist`) |
| `data/outgroup/regions/`, `data/outgroup/regions_hap2/` | ~61 MB | the 5-18 Mb NCBI stretches around each ape's SMN gene, read by steps `14`, `15` and `35`; re-fetch from NCBI with the accession and coordinates in each file name |
| `papers/supp/chen2023/mmc2.xlsx` | | Chen et al. 2023 supplementary Table S6, read by step `11`; download it from the article (PMC9943720) and put it at that path |
| publishers' PDFs of the cited papers | | subscription copies; obtain from the publisher |
| the project's working notes, slides, exploration folders and the parked side note | | internal to the study, not part of the manuscript |

## Environment

- Python 3.13 with the packages pinned in `requirements.lock`.
- PANGEA 0.27.0 (`pangea cohort bundle dist`, `pangea cohort bundle decomp`, `agc-rs get`) for steps `00`, `01`,
  `01b` and `34`.
- MAFFT v7.526 and IQ-TREE 2.4.0 for steps `02`, `03` and `04a`.
- The exact versions used are recorded in `data/SOURCES.md`.
