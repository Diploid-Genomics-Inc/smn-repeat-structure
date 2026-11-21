#!/bin/sh
# Pipeline step 03: the four great-ape SMN gene bodies (data/outgroup/<ape>_SMN_<accession>_<start>-<end>[_rc].fa,
# NCBI nuccore; orangutan already reverse-complemented to the SMN1-gene strand) added as rows to the copy alignment,
# keeping its columns (mafft --add --keeplength; ape bases that would need new columns are dropped).
# Writes results/sv/outgroups.fa and results/sv/copies_with_outgroups.aln.fa (~4 min).
set -e
D=results/sv
for ape in bonobo chimp gorilla orangutan; do
  f=$(ls data/outgroup/${ape}_SMN_*.fa)
  echo ">og_$ape"; grep -v '^>' $f
done > $D/outgroups.fa
mafft --thread 4 --quiet --add $D/outgroups.fa --keeplength $D/copies.aln.fa > $D/copies_with_outgroups.aln.fa
