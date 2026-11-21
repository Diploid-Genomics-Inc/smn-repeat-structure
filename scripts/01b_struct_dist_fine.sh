#!/bin/sh
# Pipeline step 01b: the structural distance at the FINER scale, for the supplement's two-scale comparison.
# Same command as step 01, on the clinical study's finer decomposition results/bundles/smn.bed
# (decomp -w 48 -r 4 --bundle-length-cutoff 2500 --bundle-merge-distance 10000, k 56; read-only copy from
# clinical commit 89cf028). Writes results/sv/fine/smn_struct.*. Its .dist is byte-identical to the one the
# preprint used before the switch (explore/structure_vs_sequence/out/smn_struct.dist).
# Checked 2026-11-15 (PANGEA 0.27.0): decomp with only the given settings (-w 48 -r 4 --bundle-length-cutoff 2500
# --bundle-merge-distance 10000; k 56, min-span 12, min-cov 0, min-branch-size 8 by default) on step 00's
# all_windows.fa reproduces results/bundles/smn.bed line for line. Every setting is listed by step 28i in
# results/sv/decomp_params.tsv.
set -e
mkdir -p results/sv/fine
pangea --version
pangea cohort bundle dist -b results/bundles/smn.bed -o results/sv/fine/smn_struct --top-n 50
