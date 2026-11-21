#!/bin/sh
# Runs the preprint's pipeline from its inputs (results/bundles, results/copies, tables/, data/), in order.
# Steps 01-27 (and 04a, the ML tree, which step 05 cuts the lineages from) write results/sv/ and figures/; 29 onward are added as the note needs them; the numbers script and the
# note build come last. Run from the project root. ~35 min, most of it mafft (steps 02 and 03).
set -e
for s in scripts/[0-9][0-9]*_*; do
  echo "== $s"
  case $s in
    *.sh) sh "$s" ;;
    *.py) python3 scripts/trace_io.py "$s" ;;   # records the files each step opens (results/sv/provenance_io.tsv, read by step 41)
  esac
done
