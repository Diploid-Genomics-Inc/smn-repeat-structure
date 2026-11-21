#!/bin/sh
# Pipeline step 01: structural distance between whole windows (principal-bundle composition and order, top 50
# bundles), pangea cohort bundle dist on the bundle decomposition results/decomp/smn.bed. Writes results/sv/smn_struct.*
# (from explore/structure_vs_sequence/README.md step 1).
set -e
mkdir -p results/sv
pangea --version
pangea cohort bundle dist -b results/decomp/smn.bed -o results/sv/smn_struct --top-n 50
