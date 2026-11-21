#!/bin/sh
# Pipeline step 02: multiple alignment of every counted SMN copy (SMN1, SMN2, h12, h21, del7-8; truncated pieces and
# the one fragment left out; copies joined from two pieces left out), each as gene body + 2 kb flanks on the SMN1-gene
# strand (data/derived/copy_seqs, from the clinical study's step 03). mafft FFT-NS-2 (~26 min).
# Ported from explore/structure_vs_sequence/01_gene_tree.sh; the fast IQ-TREE run there is dropped (step 28 builds
# the tree the note uses). Writes results/sv/copies.fa, copy_tips.tsv, copies.aln.fa.
set -e
D=results/sv
mkdir -p $D
python3 - <<'PY'
import pandas as pd
c = pd.read_csv('results/copies/copy_calls.tsv', sep='\t')
keep = c[c.copy_class.isin(['SMN1','SMN2','hyb_1_2','hyb_2_1','del7_8_SMN2like']) & (c.pieces == 1)]
seqs = {}
for src in sorted(keep.src.unique()):
    name = None
    for l in open(f'data/derived/copy_seqs/{src}.fa'):
        l = l.strip()
        if l.startswith('>'): name = l[1:]; seqs[name] = []
        else: seqs[name].append(l)
with open('results/sv/copies.fa', 'w') as fh:
    for i, r in enumerate(keep.sort_values('copy_id').itertuples()):
        fh.write(f'>c{i:04d}\n{"".join(seqs[r.copy_id])}\n')
keep.sort_values('copy_id').assign(tip=[f'c{i:04d}' for i in range(len(keep))])[
    ['tip','copy_id','src','copy_class','token']].to_csv('results/sv/copy_tips.tsv', sep='\t', index=False)
print(len(keep), 'copies')
PY
mafft --version 2>&1 | head -1
mafft --retree 2 --maxiterate 0 --thread 16 --quiet $D/copies.fa > $D/copies.aln.fa
