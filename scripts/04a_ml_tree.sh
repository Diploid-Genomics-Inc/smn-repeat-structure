#!/bin/sh
# Pipeline step 04a (was 28; moved so the lineages of step 05 can be cut from it): maximum-likelihood tree of the copies over the gene body, with support, rooted on the great apes.
# Rows: every copy of step 02's alignment except the del7-8 copies (misaligned, see step 07), plus the four ape SMN
# gene bodies of step 03. Columns: the gene body (step 07's rule: GRCh38 SMN1 copy's bases 2000 .. length-2000, with
# the insertion columns inside). IQ-TREE 2: GTR+G, 1000 ultrafast bootstraps (UFBoot2), seed 1, outgroup the apes.
# Writes results/sv/copies_body_og.fa (not in git, regenerable) and results/sv/mltree.* (treefile, contree, iqtree, log).
set -e
python3 - <<'PY'
import numpy as np, pandas as pd
O = 'results/sv'
names, seqs = [], []
for l in open(f'{O}/copies_with_outgroups.aln.fa'):
    l = l.strip()
    if l.startswith('>'): names.append(l[1:].split()[0]); seqs.append([])
    else: seqs[-1].append(l.upper())
A = np.array([np.frombuffer(''.join(s).encode(), dtype=np.uint8) for s in seqs])
t = pd.read_csv(f'{O}/copy_tips.tsv', sep='\t').set_index('tip')
cls = [t.copy_class.get(n, 'outgroup') for n in names]
ref = [i for i, n in enumerate(names) if cls[i] == 'SMN1' and t.src[n].startswith('GCA_000001405')][0]
rpos = np.cumsum(A[ref] != ord('-')) - 1; rlen = rpos[-1] + 1
inside = (rpos >= 2000) & (rpos < rlen - 2000); col = inside | ((A[ref] == ord('-')) & (rpos >= 1999) & (rpos < rlen - 2000))
keep = [i for i, c in enumerate(cls) if c != 'del7_8_SMN2like']
assert col.sum() == 28778, "gene-body columns as in step 07"
with open(f'{O}/copies_body_og.fa', 'w') as fh:
    for i in keep: fh.write(f'>{names[i]}\n{A[i, col].tobytes().decode()}\n')
print(len(keep), 'rows', int(col.sum()), 'columns', sum(c == 'outgroup' for c in cls), 'outgroups')
PY
iqtree2 --version | head -1
iqtree2 -s results/sv/copies_body_og.fa -m GTR+G -B 1000 -T 8 -seed 1 -o og_orangutan --prefix results/sv/mltree -quiet
