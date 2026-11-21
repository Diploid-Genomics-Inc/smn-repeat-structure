# Pipeline step 08 (marker_group_distances); ported from explore/structure_vs_sequence/03_lineage_check.py (sha256 daba57c853f6, commit 5e2b3bb),
# paths rewritten to results/sv/ and figures/; run from the project root after the steps before it.
"""Exploration: mean gene-body p-distance between groups of copies, to check the SMN1 lineage that the UPGMA copy
tree places apart from both SMN1 and SMN2. Groups: SMN1 copies carrying c.*3+80T>G (results/copies/markers.tsv),
other SMN1, SMN2, h21. Same columns and copies as 02_compare.py."""
import numpy as np, pandas as pd
O = 'results/sv'
names, seqs = [], []
for l in open(f'{O}/copies.aln.fa'):
    l = l.strip()
    if l[0] == '>': names.append(l[1:]); seqs.append([])
    else: seqs[-1].append(l.upper())
A = np.array([np.frombuffer(''.join(s).encode(), dtype=np.uint8) for s in seqs])
t = pd.read_csv(f'{O}/copy_tips.tsv', sep='\t').set_index('tip').loc[names].reset_index()
ref = np.where((t.copy_class == 'SMN1') & t.src.str.startswith('GCA_000001405'))[0][0]
rpos = np.cumsum(A[ref] != ord('-')) - 1; L = rpos[-1] + 1
col = ((rpos >= 2000) & (rpos < L - 2000)) | ((A[ref] == ord('-')) & (rpos >= 1999) & (rpos < L - 2000))
keep = (t.copy_class != 'del7_8_SMN2like').to_numpy(); A = A[keep][:, col]; t = t[keep].reset_index(drop=True)
V = np.isin(A, list(b'ACGT')).astype(np.float32); same = np.zeros((len(A), len(A)), np.float32)
for b in b'ACGT':
    M = (A == b).astype(np.float32); same += M @ M.T
D = 1 - same / (V @ V.T)
m = pd.read_csv('results/copies/markers.tsv', sep='\t')

mk = dict(zip(m.copy_id, m.c3p80))
cls = t.copy_class.replace({'hyb_2_1': 'h21', 'hyb_1_2': 'h12'})
g = np.where((cls == 'SMN1') & (t.copy_id.map(mk) == 'alt'), 'SMN1 c.*3+80G', cls)
t['group'] = g
print(t.group.value_counts().to_string())
groups = ['SMN1 c.*3+80G', 'SMN1', 'SMN2', 'h21']
rows = []
for a in groups:
    for b in groups:
        ia, ib = np.where(g == a)[0], np.where(g == b)[0]
        sub = D[np.ix_(ia, ib)]
        v = sub[np.triu_indices(len(ia), 1)] if a == b else sub.ravel()
        rows.append({'a': a, 'b': b, 'mean_per_kb': round(float(v.mean() * 1e3), 3), 'median_per_kb': round(float(np.median(v) * 1e3), 3)})
r = pd.DataFrame(rows); r.to_csv(f'{O}/lineage_group_distances.tsv', sep='\t', index=False)
print(r.pivot(index='a', columns='b', values='mean_per_kb').to_string())
