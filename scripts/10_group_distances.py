# Pipeline step 10 (group_distances); ported from explore/structure_vs_sequence/07_group_distances.py (sha256 6d2ce3c7d779, commit 5e2b3bb),
# paths rewritten to results/sv/ and figures/; run from the project root after the steps before it.
"""Exploration step 07: gene-body p-distance between every pair of copies of SMN1 lineage A, SMN1 lineage B and SMN2
(h21 hybrids as a side group). Same columns and distance as 03_lineage_check.py (gene body only; differences over the
bases both copies have). Lineages from doc/marker_lineage.tsv (B = in_marker_subclade).
Writes out/group_pair_distances.tsv (one row per pair of groups), results/sv/group_pair_distances_sub.tsv (B split into B1 and B2), out/group_copy_means.tsv (each copy's mean distance
to each group) and fig/group_distances.png/.pdf. Nothing here feeds the note."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42
O = 'results/sv'; F = 'figures'
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
A = A[:, col]; NCOL = int(col.sum())
ml = pd.read_csv(f'{O}/marker_lineage.tsv', sep='\t').set_index('tip')
cls = t.copy_class.replace({'hyb_2_1': 'h21', 'hyb_1_2': 'h12'})
g = []
for tip, c in zip(t.tip, cls):
    if c == 'SMN1' and tip in ml.index: g.append('B' if ml.loc[tip, 'in_marker_subclade'] else 'A')
    elif c in ('SMN2', 'h21'): g.append(c)
    else: g.append('x')
g = np.array(g); keep = g != 'x'
A = A[keep]; g = g[keep]; tk = t[keep].reset_index(drop=True)
V = np.isin(A, list(b'ACGT')).astype(np.float32); same = np.zeros((len(A), len(A)), np.float32)
for b in b'ACGT':
    M = (A == b).astype(np.float32); same += M @ M.T
n = V @ V.T; diff = n - same; D = diff / n * 1e3          # differences per kb
G = ['A', 'B', 'SMN2', 'h21']
idx = {k: np.where(g == k)[0] for k in G}
def pairs(a, b, X):
    s = X[np.ix_(idx[a], idx[b])]
    return s[np.triu_indices(len(idx[a]), 1)] if a == b else s.ravel()
rows = []
for i, a in enumerate(G):
    for b in G[i:]:
        d = pairs(a, b, D); k = pairs(a, b, diff); nn = pairs(a, b, n)
        rows.append(dict(a=a, b=b, n_a=len(idx[a]), n_b=len(idx[b]), pairs=len(d),
                         mean_per_kb=round(float(d.mean()), 4), median_per_kb=round(float(np.median(d)), 4),
                         q05_per_kb=round(float(np.quantile(d, .05)), 4), q95_per_kb=round(float(np.quantile(d, .95)), 4),
                         min_per_kb=round(float(d.min()), 4), max_per_kb=round(float(d.max()), 4),
                         mean_diffs=round(float(k.mean()), 1), mean_bases=int(round(float(nn.mean()))), columns=NCOL))
r = pd.DataFrame(rows)
m = {(x.a, x.b): x.mean_per_kb for x in r.itertuples()}; m.update({(b, a): v for (a, b), v in list(m.items())})
r['net_per_kb'] = [round(x.mean_per_kb - (m[x.a, x.a] + m[x.b, x.b]) / 2, 4) if x.a != x.b else np.nan for x in r.itertuples()]
ov = []
for x in r.itertuples():
    if x.a == x.b: ov.append(np.nan); continue
    thr = max(np.quantile(pairs(x.a, x.a, D), .95), np.quantile(pairs(x.b, x.b, D), .95))
    ov.append(round(float((pairs(x.a, x.b, D) < thr).mean()), 4))
r['share_below_within_q95'] = ov
r.to_csv(f'{O}/group_pair_distances.tsv', sep='\t', index=False)
print(r.to_string(index=False))
cm = pd.DataFrame({'tip': tk.tip, 'group': g})
for k in G:
    s = D[:, idx[k]].sum(1)
    cnt = np.full(len(g), len(idx[k]), float); cnt[idx[k]] -= 1   # a copy's distance to itself is 0; drop it from the count
    cm[f'mean_to_{k}'] = np.round(s / cnt, 3)
cm.to_csv(f'{O}/group_copy_means.tsv', sep='\t', index=False)
# sub-lineages of B (step 05: B1 = S1-9/S1-6 clade, B2 = S1-8 clade) against A, SMN2 and each other
sg = np.array([ml.loc[x, 'sublineage'] if gg == 'B' else gg for x, gg in zip(tk.tip, g)])
sidx = {k: np.where(sg == k)[0] for k in ['A', 'B1', 'B2', 'SMN2']}
sub_rows = []
for i, a in enumerate(['A', 'B1', 'B2', 'SMN2']):
    for b in ['A', 'B1', 'B2', 'SMN2'][i:]:
        X = D[np.ix_(sidx[a], sidx[b])]; d = X[np.triu_indices(len(sidx[a]), 1)] if a == b else X.ravel()
        sub_rows.append(dict(a=a, b=b, n_a=len(sidx[a]), n_b=len(sidx[b]), pairs=len(d), mean_per_kb=round(float(d.mean()), 4),
                             median_per_kb=round(float(np.median(d)), 4), min_per_kb=round(float(d.min()), 4), max_per_kb=round(float(d.max()), 4)))
pd.DataFrame(sub_rows).to_csv(f'{O}/group_pair_distances_sub.tsv', sep='\t', index=False)
print(pd.DataFrame(sub_rows).to_string(index=False))
import sys; sys.path.insert(0, 'scripts'); from smnlib import COPYCOL   # shared scheme: A red, B in B1's black, SMN2 blue, hybrids purple
C = {'A': COPYCOL['A'], 'B': COPYCOL['B1'], 'SMN2': COPYCOL['SMN2'], 'h21': COPYCOL['hybrid']}
lab = {'A': 'SMN1 A', 'B': 'SMN1 B', 'SMN2': 'SMN2', 'h21': 'h21'}
P6 = [('A', 'A'), ('B', 'B'), ('SMN2', 'SMN2'), ('A', 'B'), ('A', 'SMN2'), ('B', 'SMN2')]
fig = plt.figure(figsize=(11, 8.6))
ax = fig.add_axes([0.07, 0.57, 0.58, 0.37])
data = [pairs(a, b, D) for a, b in P6]
vp = ax.violinplot(data, showextrema=False, widths=0.85)
fc = [C['A'], C['B'], C['SMN2'], '#8a5a9a', '#4a8a8a', '#a07a3a']
for body, c in zip(vp['bodies'], fc): body.set_facecolor(c); body.set_alpha(0.6); body.set_edgecolor('k'); body.set_linewidth(0.4)
ytop = max(d.max() for d in data) + 0.45
for i, d in enumerate(data, 1):
    q1, med, q3 = np.quantile(d, [.25, .5, .75])
    ax.plot([i, i], [q1, q3], 'k-', lw=3); ax.plot(i, med, 'wo', ms=4, mec='k', mew=0.6)
    ax.text(i, d.max() + 0.05, f'{d.mean():.2f}', ha='center', va='bottom', fontsize=8)
ax.set_xticks(range(1, 7)); ax.set_xticklabels([f'{lab[a]}\nvs {lab[b]}\n({len(d):,} pairs)' for (a, b), d in zip(P6, data)], fontsize=8)
ax.set_ylabel('gene-body differences per kb'); ax.set_ylim(0, ytop)
ax.set_title('a  Every pair of copies (number above each = mean)', loc='left', fontsize=10)
ax.axvline(3.5, color='0.6', lw=0.6, ls=':')
ax.text(2, ytop * 0.98, 'within a group', ha='center', va='top', fontsize=8, color='0.35')
ax.text(5, ytop * 0.98, 'between groups', ha='center', va='top', fontsize=8, color='0.35')
ax2 = fig.add_axes([0.73, 0.59, 0.25, 0.33])
M3 = np.array([[m[a, b] for b in G] for a in G])
ax2.imshow(M3, cmap='viridis_r', vmin=0, vmax=M3.max())
for i in range(4):
    for j in range(4):
        ax2.text(j, i, f'{M3[i, j]:.2f}', ha='center', va='center', fontsize=8, color='w' if M3[i, j] > M3.max() * 0.55 else 'k')
tl = [f'{lab[k]}\n({len(idx[k])})' for k in G]
ax2.set_xticks(range(4)); ax2.set_xticklabels(tl, fontsize=8); ax2.set_yticks(range(4)); ax2.set_yticklabels(tl, fontsize=8)
ax2.set_title('b  Mean differences per kb', loc='left', fontsize=10)
for k_, (xg, yg, pos, ttl) in enumerate([('A', 'B', [0.07, 0.07, 0.37, 0.37], 'c  Each copy: mean distance to A and to B'),
                                          ('A', 'SMN2', [0.58, 0.07, 0.37, 0.37], 'd  Each copy: mean distance to A and to SMN2')]):
    a3 = fig.add_axes(pos)
    for k in G:
        s = cm[cm.group == k]
        a3.scatter(s[f'mean_to_{xg}'], s[f'mean_to_{yg}'], s=7, c=C[k], alpha=0.7, lw=0, label=f'{lab[k]} ({len(s)})')
    hi = max(cm[f'mean_to_{xg}'].max(), cm[f'mean_to_{yg}'].max()) + 0.1
    a3.plot([0, hi], [0, hi], color='0.6', lw=0.6, ls='--'); a3.set_xlim(0, hi); a3.set_ylim(0, hi)
    a3.set_xlabel(f'mean differences per kb to {lab[xg]} copies'); a3.set_ylabel(f'mean differences per kb to {lab[yg]} copies')
    a3.set_title(ttl, loc='left', fontsize=10)
    if k_ == 0: a3.legend(fontsize=8, markerscale=2, loc='lower right', frameon=False)
fig.savefig(f'{F}/group_distances.png', dpi=200)
fig.savefig(f'{F}/group_distances.pdf', metadata={'CreationDate': None})
print(cm.groupby('group')[[f'mean_to_{k}' for k in G]].mean().round(3))
print('columns', NCOL)
