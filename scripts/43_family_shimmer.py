"""Pipeline step 43: sequence similarity between repeat-unit families, by shimmers (supplement S1).
Sample: for each family with >= 10 donor occurrences in results/decomp/smn.repeat_units.tsv, up to 12 units,
one per donor haplotype (seed 1); unit sequence from data/derived/windows/<row>.fa (step 28).
For each pair of units from different haplotypes, PANGEA shimmer matches (w24 k48 r2 min_span 12):
  n_i        distinct shimmer positions of unit i (self-match);
  shared     min(distinct matched positions in i, in j);
  contain    matched_i / n_i (directional: share of i's shimmers found in j);
  jaccard    shared / (n_i + n_j - shared)                       (whole units);
  overlap    shared / min(n_i, n_j)                              (normalised by the shorter unit);
  span_jacc  Jaccard restricted to the matched span: shimmers of i between its first and last matched
             position, likewise for j (normalised to the part the two units share).
Cell = median over pairs. Writes results/sv/family_shimmer_pairs.tsv, family_shimmer_<measure>.tsv,
family_shimmer.json and figures/family_shimmer.png/.pdf."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd, numpy as np, random, itertools, json, bisect
import pangeatk as ptk
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from smnlib import famlab
W, K, R, MS = 24, 48, 2, 12; N = 12; MINOCC = 10
O = 'results/sv'
u = pd.read_csv('results/decomp/smn.repeat_units.tsv', sep='\t', comment='#')
u = u[~u.row.str.startswith(('GCA_', 'chm13'))]
cnt = u.family.value_counts(); fams = [f for f in sorted(cnt.index) if cnt[f] >= MINOCC]
rnd = random.Random(1); units = []
for f in fams:
    rows = sorted(u[u.family == f].row.unique()); rnd.shuffle(rows)
    for r in rows[:N]: units.append((f, u[(u.family == f) & (u.row == r)].sample(1, random_state=1).iloc[0]))
cache = {}
def seq(x):
    if x.row not in cache: cache[x.row] = ''.join(l.strip() for l in open(f'data/derived/windows/{x.row}.fa') if not l.startswith('>'))
    return cache[x.row][int(x.bgn):int(x.end)].encode()
S = [seq(x) for _, x in units]
P = [sorted(set(ptk.get_shmmr_dots(s, s, W, K, R, MS)[0])) for s in S]
n = [len(p) for p in P]
def within(p, a, b): return bisect.bisect_right(p, b) - bisect.bisect_left(p, a)
rec = []
for i, j in itertools.combinations(range(len(units)), 2):
    if units[i][1].row == units[j][1].row or n[i] == 0 or n[j] == 0: continue
    xs, ys = ptk.get_shmmr_dots(S[i], S[j], W, K, R, MS)
    mi, mj = len(set(xs)), len(set(ys)); sh = min(mi, mj)
    if sh:
        si, sj = within(P[i], min(xs), max(xs)), within(P[j], min(ys), max(ys))
        sjac = sh / (si + sj - sh)
    else: sjac = 0.0
    base = dict(jaccard=sh / (n[i] + n[j] - sh), overlap=sh / min(n[i], n[j]), span_jacc=sjac)
    rec.append(dict(fam_a=units[i][0], fam_b=units[j][0], row_a=units[i][1].row, row_b=units[j][1].row, contain=mi / n[i], **base))
    rec.append(dict(fam_a=units[j][0], fam_b=units[i][0], row_a=units[j][1].row, row_b=units[i][1].row, contain=mj / n[j], **base))
d = pd.DataFrame(rec); d.round(4).to_csv(f'{O}/family_shimmer_pairs.tsv', sep='\t', index=False)
M = {m: d.groupby(['fam_a', 'fam_b'])[m].median().unstack().reindex(index=fams, columns=fams) for m in ['contain', 'jaccard', 'overlap', 'span_jacc']}
for m in ['jaccard', 'overlap', 'span_jacc']: assert np.allclose(M[m].values, M[m].values.T, equal_nan=True), m
for m, x in M.items(): x.round(3).to_csv(f'{O}/family_shimmer_{m}.tsv', sep='\t')
kb = {f: float(np.median([len(S[k]) for k in range(len(units)) if units[k][0] == f]) / 1e3) for f in fams}
npf = {f: sum(1 for k in units if k[0] == f) for f in fams}
def offmax(m):
    x = pd.DataFrame(M[m].to_numpy(copy=True), index=fams, columns=fams); v = np.array(x.to_numpy(), dtype=float, copy=True); np.fill_diagonal(v, np.nan); x = pd.DataFrame(v, index=fams, columns=fams); s = x.stack(); k = s.idxmax(); return dict(pair=list(k), value=round(float(s[k]), 3))
J = dict(families=fams, occurrences={f: int(cnt[f]) for f in fams}, units_per_family=npf, median_unit_kb={f: round(kb[f], 1) for f in fams},
         n_units=len(units), n_pairs=len(d) // 2, w=W, k=K, r=R, min_span=MS, max_units=N, min_occ=MINOCC,
         offdiag_max={m: offmax(m) for m in ['jaccard', 'overlap', 'span_jacc']},
         cells={m: {f'{a}|{b}': round(float(M[m].loc[a, b]), 3) for a in fams for b in fams} for m in M})
json.dump(J, open(f'{O}/family_shimmer.json', 'w'), indent=1, sort_keys=True)
fig, ax = plt.subplots(figsize=(7.2, 6))
X = M['jaccard'].values; im = ax.imshow(X, cmap='viridis', vmin=0, vmax=1)
for a in range(len(fams)):
    for b in range(len(fams)):
        v = X[a, b]
        if not np.isnan(v): ax.text(b, a, f'{v:.2f}', ha='center', va='center', fontsize=10, color='w' if v < 0.6 else 'k')
ax.set_xticks(range(len(fams))); ax.set_xticklabels([famlab(f) for f in fams], fontsize=11)
ax.set_yticks(range(len(fams))); ax.set_yticklabels([f'{famlab(f)} ({kb[f]:.0f} kb)' for f in fams], fontsize=11)
fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03, label='Jaccard similarity of shimmer sets')
fig.savefig('figures/family_shimmer.png', dpi=150, bbox_inches='tight', pad_inches=0.15)
fig.savefig('figures/family_shimmer.pdf', bbox_inches='tight', pad_inches=0.15, metadata={'CreationDate': None})
print('units', len(units), 'pairs', len(d) // 2)
for m in M: print(m); print(M[m].round(2).to_string())
print(json.dumps(J['offdiag_max']))
