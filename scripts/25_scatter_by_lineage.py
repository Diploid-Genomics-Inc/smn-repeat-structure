# Pipeline step 25 (scatter_by_lineage); ported from explore/structure_vs_sequence/22_scatter_by_lineage.py (sha256 a50164c62e73, commit 5e2b3bb),
# paths rewritten to results/sv/ and figures/; run from the project root after the steps before it.
"""Exploration step 22 (not a pipeline step). Does the SMN1 lineage make-up of a haplotype (A, AA, AB, B, BB;
from out/lineage_haplotypes.tsv, step 06) shape the clustering of Fig. 1 (all pairs) and Fig. 2 (pairs of the same
arrangement type)?
Tests, for structural distance and sequence distance (all copies; SMN1 copies) separately:
  (1) all pairs: same make-up vs different make-up (permutation: make-up labels shuffled across haplotypes,
      and shuffled within arrangement type, so the type effect is held fixed);
  (2) within each arrangement type with >= 10 haplotypes and >= 2 make-ups of >= 2 haplotypes:
      same-make-up vs different-make-up pairs, labels shuffled within the type;
  (3) within each type, structure-vs-sequence Spearman over same-make-up pairs only, against all pairs of the type.
1999 shuffles, seed 1. Haplotypes with no SMN1 copy in the tree (4) have no make-up and are left out.
Outputs: out/lineage_scatter_tests.tsv, out/lineage_scatter_within.tsv, out/lineage_scatter_summary.json.
The figure is drawn by step 23."""
import numpy as np, pandas as pd, json
from scipy.stats import spearmanr
O = 'results/sv'; F = 'figures'
P = pd.read_csv(f'{O}/pair_distances.tsv.gz', sep='\t')
L = pd.read_csv(f'{O}/lineage_haplotypes.tsv', sep='\t')
hu = pd.read_csv('results/decomp/haplotype_units.tsv', sep='\t'); arr = dict(zip(hu.src, hu.arrangement))
comp = dict(zip(L.src, L.composition))
P['ca'] = P.a.map(comp); P['cb'] = P.b.map(comp); P['ta'] = P.a.map(arr); P['tb'] = P.b.map(arr)
n_drop = len((set(P.a) | set(P.b)) - set(comp))
P = P[P.ca.notna() & P.cb.notna()].copy().reset_index(drop=True)
haps = sorted(set(P.a) | set(P.b)); hi = {h: k for k, h in enumerate(haps)}
ia = P.a.map(hi).to_numpy(); ib = P.b.map(hi).to_numpy()
C = np.array([comp[h] for h in haps]); T = np.array([arr[h] for h in haps])
cols = {'struct': 'struct', 'seq_all': 'seq_all_copies', 'seq_SMN1': 'seq_SMN1'}
rng = np.random.default_rng(1); NS = 1999
def diff(lab, v, mask):
    same = lab[ia] == lab[ib]; m = mask & ~np.isnan(v)
    return np.mean(v[m & ~same]) - np.mean(v[m & same])
def shuffle_within(groups):
    lab = C.copy()
    for g in np.unique(groups):
        k = np.where(groups == g)[0]
        if len(k) > 1: lab[k] = C[rng.permutation(k)]
    return lab
rows = []; allm = np.ones(len(P), bool); same = C[ia] == C[ib]
for name, c in cols.items():
    v = P[c].to_numpy(float); ok = ~np.isnan(v); obs = diff(C, v, allm)
    for scheme, g in [('across all haplotypes', np.zeros(len(haps), int)), ('within arrangement type', T)]:
        ge = sum(diff(shuffle_within(g), v, allm) >= obs for _ in range(NS))
        rows.append(dict(scope='all pairs', distance=name, shuffle=scheme, pairs_same=int((same & ok).sum()),
            pairs_diff=int((~same & ok).sum()), median_same=np.median(v[same & ok]), median_diff=np.median(v[~same & ok]),
            mean_diff_minus_same=obs, p=(ge + 1) / (NS + 1)))
cnt = pd.Series(T).value_counts(); W = []
for t, n in cnt.items():
    if n < 10: continue
    k = T == t; cc = pd.Series(C[k]).value_counts(); tm = (T[ia] == t) & (T[ib] == t)
    rec = dict(arrangement=t, haplotypes=int(n), makeups=' '.join(f'{a}:{b}' for a, b in cc.items()))
    s = P.struct.to_numpy(float)
    for nm, c in [('seq_all', 'seq_all_copies'), ('seq_SMN1', 'seq_SMN1')]:
        g = P[c].to_numpy(float); m = tm & ~np.isnan(g)
        rec[f'r_{nm}_allpairs'] = spearmanr(s[m], g[m]).statistic
        ms = m & same; rec[f'pairs_same_{nm}'] = int(ms.sum())
        rec[f'r_{nm}_samemakeup'] = spearmanr(s[ms], g[ms]).statistic if ms.sum() >= 10 else np.nan
    if (cc >= 2).sum() >= 2:
        g = np.where(k, 0, np.arange(len(haps)) + 1)
        for name, c in cols.items():
            v = P[c].to_numpy(float); ok = ~np.isnan(v); obs = diff(C, v, tm)
            ge = sum(diff(shuffle_within(g), v, tm) >= obs for _ in range(NS))
            rec[f'{name}_median_same'] = np.median(v[tm & same & ok]); rec[f'{name}_median_diff'] = np.median(v[tm & ~same & ok])
            rec[f'{name}_p'] = (ge + 1) / (NS + 1)
    W.append(rec)
R = pd.DataFrame(rows); Wd = pd.DataFrame(W)
R.to_csv(f'{O}/lineage_scatter_tests.tsv', sep='\t', index=False, float_format='%.4g')
Wd.to_csv(f'{O}/lineage_scatter_within.tsv', sep='\t', index=False, float_format='%.4g')
sm = dict(haplotypes=len(haps), dropped_no_SMN1_in_tree=n_drop, pairs=len(P), shuffles=NS,
          makeups={k: int(v) for k, v in pd.Series(C).value_counts().items()})
json.dump(sm, open(f'{O}/lineage_scatter_summary.json', 'w'), indent=1)
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 40)
print(json.dumps(sm)); print(R.round(5).to_string(index=False)); print(Wd.round(4).T.to_string())
# Figures: drawn by step 23 (fig/hap_scatter_lineage.png), which consolidates step 22's two earlier figures.
