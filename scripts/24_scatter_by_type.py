# Pipeline step 24 (scatter_by_type); ported from explore/structure_vs_sequence/21_scatter_by_type.py (sha256 cb3e20b75efc, commit 5e2b3bb),
# paths rewritten to results/sv/ and figures/; run from the project root after the steps before it.
"""Exploration step 21 (not a pipeline step). Fig. 1 of the document redrawn with same-type pairs coloured by
arrangement type, and the within-type questions it raises:
  (1) within one arrangement type, do structural and sequence distance still go together (Spearman, Mantel)?
  (2) which types are more diverse, structurally and in sequence?
  (3) is the structural spread within a type driven by window length (and number of repeat units)?
Inputs: out/pair_distances.tsv.gz (02_compare.py), results/decomp/haplotype_units.tsv (arrangement, units),
results/copies/windows.tsv (window length; resolved haplotypes have one copy-bearing window).
Outputs: out/scatter_type_stats.tsv, out/scatter_type_summary.json, fig/hap_scatter_bytype.png/.pdf.
Types with >= 10 haplotypes are tested; Mantel: 999 shuffles of haplotype labels within the type, seed 1."""
import numpy as np, pandas as pd, json
from scipy.stats import spearmanr
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
O = 'results/sv'; F = 'figures'
P = pd.read_csv(f'{O}/pair_distances.tsv.gz', sep='\t')
hu = pd.read_csv('results/decomp/haplotype_units.tsv', sep='\t')
arr = dict(zip(hu.src, hu.arrangement)); nunit = dict(zip(hu.src, hu.units.fillna('').str.split().str.len()))
W = pd.read_csv('results/copies/windows.tsv', sep='\t'); W = W[W.n_copies > 0]
wlen = W.groupby('src').length.max().to_dict(); nwin = W.groupby('src').size().to_dict()
P['ta'] = P.a.map(arr); P['tb'] = P.b.map(arr)
assert P.ta.notna().all() and P.tb.notna().all()
assert (P.same_arrangement == (P.ta == P.tb)).all()
haps = sorted(set(P.a) | set(P.b)); cnt = pd.Series({h: arr[h] for h in haps}).value_counts()
types = [t for t, n in cnt.items() if n >= 10]
rng = np.random.default_rng(1)
def mat(sub, hs, col):
    i = {h: k for k, h in enumerate(hs)}; M = np.full((len(hs), len(hs)), np.nan)
    for a, b, v in sub[['a', 'b', col]].itertuples(index=False): M[i[a], i[b]] = M[i[b], i[a]] = v
    return M
def mantel(S, G, n=999):
    iu = np.triu_indices(len(S), 1); ok = ~np.isnan(G[iu])
    if ok.sum() < 10: return np.nan, np.nan
    r = spearmanr(S[iu][ok], G[iu][ok]).statistic; ge = 0
    for _ in range(n):
        p = rng.permutation(len(S)); Gp = G[np.ix_(p, p)]; ok2 = ~np.isnan(Gp[iu])
        ge += spearmanr(S[iu][ok2], Gp[iu][ok2]).statistic >= r
    return r, (ge + 1) / (n + 1)
rows = []
for t in types:
    hs = sorted(h for h in haps if arr[h] == t); sub = P[(P.ta == t) & (P.tb == t)]
    S = mat(sub, hs, 'struct'); G = mat(sub, hs, 'seq_all_copies'); G1 = mat(sub, hs, 'seq_SMN1')
    rA, pA = mantel(S, G); r1, p1 = mantel(S, G1)
    L = np.array([wlen[h] for h in hs], float); U = np.array([nunit[h] for h in hs], float)
    iu = np.triu_indices(len(hs), 1); dL = np.abs(L[:, None] - L[None, :])[iu]; dU = np.abs(U[:, None] - U[None, :])[iu]
    rL = spearmanr(S[iu], dL).statistic; rU = spearmanr(S[iu], dU).statistic if dU.std() > 0 else np.nan
    rLG = spearmanr(G[iu], dL).statistic
    # partial Spearman of structure vs sequence (all copies), window-length difference held fixed
    okp = ~np.isnan(G[iu]); rs = lambda x, y: spearmanr(x[okp], y[okp]).statistic
    a_, b_, c_ = rs(S[iu], G[iu]), rs(S[iu], dL), rs(G[iu], dL)
    rP = (a_ - b_ * c_) / np.sqrt((1 - b_ ** 2) * (1 - c_ ** 2))
    rows.append(dict(arrangement=t, haplotypes=len(hs), pairs=len(sub), multi_window=int(sum(nwin.get(h, 0) > 1 for h in hs)),
        struct_median=sub.struct.median(), struct_iqr=sub.struct.quantile(.75) - sub.struct.quantile(.25),
        seq_all_median_perkb=sub.seq_all_copies.median() * 1e3, seq_all_iqr_perkb=(sub.seq_all_copies.quantile(.75) - sub.seq_all_copies.quantile(.25)) * 1e3,
        seq_SMN1_median_perkb=sub.seq_SMN1.median() * 1e3,
        window_len_median_kb=np.median(L) / 1e3, window_len_iqr_kb=(np.quantile(L, .75) - np.quantile(L, .25)) / 1e3,
        units_median=np.median(U), units_range=f'{int(U.min())}-{int(U.max())}',
        mantel_all_r=rA, mantel_all_p=pA, mantel_SMN1_r=r1, mantel_SMN1_p=p1,
        struct_vs_lendiff_r=rL, struct_seq_partial_lendiff_r=rP, struct_vs_unitdiff_r=rU, seq_vs_lendiff_r=rLG))
T = pd.DataFrame(rows); T.to_csv(f'{O}/scatter_type_stats.tsv', sep='\t', index=False, float_format='%.4g')
# all pairs: same vs different type, and length
La = P.a.map(wlen); Lb = P.b.map(wlen); P['dL'] = (La - Lb).abs()
sm = dict(types_tested=len(types), haplotypes_in_tested=int(T.haplotypes.sum()), haplotypes=len(haps),
    same_pairs=int(P.same_arrangement.sum()), all_pairs=len(P),
    struct_vs_lendiff_all=round(float(spearmanr(P.struct, P.dL).statistic), 4),
    seq_vs_lendiff_all=round(float(spearmanr(P.seq_all_copies, P.dL).statistic), 4),
    struct_vs_lendiff_samepairs=round(float(spearmanr(P.struct[P.same_arrangement], P.dL[P.same_arrangement]).statistic), 4))
json.dump(sm, open(f'{O}/scatter_type_summary.json', 'w'), indent=1)
print(json.dumps(sm)); pd.set_option('display.width', 250); pd.set_option('display.max_columns', 30)
print(T.round(3).to_string(index=False))
# figure
pal = ['#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02', '#a6761d', '#1f78b4', '#b15928', '#6a3d9a', '#fb9a99', '#cab2d6']
col = {t: pal[i % len(pal)] for i, t in enumerate(types)}
fig = plt.figure(figsize=(14, 13)); gs = fig.add_gridspec(4, 4)
Q = lambda s, q: s.quantile(q)
for k, t in enumerate(types):
    ax = fig.add_subplot(gs[k // 4, k % 4])
    d = P[~P.same_arrangement]
    ax.scatter(d.struct, d.seq_all_copies * 1e3, s=1, c=[(.88, .88, .88)], alpha=.15, rasterized=True)
    d = P[(P.ta == t) & (P.tb == t)]
    ax.scatter(d.struct, d.seq_all_copies * 1e3, s=3, c=col[t], alpha=.6, rasterized=True)
    r = T.set_index('arrangement').loc[t]
    ax.set_title(f'{chr(97 + k)}  "{t}": {cnt[t]} haplotypes, {len(d)} pairs\nwithin-type Spearman {r.mantel_all_r:.2f} (Mantel p {r.mantel_all_p:.3g})', fontsize=8, loc='left')
    ax.set_xlim(-.02, .85); ax.set_ylim(-.05, 2.05)
    if k % 4 == 0: ax.set_ylabel('sequence distance, all copies\n(differences per kb)', fontsize=8)
    if k >= 3: ax.set_xlabel('structural distance', fontsize=8)
    ax.tick_params(labelsize=7)
ax = fig.add_subplot(gs[2:, :2])
for _, r in T.iterrows():
    d = P[(P.ta == r.arrangement) & (P.tb == r.arrangement)]
    x0, x1, xm = Q(d.struct, .25), Q(d.struct, .75), d.struct.median(); y = d.seq_all_copies * 1e3; y0, y1, ym = Q(y, .25), Q(y, .75), y.median()
    ax.errorbar(xm, ym, xerr=[[xm - x0], [x1 - xm]], yerr=[[ym - y0], [y1 - ym]], fmt='o', color=col[r.arrangement], ms=4 + r.haplotypes ** .5, capsize=2)
    ax.annotate(f'"{r.arrangement}"', (xm, ym), fontsize=8, xytext=(5, 5), textcoords='offset points')
ax.set_xlabel('within-type structural distance (median; bars: quartiles)'); ax.set_ylabel('within-type sequence distance, all copies\n(differences per kb; median, quartiles)')
ax.set_title(f'{chr(97 + len(types))}  diversity within each type (point size: haplotypes)', fontsize=9, loc='left')
ax = fig.add_subplot(gs[2:, 2:])
for t in types:
    d = P[(P.ta == t) & (P.tb == t)]
    ax.scatter(d.dL / 1e3, d.struct, s=3, c=col[t], alpha=.5, rasterized=True, label=f'"{t}"')
ax.legend(markerscale=4, fontsize=7, frameon=False, loc='lower right')
ax.set_xlabel('difference in window length between the two haplotypes (kb)'); ax.set_ylabel('structural distance')
ax.set_title(f'{chr(98 + len(types))}  same-type pairs: structural distance against window-length difference (Spearman {sm["struct_vs_lendiff_samepairs"]:.2f})', fontsize=9, loc='left')
fig.tight_layout(); fig.savefig(f'{F}/hap_scatter_bytype.png', dpi=200); fig.savefig(f'{F}/hap_scatter_bytype.pdf', metadata={'CreationDate': None}); plt.close(fig)
