# Pipeline step 26 (scatter_lineage_haplogroup); ported from explore/structure_vs_sequence/23_scatter_lineage_haplogroup.py (sha256 e3ecce372474, commit 5e2b3bb),
# paths rewritten to results/sv/ and figures/; run from the project root after the steps before it.
"""Exploration step 23 (not a pipeline step). Do Chen 2023's SMN1 haplogroups explain the clusters left inside
arrangement types once lineage A/B is accounted for (step 22)? And one consolidated figure for steps 22-23.
Haplogroup make-up of a haplotype = its SMN1 copies' haplogroup families (out/chen_hg_copies.tsv, step 08),
sorted and joined with '+' (e.g. 'S1-8+S1-9'); a copy tied across families is 'tie'; the 2 S1-6 copies of lineage A are 'S1-6(A)'.
Tests, within each arrangement type with >= 10 haplotypes and >= 2 make-ups of >= 2 haplotypes:
  same-make-up vs different-make-up pairs (medians; permutation of make-up labels within the type, 1999, seed 1)
  for structural distance and sequence distance (all copies), and Spearman structure-vs-sequence over
  same-make-up pairs, beside all pairs and same-lineage-make-up pairs (step 22).
Outputs: out/hg_scatter_within.tsv, out/hg_haplotypes.tsv, out/hg_scatter_summary.json,
fig/hap_scatter_lineage.png/.pdf (replaces step 22's two figures)."""
import numpy as np, pandas as pd, json, textwrap
from scipy.stats import spearmanr
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
O = 'results/sv'; F = 'figures'
P = pd.read_csv(f'{O}/pair_distances.tsv.gz', sep='\t')
L = pd.read_csv(f'{O}/lineage_haplotypes.tsv', sep='\t')
H = pd.read_csv(f'{O}/chen_hg_copies.tsv', sep='\t'); H['src'] = H.copy_id.str.split('|').str[0]
# S1-6 copies fall in both lineages (10 B, 2 A): label the A ones apart so a haplogroup make-up implies a lineage make-up
H.loc[(H.family == 'S1-6') & (H.lineage == 'A'), 'family'] = 'S1-6(A)'
assert (H.groupby('family').lineage.nunique() == 1).all()
hgm = H.groupby('src').family.apply(lambda s: '+'.join(sorted(s))).to_dict()
hu = pd.read_csv('results/decomp/haplotype_units.tsv', sep='\t'); arr = dict(zip(hu.src, hu.arrangement))
comp = dict(zip(L.src, L.composition))
P = P[P.a.isin(comp) & P.b.isin(comp)].copy().reset_index(drop=True)
haps = sorted(set(P.a) | set(P.b)); hi = {h: k for k, h in enumerate(haps)}
assert all(h in hgm for h in haps)
ia = P.a.map(hi).to_numpy(); ib = P.b.map(hi).to_numpy()
C = np.array([comp[h] for h in haps]); G = np.array([hgm[h] for h in haps]); T = np.array([arr[h] for h in haps])
pd.DataFrame(dict(src=haps, arrangement=T, lineage_makeup=C, hg_makeup=G)).to_csv(f'{O}/hg_haplotypes.tsv', sep='\t', index=False)
sameC = C[ia] == C[ib]; sameG = G[ia] == G[ib]
assert not np.any(sameG & ~sameC)  # same haplogroup make-up implies same lineage make-up
rng = np.random.default_rng(1); NS = 1999
s = P.struct.to_numpy(float); q = P.seq_all_copies.to_numpy(float) * 1e3
cnt = pd.Series(T).value_counts(); types10 = [t for t, n in cnt.items() if n >= 10]
W = []
for t in types10:
    k = np.where(T == t)[0]; tm = (T[ia] == t) & (T[ib] == t); mk = pd.Series(G[k]).value_counts()
    rec = dict(arrangement=t, haplotypes=len(k), hg_makeups=' '.join(f'{a}:{b}' for a, b in mk.items()),
               pairs=int(tm.sum()), pairs_same_hg=int((tm & sameG).sum()), pairs_same_lineage=int((tm & sameC).sum()),
               r_all=spearmanr(s[tm], q[tm]).statistic)
    for nm, sm in [('lineage', sameC), ('hg', sameG)]:
        m = tm & sm
        rec[f'r_same_{nm}'] = spearmanr(s[m], q[m]).statistic if m.sum() >= 10 and np.ptp(s[m]) > 0 else np.nan
    testable = (mk >= 2).sum() >= 2
    for dn, v in [('struct', s), ('seq', q)]:
        rec[f'{dn}_median_same_hg'] = np.median(v[tm & sameG]) if (tm & sameG).any() else np.nan
        rec[f'{dn}_median_diff_hg'] = np.median(v[tm & ~sameG]) if (tm & ~sameG).any() else np.nan
        if testable:
            def d(lab):
                sg = lab[ia] == lab[ib]; return v[tm & ~sg].mean() - v[tm & sg].mean()
            obs = d(G); ge = 0
            for _ in range(NS):
                lab = G.copy(); lab[k] = G[rng.permutation(k)]; ge += d(lab) >= obs
            rec[f'{dn}_p_hg'] = (ge + 1) / (NS + 1)
        else: rec[f'{dn}_p_hg'] = np.nan
    # remaining spread: median sequence distance among same-hg pairs, by make-up
    W.append(rec)
Wd = pd.DataFrame(W); Wd.to_csv(f'{O}/hg_scatter_within.tsv', sep='\t', index=False, float_format='%.4g')
ST = (T[ia] == T[ib])
summ = dict(haplotypes=len(haps), same_type_pairs=int(ST.sum()), same_type_same_hg=int((ST & sameG).sum()),
            same_type_same_lineage=int((ST & sameC).sum()), shuffles=NS,
            seq_median_same_type_same_hg=float(np.median(q[ST & sameG])), seq_median_same_type_diff_hg=float(np.median(q[ST & ~sameG])),
            struct_median_same_type_same_hg=float(np.median(s[ST & sameG])), struct_median_same_type_diff_hg=float(np.median(s[ST & ~sameG])))
json.dump(summ, open(f'{O}/hg_scatter_summary.json', 'w'), indent=1)
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 40)
print(json.dumps(summ)); print(Wd.round(4).T.to_string())
# ---- consolidated figure ----
def klass(x, y):
    st = {x, y}
    if st <= {'A', 'AA'}: return 'both A-only'
    if x == y: return f'both {x}'
    if st & {'A', 'AA'}: return 'A-only vs B-carrying'
    return 'two different B-carrying'
kl = np.array([klass(C[i], C[j]) for i, j in zip(ia, ib)])
kc = {'both A-only': (.8, .8, .8), 'A-only vs B-carrying': '#4393c3', 'both AB': '#d6604d', 'both B': '#f4a582',
      'both BB': '#67001f', 'two different B-carrying': '#8073ac'}
top = pd.Series(G).value_counts(); top = [g for g in top.index if top[g] >= 5 and g != 'tie']
pal = ['#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02', '#a6761d', '#17becf', '#666666']
hc = {g: pal[i % len(pal)] for i, g in enumerate(top)}
fig = plt.figure(figsize=(21, 15)); gs = fig.add_gridspec(3, 4)
axs = [fig.add_subplot(gs[r, c]) for r in range(3) for c in range(4)]
YL = 'sequence distance, all copies (per kb)'
def by_lineage(ax, m, ycol, title, small):
    for kk in kc:
        e = m & (kl == kk)
        if e.any():
            a0 = kk == 'both A-only'
            ax.scatter(s[e], (P[ycol].to_numpy(float) * 1e3)[e], s=1.5 if (a0 and small) else 3, c=[kc[kk]],
                       alpha=.3 if a0 else .75, label=f'{kk} ({int(e.sum())})', rasterized=True)
    ax.set_title(title, loc='left', fontsize=10); ax.set_xlabel('structural distance')
    ax.legend(fontsize=6.5, markerscale=3, title='pair class by SMN1 lineage (pairs)', title_fontsize=6.5, loc='upper left')
def by_hg(ax, m, title):
    d = m & ~sameG; dl = d & (kl == 'A-only vs B-carrying')
    ax.scatter(s[d & ~dl], q[d & ~dl], s=2.5, c=[(.82, .82, .82)], alpha=.5, label=f'different haplogroup, same lineage class ({int((d & ~dl).sum())})', rasterized=True)
    ax.scatter(s[dl], q[dl], s=2.5, c=['#4393c3'], alpha=.35, label=f'different, A-only vs B-carrying ({int(dl.sum())})', rasterized=True)
    for g in list(hc) + ['other']:
        e = m & sameG & (np.isin(G[ia], list(hc), invert=True) if g == 'other' else (G[ia] == g))
        if e.any(): ax.scatter(s[e], q[e], s=5, c=[hc.get(g, 'k')], alpha=.8, label=f'both {g} ({int(e.sum())})', rasterized=True)
    ax.set_title(title, loc='left', fontsize=10); ax.set_xlabel('structural distance'); ax.set_ylabel(YL)
    ax.legend(fontsize=6, markerscale=2, title='pair class by SMN1 haplogroup make-up (pairs)', title_fontsize=6, loc='upper left')
allm = np.ones(len(P), bool)
by_lineage(axs[0], allm, 'seq_all_copies', f'a  all pairs ({len(P)}), by lineage', True); axs[0].set_ylabel(YL)
by_lineage(axs[1], allm, 'seq_SMN1', 'b  all pairs, SMN1 copies only, by lineage', True); axs[1].set_ylabel('sequence distance, SMN1 copies (per kb)')
by_lineage(axs[2], ST, 'seq_all_copies', f'c  same-type pairs ({int(ST.sum())}), by lineage', False); axs[2].set_ylabel(YL)
by_hg(axs[3], ST, f'd  same-type pairs, by haplogroup make-up')
for j, t in enumerate(types10):
    k = T == t; mk = pd.Series(G[k]).value_counts()
    by_hg(axs[4 + j], (T[ia] == t) & (T[ib] == t),
          f'{"efghijk"[j]}  "{t}": {int(k.sum())} haplotypes\n' + '\n'.join(textwrap.wrap(', '.join(f'{a} {b}' for a, b in mk.items()), 62)))
ax = axs[11]; x = np.arange(len(Wd)); w = .27
for o, (c, lab, col) in enumerate([('r_all', 'all pairs of the type', '#999999'), ('r_same_lineage', 'same lineage make-up', '#d6604d'), ('r_same_hg', 'same haplogroup make-up', '#1b9e77')]):
    ax.bar(x + (o - 1) * w, Wd[c].fillna(0), w, color=col, label=lab)
ax.set_xticks(x); ax.set_xticklabels([f'"{t}"' for t in Wd.arrangement], rotation=30, ha='right', fontsize=8)
ax.axhline(0, c='k', lw=.6); ax.set_ylabel('Spearman, structure vs sequence'); ax.legend(fontsize=7)
ax.set_title('l  within-type correlation, by pairs kept', loc='left', fontsize=10)
fig.tight_layout(); fig.savefig(f'{F}/hap_scatter_lineage.png', dpi=130)
fig.savefig(f'{F}/hap_scatter_lineage.pdf', metadata={'CreationDate': None})

# ---- supplement figure: the haplogroup panels only (same-type pairs, each type, within-type correlation)
fig2, ax2 = plt.subplots(3, 3, figsize=(16, 14)); ax2 = ax2.ravel()
by_hg(ax2[0], ST, f'a  same-type pairs ({int(ST.sum())}), by haplogroup make-up')
for j, t in enumerate(types10):
    k = T == t; mk = pd.Series(G[k]).value_counts()
    by_hg(ax2[1 + j], (T[ia] == t) & (T[ib] == t),
          f'{"bcdefgh"[j]}  "{t}": {int(k.sum())} haplotypes\n' + '\n'.join(textwrap.wrap(', '.join(f'{a} {b}' for a, b in mk.items()), 62)))
ax = ax2[8]
for o, (c, lab, col) in enumerate([('r_all', 'all pairs of the type', '#999999'), ('r_same_lineage', 'same lineage make-up', '#d6604d'), ('r_same_hg', 'same haplogroup make-up', '#1b9e77')]):
    ax.bar(x + (o - 1) * w, Wd[c].fillna(0), w, color=col, label=lab)
ax.set_xticks(x); ax.set_xticklabels([f'"{t}"' for t in Wd.arrangement], rotation=30, ha='right', fontsize=8)
ax.axhline(0, c='k', lw=.6); ax.set_ylabel('Spearman, structure vs sequence'); ax.legend(fontsize=7)
ax.set_title('i  within-type correlation, by pairs kept', loc='left', fontsize=10)
assert len(types10) == 7
fig2.tight_layout(); fig2.savefig(f'{F}/hap_scatter_hg.png', dpi=130)
fig2.savefig(f'{F}/hap_scatter_hg.pdf', metadata={'CreationDate': None})
