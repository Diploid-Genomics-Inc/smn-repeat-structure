# Pipeline step 09 (ancestry_structure); ported from explore/structure_vs_sequence/05_ancestry_structure.py (sha256 80534283bb66, commit 5e2b3bb),
# paths rewritten to results/sv/ and figures/; run from the project root after the steps before it.
"""Exploration (not a pipeline step; nothing here feeds the note).
Does the large-scale structure of the SMN region differ between continental groups?
Unit: resolved cohort haplotypes in AFR, AMR, EAS, EUR, SAS (ASL and unlabelled left out; groups as in 02_compare.py).
Tests (all by permuting the group label of whole PEOPLE, both haplotypes together, 9999 permutations, seed 1,
so the two haplotypes of a person are not treated as independent):
  * arrangement type (types on <5 such haplotypes pooled as 'rare'): chi-square statistic, Cramer's V;
  * per-haplotype counts (SMN copies, SMN1 and SMN2 by exon 7, SERF1, full NAIP, GTF2H2, large repeat units):
    Kruskal-Wallis H;
  * structural distance (pangea bundle dist, out/smn_struct.dist): PERMANOVA pseudo-F, all five groups and the
    four groups other than AFR.
"""
import numpy as np, pandas as pd, json, collections
from scipy.stats import kruskal, chi2_contingency
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
O = 'results/sv'; F = 'figures'
AFR = {'ACB', 'ASW', 'ESN', 'GWD', 'LWK', 'MSL', 'YRI', 'MKK'}
SUPER = {**{p: 'AFR' for p in AFR}, **{p: 'AMR' for p in ['PUR', 'PEL', 'CLM', 'MXL']},
         **{p: 'EAS' for p in ['JPT', 'KHV', 'CHS', 'CDX', 'CHB']}, **{p: 'EUR' for p in ['GBR', 'FIN', 'TSI', 'IBS']},
         **{p: 'SAS' for p in ['BEB', 'PJL', 'GIH', 'STU', 'ITU']}}
COL = {'AFR': '#8c510a', 'AMR': '#d73027', 'EAS': '#4575b4', 'EUR': '#1a9850', 'SAS': '#984ea3'}
G = list(COL); NPERM = 9999; rng = np.random.default_rng(1)

hc = pd.read_csv('tables/haplotype_copy_counts.tsv', sep='\t')
hc = hc[hc.resolved].copy(); hc['grp'] = hc.population.map(SUPER); hc = hc[hc.grp.notna()]
nh = pd.read_csv('tables/neighbour_haplotypes.tsv', sep='\t').set_index('src')
hu = pd.read_csv('results/decomp/haplotype_units.tsv', sep='\t').set_index('src')
ru = pd.read_csv('results/decomp/smn.repeat_units.tsv', sep='\t')
H = hc.set_index('src')
for c in ['n_smn', 'n_serf1', 'n_naip_full', 'n_gtf2h2']: H[c] = nh.loc[H.index, c]
H['n_units'] = hu.loc[H.index, 'units'].fillna('').str.split().str.len()
H = H.reset_index()
assert set(H.arrangement) == set(hu.loc[H.src, 'arrangement'])
samples = H['sample'].unique(); sgrp = H.groupby('sample').grp.first()
assert (H.groupby('sample').grp.nunique() == 1).all()
hs = H['sample'].to_numpy()
def perms():
    lab = sgrp.loc[samples].to_numpy()
    yield dict(zip(samples, lab))
    for _ in range(NPERM): yield dict(zip(samples, rng.permutation(lab)))
PERM = [np.array([m[s] for s in hs]) for m in perms()]
summary = {'haplotypes': int(len(H)), 'people': int(len(samples)), 'permutations': NPERM,
           'haplotypes_by_group': {g: int((H.grp == g).sum()) for g in G}}

# arrangement types
cnt = H.arrangement.value_counts(); keep = set(cnt[cnt >= 5].index)
H['atype'] = np.where(H.arrangement.isin(keep), H.arrangement, 'rare')
at = H.atype.to_numpy()
def chi(lab):
    t = pd.crosstab(at, lab); return chi2_contingency(t, correction=False)[0]
obs = chi(PERM[0]); null = np.array([chi(l) for l in PERM[1:]])
n = len(H); tab = pd.crosstab(H.atype, H.grp)[G]
summary['arrangement'] = {'types_tested': int(len(keep)), 'rare_pooled_haplotypes': int((at == 'rare').sum()),
    'chi2': float(obs), 'cramers_v': float(np.sqrt(obs / (n * (min(tab.shape) - 1)))),
    'p_perm': float((1 + (null >= obs).sum()) / (NPERM + 1))}
(tab / tab.sum(axis=0)).round(3).assign(haplotypes=tab.sum(axis=1)).to_csv(f'{O}/ancestry_arrangement_shares.tsv', sep='\t')

# counts
meas = {'n_smn': 'SMN copies', 'SMN1_ex7': 'SMN1 (exon 7)', 'SMN2_ex7': 'SMN2 (exon 7)', 'n_serf1': 'SERF1',
        'n_naip_full': 'full-length NAIP', 'n_gtf2h2': 'GTF2H2', 'n_units': 'large repeat units'}
rows = []
for c, lab in meas.items():
    x = H[c].to_numpy(float)
    def kw(l): return kruskal(*[x[l == g] for g in G]).statistic if len(set(x)) > 1 else 0.0
    o = kw(PERM[0]); nl = np.array([kw(l) for l in PERM[1:2000]])
    r = {'measure': lab, 'column': c, 'H': round(o, 2), 'p_perm_1999': (1 + (nl >= o).sum()) / 2000}
    for g in G: r[f'mean_{g}'] = round(x[PERM[0] == g].mean(), 2)
    rows.append(r)
CT = pd.DataFrame(rows); CT.to_csv(f'{O}/ancestry_counts.tsv', sep='\t', index=False)

# structural distance
st = pd.read_csv(f'{O}/smn_struct.dist', sep=' ', header=None, usecols=[0, 1, 2], names=['a', 'b', 'd'])
src = H.src.tolist(); ri = {r: i for i, r in enumerate(src)}
S = np.full((len(src),) * 2, np.nan)
for a, b, d in st[st.a.isin(ri) & st.b.isin(ri)].itertuples(index=False): S[ri[a], ri[b]] = d
assert not np.isnan(S).any(); S = (S + S.T) / 2; np.fill_diagonal(S, 0)
D2 = S ** 2
def pseudoF(lab, idx):
    d2 = D2[np.ix_(idx, idx)]; l = lab[idx]; N = len(idx); gs = sorted(set(l))
    sst = d2.sum() / (2 * N); ssw = sum(d2[np.ix_(l == g, l == g)].sum() / (2 * (l == g).sum()) for g in gs)
    return ((sst - ssw) / (len(gs) - 1)) / (ssw / (N - len(gs))), (sst - ssw) / sst
for name, idx in [('all5', np.arange(len(src))), ('non_AFR4', np.where(PERM[0] != 'AFR')[0])]:
    # for the non-AFR test, permute among non-AFR people only
    if name == 'non_AFR4':
        ss = [s for s in samples if sgrp[s] != 'AFR']; lab0 = sgrp.loc[ss].to_numpy()
        PL = [np.array([dict(zip(ss, lab0)).get(s, 'AFR') for s in hs])] + \
             [np.array([dict(zip(ss, rng.permutation(lab0))).get(s, 'AFR') for s in hs]) for _ in range(1999)]
    else: PL = PERM[:2000]
    o, r2 = pseudoF(PL[0], idx); nl = np.array([pseudoF(l, idx)[0] for l in PL[1:]])
    summary[f'permanova_{name}'] = {'pseudoF': float(o), 'R2': float(r2), 'p_perm_1999': float((1 + (nl >= o).sum()) / 2000)}
M = pd.DataFrame([[np.median(S[np.ix_(PERM[0] == a, PERM[0] == b)][np.triu_indices((PERM[0] == a).sum(), 1)] if a == b
                   else S[np.ix_(PERM[0] == a, PERM[0] == b)]) for b in G] for a in G], index=G, columns=G).round(3)
M.to_csv(f'{O}/ancestry_struct_median_dist.tsv', sep='\t')
summary['median_within'] = {g: float(M.loc[g, g]) for g in G}
json.dump(summary, open(f'{O}/ancestry_summary.json', 'w'), indent=1)

# figure
n_ = len(S); J = np.eye(n_) - 1 / n_; B = -0.5 * J @ D2 @ J; w, v = np.linalg.eigh(B); o_ = np.argsort(w)[::-1]
X = v[:, o_[:2]] * np.sqrt(np.maximum(w[o_[:2]], 0)); ve = w[o_[:2]] / w[w > 0].sum()
fig = plt.figure(figsize=(15, 5.2))
ax = fig.add_axes([0.04, 0.12, 0.30, 0.78])
order = [t for t in tab.sum(axis=1).sort_values(ascending=False).index if t != 'rare'] + ['rare']
sh = (tab / tab.sum(axis=0)).loc[order]; cm = plt.get_cmap('tab20'); bot = np.zeros(len(G))
for k, t in enumerate(order):
    ax.bar(G, sh.loc[t], bottom=bot, color=('#dddddd' if t == 'rare' else cm(k % 20)), edgecolor='w', lw=0.5, label=f'{t} ({tab.loc[t].sum()})'); bot += sh.loc[t].to_numpy()
ax.set_xticks(range(len(G))); ax.set_xticklabels([f'{g}\n({tab[g].sum()})' for g in G])
ax.set_ylabel('share of resolved haplotypes'); ax.set_title('a  arrangement type by group', loc='left')
ax.legend(fontsize=7, bbox_to_anchor=(1.0, 1.0), loc='upper left', title='type (haplotypes)', title_fontsize=7, frameon=False)
ax = fig.add_axes([0.46, 0.12, 0.20, 0.78])
im = ax.imshow(M.to_numpy(), cmap='viridis'); ax.set_xticks(range(5)); ax.set_yticks(range(5)); ax.set_xticklabels(G); ax.set_yticklabels(G)
for i in range(5):
    for j in range(5): ax.text(j, i, f'{M.iat[i, j]:.2f}', ha='center', va='center', fontsize=8, color='w' if M.iat[i, j] < M.to_numpy().mean() else 'k')
ax.set_title('b  median structural distance', loc='left'); fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
ax = fig.add_axes([0.78, 0.12, 0.20, 0.78])
for g in G[::-1]:
    m = PERM[0] == g; ax.scatter(X[m, 0], X[m, 1], s=9, color=COL[g], alpha=0.7, label=f'{g} ({m.sum()})', lw=0)
ax.set_xlabel(f'axis 1 ({ve[0]*100:.0f}%)'); ax.set_ylabel(f'axis 2 ({ve[1]*100:.0f}%)')
ax.set_title('c  structural distance, PCoA', loc='left'); ax.legend(fontsize=7, frameon=False)
fig.savefig(f'{F}/ancestry_structure.png', dpi=200); fig.savefig(f'{F}/ancestry_structure.pdf', metadata={'CreationDate': None})
print(json.dumps(summary, indent=1)); print(CT.to_string()); print(M); print((tab).to_string())
