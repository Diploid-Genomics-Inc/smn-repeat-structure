# Pipeline step 07 (compare); ported from explore/structure_vs_sequence/02_compare.py (sha256 5fd99eda99e1, commit 5e2b3bb),
# paths rewritten to results/sv/ and figures/; run from the project root after the steps before it.
"""Exploration (not a pipeline step; nothing here feeds the note).
Compare the large-scale structure of each resolved haplotype with the sequence of its SMN copies.

Structure: pangea cohort bundle dist on results/decomp/smn.bed (principal-bundle composition and order of each
whole window, top 50 bundles), out/smn_struct.dist; 438 resolved haplotypes + GRCh38 + CHM13.
Sequence: p-distance between SMN copies (gene body only; flanks and del7-8 copies dropped, see below) from out/copies.aln.fa (01_gene_tree.sh), ignoring
columns where either copy has a gap or N. Haplotype sequence distance = mean best match between the two copy
sets (each copy to its nearest copy on the other haplotype, both ways, averaged).
Comparisons: Mantel test (Spearman, 999 permutations, seed 1); UPGMA trees of both matrices drawn as a
tanglegram; copy-level UPGMA tree with structural labels; nearest-sequence-neighbor agreement of each
structural label, within SMN1 and within SMN2, against label permutations (999, seed 1).
"""
import numpy as np, pandas as pd, collections, json, os
from scipy.cluster.hierarchy import average, dendrogram, cophenet
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
O = 'results/sv'; F = 'figures'
AFR = {'ACB', 'ASW', 'ESN', 'GWD', 'LWK', 'MSL', 'YRI', 'MKK'}
# continental groups for the plot only (1000 Genomes superpopulations); AFR is the study's list above,
# so ASL (HG06807) stays outside it, as in the note, and is drawn with the unlabelled samples
SUPER = {**{p: 'AFR' for p in AFR}, **{p: 'AMR' for p in ['PUR', 'PEL', 'CLM', 'MXL']},
         **{p: 'EAS' for p in ['JPT', 'KHV', 'CHS', 'CDX', 'CHB']}, **{p: 'EUR' for p in ['GBR', 'FIN', 'TSI', 'IBS']},
         **{p: 'SAS' for p in ['BEB', 'PJL', 'GIH', 'STU', 'ITU']}}
SUPER_COL = {'AFR': '#8c510a', 'AMR': '#d73027', 'EAS': '#4575b4', 'EUR': '#1a9850', 'SAS': '#984ea3', 'other': '#bdbdbd'}
rng = np.random.default_rng(1)
summary = {}

# ---- copy p-distances
names, seqs = [], []
for l in open(f'{O}/copies.aln.fa'):
    l = l.strip()
    if l.startswith('>'): names.append(l[1:]); seqs.append([])
    else: seqs[-1].append(l.upper())
A = np.array([np.frombuffer(''.join(s).encode(), dtype=np.uint8) for s in seqs])
# Checked in exploration: (1) mafft misplaces the 6.3-kb deletion of the 8 del7-8 copies (~508 mismatches to every
# copy instead of a gap), so they are left out; (2) the 2 kb flanks carry flank-specific differences (HG02965 hap2's
# SMN1 copy: 72 of 101 differences in its last 1 kb), so only columns of the gene body are used, located by the
# GRCh38 SMN1 copy (ungapped positions 2000 .. length-2000).
_t = pd.read_csv(f'{O}/copy_tips.tsv', sep='\t').set_index('tip').loc[names]
keep_rows = (_t.copy_class != 'del7_8_SMN2like').to_numpy()
ref = np.where((_t.copy_class == 'SMN1').to_numpy() & _t.src.str.startswith('GCA_000001405').to_numpy())[0][0]
rpos = np.cumsum(A[ref] != ord('-')) - 1; rlen = rpos[-1] + 1
colkeep = (A[ref] != ord('-')) & (rpos >= 2000) & (rpos < rlen - 2000)
# keep insertion columns inside the gene body too
inside = (rpos >= 2000) & (rpos < rlen - 2000); colkeep = inside | ((A[ref] == ord('-')) & (rpos >= 1999) & (rpos < rlen - 2000))
A = A[keep_rows][:, colkeep]; names = [n for n, k in zip(names, keep_rows) if k]
summary['del7_8_copies_left_out'] = int((~keep_rows).sum()); summary['gene_body_columns'] = int(colkeep.sum())
valid = np.isin(A, list(b'ACGT'))
V = valid.astype(np.float32); nvalid = V @ V.T
var = np.array([len(set(A[valid[:, j], j])) > 1 for j in range(A.shape[1])])
same = np.zeros_like(nvalid)
for b in b'ACGT':
    M = (A[:, var] == b).astype(np.float32); same += M @ M.T
nvar_valid = V[:, var] @ V[:, var].T
D = (nvar_valid - same) / nvalid
np.fill_diagonal(D, 0)
summary['copies'] = len(names); summary['aln_columns_used'] = int(A.shape[1]); summary['variable_columns'] = int(var.sum())

tips = pd.read_csv(f'{O}/copy_tips.tsv', sep='\t').set_index('tip').loc[names].reset_index()
nb = pd.read_csv('tables/neighbour_copies.tsv', sep='\t').set_index('copy_id')
hc = pd.read_csv('tables/haplotype_copy_counts.tsv', sep='\t')
res = set(hc[hc.resolved].src) | {s for s in nb.index.str.split('|').str[0] if 'PanSN' in s}
pop = dict(zip(hc.src, hc.population))
hu = pd.read_csv('results/decomp/haplotype_units.tsv', sep='\t'); arr = dict(zip(hu.src, hu.arrangement))

# ---- haplotype level
st = pd.read_csv(f'{O}/smn_struct.dist', sep=' ', header=None, usecols=[0, 1, 2], names=['a', 'b', 'd'])
rows = sorted(set(st.a) & res & set(tips.src))
ri = {r: i for i, r in enumerate(rows)}
S = np.full((len(rows), len(rows)), np.nan)
for a, b, d in st[st.a.isin(ri) & st.b.isin(ri)].itertuples(index=False): S[ri[a], ri[b]] = d
assert not np.isnan(S).any(); S = (S + S.T) / 2; np.fill_diagonal(S, 0)
tsrc = tips.src.to_numpy(dtype=object); tcls = tips.copy_class.to_numpy(dtype=object)
def setdist(classes):
    ix = {r: np.where((tsrc == r) & np.isin(tcls, classes))[0] for r in rows}
    M = np.full((len(rows), len(rows)), np.nan); np.fill_diagonal(M, 0)
    for i, a in enumerate(rows):
        if not len(ix[a]): continue
        for j in range(i + 1, len(rows)):
            if not len(ix[rows[j]]): continue
            sub = D[np.ix_(ix[a], ix[rows[j]])]
            M[i, j] = M[j, i] = (sub.min(1).mean() + sub.min(0).mean()) / 2
    return M
G = setdist(['SMN1', 'SMN2', 'hyb_1_2', 'hyb_2_1'])   # all copies: mixes copy content with lineage
G1 = setdist(['SMN1'])                                                       # SMN1 copies only: lineage
G2 = setdist(['SMN2'])                                                       # SMN2 copies only: lineage
summary['haplotypes'] = len(rows)
iu = np.triu_indices(len(rows), 1)
r0 = spearmanr(S[iu], G[iu]).statistic
perm = []
for _ in range(999):
    p = rng.permutation(len(rows)); perm.append(spearmanr(S[iu], G[np.ix_(p, p)][iu]).statistic)
summary['mantel_spearman'] = round(float(r0), 4)
summary['mantel_p'] = (1 + sum(x >= r0 for x in perm)) / 1000
afrrow = np.array([pop.get(r) in AFR for r in rows])
def mantel(M, keep, nperm=999, strata=None):
    k = np.where(keep & ((~np.isnan(M)).sum(1) > 1))[0]
    Sk, Mk = S[np.ix_(k, k)], M[np.ix_(k, k)]; u = np.triu_indices(len(k), 1)
    obs = spearmanr(Sk[u], Mk[u]).statistic; null = []
    for _ in range(nperm):
        if strata is None: p = rng.permutation(len(k))
        else:
            p = np.arange(len(k)); st_ = strata[k]
            for g in (True, False):
                w = np.where(st_ == g)[0]; p[w] = w[rng.permutation(len(w))]
        null.append(spearmanr(Sk[u], Mk[np.ix_(p, p)][u]).statistic)
    return {'haplotypes': int(len(k)), 'spearman': round(float(obs), 4), 'p': (1 + sum(x >= obs for x in null)) / (nperm + 1)}
allrows = np.ones(len(rows), bool)
summary['mantel_SMN1_only'] = mantel(G1, allrows)
summary['mantel_SMN1_only_perm_within_ancestry'] = mantel(G1, allrows, strata=afrrow)
summary['mantel_SMN1_only_AFR'] = mantel(G1, afrrow)
summary['mantel_SMN1_only_nonAFR'] = mantel(G1, ~afrrow)
summary['mantel_SMN2_only'] = mantel(G2, allrows)
summary['mantel_SMN2_only_perm_within_ancestry'] = mantel(G2, allrows, strata=afrrow)
lab = np.array([arr.get(r, '?') for r in rows])
sameA = (lab[:, None] == lab[None, :])[iu]
summary['struct_dist_same_vs_diff_arrangement'] = [round(float(np.median(S[iu][sameA])), 4), round(float(np.median(S[iu][~sameA])), 4)]
summary['seq_dist_same_vs_diff_arrangement'] = [round(float(np.median(G[iu][sameA])), 6), round(float(np.median(G[iu][~sameA])), 6)]
from scipy.cluster.hierarchy import optimal_leaf_ordering
k1 = np.where((~np.isnan(G1)).sum(1) > 1)[0]   # haplotypes with an SMN1 copy
summary['tanglegram_haplotypes'] = int(len(k1))
S1, G1k = S[np.ix_(k1, k1)], G1[np.ix_(k1, k1)]
Zs = optimal_leaf_ordering(average(squareform(S1)), squareform(S1)); Zg = optimal_leaf_ordering(average(squareform(G1k)), squareform(G1k))
summary['cophenetic_corr_struct_tree_vs_seq_tree'] = round(float(np.corrcoef(cophenet(Zs), cophenet(Zg))[0, 1]), 4)
pd.DataFrame({'a': [rows[i] for i in iu[0]], 'b': [rows[j] for j in iu[1]], 'struct': S[iu], 'seq_all_copies': G[iu], 'seq_SMN1': G1[iu], 'seq_SMN2': G2[iu], 'same_arrangement': sameA}).to_csv(f'{O}/pair_distances.tsv.gz', sep='\t', index=False, compression={'method': 'gzip', 'mtime': 0})

top = [a for a, _ in collections.Counter(lab).most_common(7)]
cols = dict(zip(top, ['#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02', '#a6761d'])); grey = (0.75, 0.75, 0.75)

# fig: scatter
fig, axs = plt.subplots(1, 2, figsize=(12, 5))
for ax, M, ttl in [(axs[0], G, f'a  all copies of each haplotype (Spearman {r0:.2f})'), (axs[1], G1, f'b  SMN1 copies only (Spearman {summary["mantel_SMN1_only"]["spearman"]:.2f}, {summary["mantel_SMN1_only"]["haplotypes"]} haplotypes)')]:
    ok = ~np.isnan(M[iu])
    ax.scatter(S[iu][~sameA & ok], M[iu][~sameA & ok] * 1e3, s=2, c=[grey], alpha=.25, label='different arrangement type', rasterized=True)
    ax.scatter(S[iu][sameA & ok], M[iu][sameA & ok] * 1e3, s=2, c='#d95f02', alpha=.4, label='same arrangement type', rasterized=True)
    ax.set_xlabel('structural distance (PANGEA bundle dist, whole window)'); ax.set_ylabel('copy-set sequence distance (differences per kb)')
    ax.set_title(ttl, fontsize=9, loc='left')
axs[0].legend(markerscale=5, fontsize=8, frameon=False)
fig.tight_layout()
fig.savefig(f'{F}/hap_scatter.png', dpi=200); plt.close(fig)

# fig: tanglegram
ds = dendrogram(Zs, no_plot=True); dg = dendrogram(Zg, no_plot=True)
rk = [rows[i] for i in k1]
posS = {rk[k]: i for i, k in enumerate(ds['leaves'])}; posG = {rk[k]: i for i, k in enumerate(dg['leaves'])}
fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(11, 13), gridspec_kw={'width_ratios': [1, 1.1, 1]})
dendrogram(Zs, ax=a1, orientation='left', no_labels=True, color_threshold=0, above_threshold_color='k')
dendrogram(Zg, ax=a3, orientation='right', no_labels=True, color_threshold=0, above_threshold_color='k')
n = len(rk)
for r in rk:
    a2.plot([0.13, 0.87], [posS[r], posG[r]], lw=.5, color=cols.get(arr.get(r), grey), alpha=.8)
# ancestry strip beside the leaves of each tree (continental group; GRCh38 and CHM13 in black)
def sp_of(r):
    return SUPER.get(pop.get(r), 'other') if r in pop else 'reference'
SP_COL = {**SUPER_COL, 'reference': 'k'}
for r in rk:
    c = SP_COL[sp_of(r)]
    a2.add_patch(plt.Rectangle((0.0, posS[r] - .5), 0.05, 1, color=c, lw=0))
    a2.add_patch(plt.Rectangle((0.95, posG[r] - .5), 0.05, 1, color=c, lw=0))
a2.text(0.025, n - .5 + n * .004, 'ancestry', ha='center', va='bottom', fontsize=6.5, rotation=90)
# SMN1 lineage strip (written by step 06, out/lineage_haplotypes.tsv; rerun this step after 06 to draw it)
LIN_COL = {'A': '#9ecae1', 'AA': '#2171b5', 'AB': '#8856a7', 'B': '#fc9272', 'BB': '#cb181d'}
LIN_LAB = {'A': 'one A', 'AA': 'two A', 'AB': 'one A + one B', 'B': 'one B', 'BB': 'two B'}
linf = f'{O}/lineage_haplotypes.tsv'
comp = pd.read_csv(linf, sep='\t').set_index('src').composition.to_dict() if os.path.exists(linf) else {}
nlin = collections.Counter(comp.get(r) for r in rk)
if comp:
    for r in rk:
        c = LIN_COL.get(comp.get(r), '#d9d9d9')
        a2.add_patch(plt.Rectangle((0.06, posS[r] - .5), 0.05, 1, color=c, lw=0))
        a2.add_patch(plt.Rectangle((0.89, posG[r] - .5), 0.05, 1, color=c, lw=0))
    a2.text(0.085, n - .5 + n * .004, 'SMN1 lineage', ha='center', va='bottom', fontsize=6.5, rotation=90)
    a2.text(0.915, n - .5 + n * .004, 'SMN1 lineage', ha='center', va='bottom', fontsize=6.5, rotation=90)
a2.text(0.975, n - .5 + n * .004, 'ancestry', ha='center', va='bottom', fontsize=6.5, rotation=90)
a2.set_xlim(0, 1); a2.set_ylim(-.5, n - .5); a2.axis('off')
nsp_h = collections.Counter(sp_of(r) for r in rk)
summary_tangle_superpop = dict(nsp_h)
for a in (a1, a3): a.set_ylim(0, 10 * n); a.set_xticks([]); [s.set_visible(False) for s in a.spines.values()]
a1.set_title('large-scale structure\n(bundle composition, UPGMA)', fontsize=10); a3.set_title('SMN1 copy sequences\n(SMN1 copy-set distance, UPGMA)', fontsize=10)
h = [plt.Line2D([], [], color=cols[t], lw=2, label=f'{t} ({sum(arr.get(r) == t for r in rk)})') for t in top] + [plt.Line2D([], [], color=grey, lw=2, label='other types')]
lg = a2.legend(handles=h, loc='upper center', bbox_to_anchor=(.5, 0), fontsize=8, frameon=False, title='line: arrangement type (haplotypes)', ncol=2)
a2.add_artist(lg)
lab_sp = {'other': 'ASL or no label', 'reference': 'GRCh38 / CHM13'}
h2 = [plt.Rectangle((0, 0), 1, 1, color=SP_COL[v], label=f'{lab_sp.get(v, v)} ({nsp_h[v]})') for v in ['AFR', 'AMR', 'EAS', 'EUR', 'SAS', 'other', 'reference'] if nsp_h.get(v)]
a1.legend(handles=h2, loc='upper center', bbox_to_anchor=(.5, -0.005), fontsize=8, frameon=False, title='outer strip: ancestry group (haplotypes)', ncol=2)
if comp:
    h3 = [plt.Rectangle((0, 0), 1, 1, color=LIN_COL[v], label=f'{LIN_LAB[v]} ({nlin[v]})') for v in LIN_COL if nlin.get(v)]
    a3.legend(handles=h3, loc='upper center', bbox_to_anchor=(.5, -0.005), fontsize=8, frameon=False, title='inner strip: SMN1 copies by lineage', ncol=2)
    summary_tangle_lineage = {k: int(v) for k, v in nlin.items() if k}
fig.suptitle(f'Same {n} haplotypes (those with an SMN1 copy) in both trees; cophenetic correlation {summary["cophenetic_corr_struct_tree_vs_seq_tree"]:.2f}', fontsize=10)
fig.tight_layout(); fig.savefig(f'{F}/hap_tanglegram.png', dpi=200); plt.close(fig)

# ---- copy level
labels = pd.DataFrame({'tip': names, 'copy_id': tips.copy_id, 'src': tips.src, 'cls': tips.copy_class.replace({'hyb_1_2': 'h12', 'hyb_2_1': 'h21', 'del7_8_SMN2like': 'd78'})})
def position(cid):
    if cid not in nb.index: return None
    r = nb.loc[cid]
    if not r.resolved: return None
    if r.n_in_window == 1: return 'single'
    return 'centromeric' if r.rank_from_cen == 1 else ('telomeric' if r.rank_from_cen == r.n_in_window else 'middle')
labels['position'] = [position(c) for c in labels.copy_id]
labels['strand'] = [nb.strand_grch38.get(c) for c in labels.copy_id]
labels['naip3'] = [nb.naip_3p.get(c) for c in labels.copy_id]
labels['env'] = [nb.environment.get(c) for c in labels.copy_id]
labels['arrangement'] = [arr.get(s) if s in res else None for s in labels.src]
labels['afr'] = [pop.get(s) in AFR for s in labels.src]
labels['superpop'] = [SUPER.get(pop.get(s), 'other') if s in pop else None for s in labels.src]
import io; _b = io.StringIO(); labels.to_csv(_b, sep='\t', index=False)
assert _b.getvalue() == open(f'{O}/copy_labels.tsv').read(), 'copy labels differ from step 04'

nn = []
for cls in ['SMN1', 'SMN2']:
    sel = np.where((labels.cls == cls).values)[0]
    Dc = D[np.ix_(sel, sel)].copy(); srcs = labels.src.to_numpy(dtype=object)[sel]
    Dc[srcs[:, None] == srcs[None, :]] = np.inf   # neighbor on another haplotype
    neigh = [np.where(Dc[i] == Dc[i].min())[0] for i in range(len(sel))]
    for f in ['position', 'strand', 'naip3', 'env', 'arrangement', 'afr']:
        y = labels[f].to_numpy(dtype=object)[sel]
        ok = np.array([v is not None and v == v for v in y])
        def agree(yy):
            return np.mean([np.mean(yy[ne[ok[ne]]] == yy[i]) for i, ne in enumerate(neigh) if ok[i] and ok[ne].any()])
        obs = agree(y)
        null, nullA = [], []
        afr_ = labels.afr.to_numpy(dtype=object)[sel]
        for _ in range(999):
            yp = y.copy(); yp[ok] = rng.permutation(y[ok]); null.append(agree(yp))
            if f != 'afr':
                yp = y.copy()
                for g in (True, False):
                    w = np.where(ok & (afr_ == g))[0]; yp[w] = y[w][rng.permutation(len(w))]
                nullA.append(agree(yp))
        nn.append({'class': cls, 'label': f, 'copies': int(ok.sum()), 'agreement': round(obs, 3), 'null_mean': round(float(np.mean(null)), 3),
                   'p': (1 + sum(x >= obs for x in null)) / 1000,
                   'null_mean_within_ancestry': round(float(np.mean(nullA)), 3) if nullA else None,
                   'p_within_ancestry': (1 + sum(x >= obs for x in nullA)) / 1000 if nullA else None, 'levels': json.dumps(collections.Counter(map(str, y[ok])).most_common(6))})
nn = pd.DataFrame(nn); nn.to_csv(f'{O}/copy_neighbor_agreement.tsv', sep='\t', index=False)

# fig: copy tree with strips
Zc = average(squareform(D, checks=False)); dc = dendrogram(Zc, no_plot=True); order = dc['leaves']
strips = [('class', 'cls', {'SMN1': 'C0', 'SMN2': 'C1', 'h12': 'C2', 'h21': 'C3', 'd78': 'k'}),
          ('position', 'position', {'centromeric': 'C4', 'middle': 'C8', 'telomeric': 'C9', 'single': 'C7'}),
          ('strand', 'strand', {'+': 'C2', '-': 'C5'}),
          ('NAIP 3\'', 'naip3', {'full': 'C0', 'trunc5': 'C1'}),
          ('environment', 'env', {'SMN1': 'C0', 'SMN2': 'C1'}),
          ('arrangement', 'arrangement', cols),
          ('ancestry', 'superpop', SUPER_COL)]
fig = plt.figure(figsize=(9, 14)); gs = fig.add_gridspec(1, 1 + len(strips), width_ratios=[5] + [.45] * len(strips), wspace=.05)
ax = fig.add_subplot(gs[0]); dendrogram(Zc, ax=ax, orientation='left', no_labels=True, color_threshold=0, above_threshold_color='k', link_color_func=lambda k: 'k')
ax.set_xlabel('p-distance (UPGMA height)')
for s in ['top', 'right', 'left']: ax.spines[s].set_visible(False)
ax.set_yticks([])
for k, (title, col, cmap) in enumerate(strips):
    a = fig.add_subplot(gs[k + 1], sharey=ax)
    for y, i in enumerate(order):
        v = labels[col].iloc[i]
        c = cmap.get(v, (1, 1, 1) if v is None or v is False or v != v else grey)
        a.add_patch(plt.Rectangle((0, 10 * y), 1, 10, color=c, lw=0))
    a.set_xlim(0, 1); a.set_xticks([]); a.set_yticks([]); a.set_title(title, fontsize=8, rotation=60, ha='left', va='bottom')
    for s in a.spines.values(): s.set_visible(False)
h = []
for title, col, cmap in strips[:5]:
    h += [plt.Rectangle((0, 0), 1, 1, color=c, label=f'{title}: {v}') for v, c in cmap.items()]
h += [plt.Rectangle((0, 0), 1, 1, color=cols[t], label=f'arrangement: {t}') for t in top] + [plt.Rectangle((0, 0), 1, 1, color=grey, label='arrangement: other')]
nsp = collections.Counter(labels.superpop.dropna())
h += [plt.Rectangle((0, 0), 1, 1, color=c, label=(f'ancestry: {v}' if v != 'other' else 'ancestry: ASL or no label') + f' ({nsp.get(v, 0)} copies)') for v, c in SUPER_COL.items()]
fig.legend(handles=h, loc='lower center', ncol=4, fontsize=6.5, frameon=False)
fig.suptitle(f'{len(names)} SMN copies: sequence tree (UPGMA on p-distance) beside each copy\'s structural context\n(white = not assessed: copy on an unresolved haplotype)', fontsize=9, y=.985)
fig.subplots_adjust(bottom=.135, top=.89); fig.savefig(f'{F}/copy_tree_strips.png', dpi=200); plt.close(fig)

json.dump(summary, open(f'{O}/summary.json', 'w'), indent=1)
print(json.dumps(summary, indent=1)); print(nn.drop(columns='levels').to_string())
