# Pipeline step 19 (cluster_correspondence); ported from explore/structure_vs_sequence/16_cluster_correspondence.py (sha256 6226e016428a, commit 5e2b3bb),
# paths rewritten to results/sv/ and figures/; run from the project root after the steps before it.
"""Step 16: which clusters of the structure tree correspond to which clusters of the SMN1 sequence tree (Fig. 2)?
Trees rebuilt exactly as steps 02/14 (UPGMA + optimal leaf ordering, same 436 haplotypes, out/pair_distances.tsv.gz).
Each tree is cut into k clusters (k = 8, 12, 16, 24); every structure-cluster x sequence-cluster pair is tested for
excess overlap (hypergeometric, Bonferroni over the k*k pairs). Pairs are annotated by arrangement type, ancestry and
SMN1 lineage make-up (out/lineage_haplotypes.tsv), and by where they sit in each tree (leaf position, 0 = bottom,
1 = top of the drawn tree).
Writes out/cluster_pairs.tsv, out/cluster_pairs_members.tsv, fig/tanglegram_clusters.png/.pdf."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import average, optimal_leaf_ordering, dendrogram, fcluster
from scipy.spatial.distance import squareform
from scipy.stats import hypergeom
O, F = 'results/sv', 'figures'
P = pd.read_csv(f'{O}/pair_distances.tsv.gz', sep='\t')
rows = sorted(set(P.a) | set(P.b)); ri = {r: i for i, r in enumerate(rows)}; n0 = len(rows)
def mat(col):
    M = np.full((n0, n0), np.nan); np.fill_diagonal(M, 0)
    ia, ib = P.a.map(ri).to_numpy(), P.b.map(ri).to_numpy(); M[ia, ib] = M[ib, ia] = P[col].to_numpy(); return M
S, G1 = mat('struct'), mat('seq_SMN1')
k1 = np.where((~np.isnan(G1)).sum(1) > 1)[0]
S1, G1k = S[np.ix_(k1, k1)], G1[np.ix_(k1, k1)]
rk = [rows[i] for i in k1]; n = len(rk)
Zs = optimal_leaf_ordering(average(squareform(S1)), squareform(S1)); Zg = optimal_leaf_ordering(average(squareform(G1k)), squareform(G1k))
leafS = dendrogram(Zs, no_plot=True)['leaves']; leafG = dendrogram(Zg, no_plot=True)['leaves']
posS = np.empty(n); posS[leafS] = np.arange(n) / (n - 1); posG = np.empty(n); posG[leafG] = np.arange(n) / (n - 1)
L = pd.read_csv(f'{O}/lineage_haplotypes.tsv', sep='\t').set_index('src')
arr = np.array([L.loc[r, 'arrangement'] if r in L.index else 'NA' for r in rk], dtype=object)
sp = np.array([L.loc[r, 'superpop'] if r in L.index else 'reference' for r in rk], dtype=object)
mk = np.array([L.loc[r, 'composition'] if r in L.index else 'NA' for r in rk], dtype=object)
def top(v, m=3):
    s = pd.Series(v).value_counts(); return '; '.join(f'{a} {b}' for a, b in s.head(m).items())
out, mem = [], []
for k in (8, 12, 16, 24):
    cs, cg = fcluster(Zs, k, 'maxclust'), fcluster(Zg, k, 'maxclust')
    npairs = len(set(cs)) * len(set(cg))
    for i in sorted(set(cs)):
        for j in sorted(set(cg)):
            m = (cs == i) & (cg == j); o = int(m.sum())
            if o < 5: continue
            a, b = int((cs == i).sum()), int((cg == j).sum())
            p = float(hypergeom.sf(o - 1, n, a, b))
            out.append({'k': k, 'struct_cluster': i, 'seq_cluster': j, 'struct_size': a, 'seq_size': b, 'overlap': o,
                        'expected': round(a * b / n, 1), 'frac_of_struct': round(o / a, 3), 'frac_of_seq': round(o / b, 3),
                        'jaccard': round(o / (a + b - o), 3), 'p': p, 'p_bonf': min(1.0, p * npairs),
                        'struct_pos_median': round(float(np.median(posS[m])), 3), 'seq_pos_median': round(float(np.median(posG[m])), 3),
                        'arrangements': top(arr[m]), 'superpop': top(sp[m], 5), 'SMN1_lineages': top(mk[m])})
            for x in np.where(m)[0]: mem.append({'k': k, 'struct_cluster': i, 'seq_cluster': j, 'src': rk[x], 'arrangement': arr[x], 'superpop': sp[x], 'SMN1_lineages': mk[x]})
T = pd.DataFrame(out).sort_values(['k', 'p']); T.to_csv(f'{O}/cluster_pairs.tsv', sep='\t', index=False)
pd.DataFrame(mem).to_csv(f'{O}/cluster_pairs_members.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 60)
for k in (8, 16):
    print(f'--- k={k}'); print(T[(T.k == k) & (T.p_bonf < 1e-3)].drop(columns=['p']).to_string(index=False))
# figure: k = 16, the strongest significant pairs with >= 10 haplotypes. Since the researcher's decision of 2026-11-14 the
# main-text bands come from step 19b (ML clades); this UPGMA-on-both-sides figure is the supplement's check, and it also
# carries the strips the earlier supplementary tanglegrams drew (ancestry, arrangement, SMN1 lineage make-up).
k = 16; cs, cg = fcluster(Zs, k, 'maxclust'), fcluster(Zg, k, 'maxclust')
H = T[(T.k == k) & (T.p_bonf < 1e-3) & (T.overlap >= 10)].sort_values('overlap', ascending=False).head(6).reset_index(drop=True)
pal = ['#d73027', '#4575b4', '#1a9850', '#984ea3', '#ff7f00', '#a65628']
SUPER_COL = {'AFR': '#8c510a', 'AMR': '#d73027', 'EAS': '#4575b4', 'EUR': '#1a9850', 'SAS': '#984ea3', 'other': '#bdbdbd', 'reference': '#000000'}
import collections
TOPARR = [a for a, _ in collections.Counter(arr).most_common(7) if a != 'NA'][:6]
ARRCOL = dict(zip(TOPARR, ['#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02']))
LIN_COL = {'A': '#9ecae1', 'AA': '#2171b5', 'AB': '#8856a7', 'B': '#fc9272', 'BB': '#cb181d'}
LIN_LAB = {'A': 'one A', 'AA': 'two A', 'AB': 'one A + one B', 'B': 'one B', 'BB': 'two B'}
plt.rcParams['lines.linewidth'] = 0.5
fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(11, 13.5), gridspec_kw={'width_ratios': [1, 1.2, 1]})
dendrogram(Zs, ax=a1, orientation='left', no_labels=True, color_threshold=0, above_threshold_color='k')
dendrogram(Zg, ax=a3, orientation='right', no_labels=True, color_threshold=0, above_threshold_color='k')
yS, yG = posS * (n - 1), posG * (n - 1)
for x in range(n): a2.plot([0.13, 0.87], [yS[x], yG[x]], lw=.3, color=(.85, .85, .85))
hand = []
for h, r in H.iterrows():
    m = (cs == r.struct_cluster) & (cg == r.seq_cluster)
    for x in np.where(m)[0]: a2.plot([0.13, 0.87], [yS[x], yG[x]], lw=.8, color=pal[h])
    hand.append(plt.Line2D([], [], color=pal[h], lw=2, label=f'{r.overlap} haplotypes ({r.expected} expected): {r.arrangements.split(";")[0]}'))
hand.append(plt.Line2D([], [], color=(.75, .75, .75), lw=2, label='in no such pair'))
for x in range(n):
    a2.add_patch(plt.Rectangle((0, yS[x] - .5), .055, 1, color=SUPER_COL.get(sp[x], '#bdbdbd'), lw=0))
    a2.add_patch(plt.Rectangle((.065, yS[x] - .5), .055, 1, color=ARRCOL.get(arr[x], '#e0e0e0'), lw=0))
    a2.add_patch(plt.Rectangle((.88, yG[x] - .5), .055, 1, color=LIN_COL.get(mk[x], '#e0e0e0'), lw=0))
    a2.add_patch(plt.Rectangle((.945, yG[x] - .5), .055, 1, color=SUPER_COL.get(sp[x], '#bdbdbd'), lw=0))
for xx, t_ in [(.0275, 'ancestry'), (.0925, 'arrangement'), (.9075, 'SMN1 lineage'), (.9725, 'ancestry')]:
    a2.text(xx, n - .5 + n * .004, t_, rotation=90, ha='center', va='bottom', fontsize=7)
a2.set_xlim(0, 1); a2.set_ylim(-0.5, n - 0.5); a2.axis('off')
for a in (a1, a3): a.set_ylim(0, 10 * n); a.set_yticks([]); [a.spines[s_].set_visible(False) for s_ in ('top', 'right', 'left')]
a1.set_title('large-scale structure\n(PANGEA bundle distance, UPGMA)', fontsize=10); a3.set_title('SMN1 copy sequences\n(SMN1 copy-set distance, UPGMA)', fontsize=10)
from matplotlib.ticker import FuncFormatter, MaxNLocator
a3.xaxis.set_major_locator(MaxNLocator(4)); a3.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v*1e3:.2g}')); a3.set_xlabel('differences per kb', fontsize=8); a1.set_xlabel('bundle distance', fontsize=8)
a2.legend(handles=hand, loc='upper center', bbox_to_anchor=(.5, -0.005), fontsize=7.5, frameon=False, ncol=1,
          title='line: same cluster in both trees (k = 16; top arrangement type)', title_fontsize=8)
nsp = collections.Counter(sp); narr = collections.Counter(arr); nlin = collections.Counter(mk)
h2 = [plt.Rectangle((0, 0), 1, 1, color=SUPER_COL[v], label=f'{"ASL or no label" if v == "other" else "GRCh38/CHM13" if v == "reference" else v} ({nsp[v]})') for v in SUPER_COL if nsp.get(v)]
h2 += [plt.Rectangle((0, 0), 1, 1, color=ARRCOL[a], label=f'"{a}" ({narr[a]})') for a in TOPARR] + [plt.Rectangle((0, 0), 1, 1, color='#e0e0e0', label='other types')]
a1.legend(handles=h2, loc='upper center', bbox_to_anchor=(.5, -0.03), fontsize=7, frameon=False, ncol=2, title='strips: ancestry, arrangement (haplotypes)', title_fontsize=8)
h3 = [plt.Rectangle((0, 0), 1, 1, color=LIN_COL[v], label=f'{LIN_LAB[v]} ({nlin[v]})') for v in LIN_COL if nlin.get(v)]
a3.legend(handles=h3, loc='upper center', bbox_to_anchor=(.5, -0.03), fontsize=7, frameon=False, ncol=2, title='strip: SMN1 copies by lineage', title_fontsize=8)
fig.tight_layout(); fig.savefig(f'{F}/tanglegram_clusters.png', dpi=200); fig.savefig(f'{F}/tanglegram_clusters.pdf', metadata={'CreationDate': None}); plt.close(fig)
print(H[['struct_cluster', 'seq_cluster', 'overlap', 'expected', 'struct_pos_median', 'seq_pos_median', 'arrangements']].to_string())
