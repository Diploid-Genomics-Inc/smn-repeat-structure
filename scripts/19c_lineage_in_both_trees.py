# Pipeline step 19c (lineage_in_both_trees; was step 18, renumbered to run after step 19b, whose ML clusters it reads); ported from explore/structure_vs_sequence/14_lineage_in_both_trees.py (sha256 2a22b24f63ef, commit 5e2b3bb),
# paths rewritten to results/sv/ and figures/; run from the project root after the steps before it.
"""Step 14: do haplotypes of the same SMN1 lineage make-up (A only, one A + one B, B only, two B) sit together
in BOTH trees of the tanglegram: the structure tree (pangea bundle distance) and the SMN1 sequence tree?
Reads out/pair_distances.tsv.gz (step 02) and out/lineage_haplotypes.tsv (step 06). Structure tree: UPGMA +
optimal leaf ordering on the bundle distance (as step 07). Sequence tree: the gene-body maximum-likelihood tree (step 04a),
read at the haplotype level by smnlib.ml_smn1_hap (patristic distance between SMN1 copy sets; leaf order = the tree's tip
order); its clusters are step 19b's ML clusters (clades of support >= 95 holding >= 20 SMN1 copies).
Writes results/sv/lineage_tree_match.tsv and .json."""
import json, sys, numpy as np, pandas as pd
sys.path.insert(0, 'scripts'); from smnlib import ml_smn1_hap
from scipy.cluster.hierarchy import average, optimal_leaf_ordering, dendrogram, cophenet
from scipy.spatial.distance import squareform
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
Zs = optimal_leaf_ordering(average(squareform(S1)), squareform(S1))
H, leafG = ml_smn1_hap(rk, O)   # sequence side: the ML tree
L = pd.read_csv(f'{O}/lineage_haplotypes.tsv', sep='\t').set_index('src')
def make_up(r):
    if r not in L.index: return 'none'
    a, b = int(L.loc[r, 'nA']), int(L.loc[r, 'nB'])
    return {(1, 0): 'A', (2, 0): 'AA', (1, 1): 'AB', (0, 1): 'B', (0, 2): 'BB'}.get((a, b), f'A{a}B{b}')
lab = np.array([make_up(r) for r in rk]); sp = np.array([L.loc[r, 'superpop'] if r in L.index else 'reference' for r in rk])
afr = sp == 'AFR'
rng = np.random.default_rng(14)
Cs, Cg = squareform(cophenet(Zs)), H
leafS = dendrogram(Zs, no_plot=True)['leaves']
def runs(order, mask):
    seq = mask[order]; return int(seq[0]) + int(((seq[1:]) & (~seq[:-1])).sum())
def nn_agree(M, mask):
    Mm = M.copy(); np.fill_diagonal(Mm, np.inf); nn = np.nanargmin(Mm, 1); idx = np.where(mask)[0]
    return float(mask[nn[idx]].mean())
def within(C, mask):
    i = np.where(mask)[0]; u = np.triu_indices(len(i), 1); return float(np.median(C[np.ix_(i, i)][u]))
out = []; NP = 1999
for g in ['A', 'AA', 'AB', 'B', 'BB']:
    m = lab == g
    if m.sum() < 2: continue
    for tree, C, M, order in [('structure', Cs, S1, leafS), ('SMN1 sequence', Cg, H, leafG)]:
        obs_w, obs_nn, obs_runs = within(C, m), nn_agree(M, m), runs(order, m)
        # null 1: shuffle make-up labels among all haplotypes; null 2: shuffle only within AFR / non-AFR
        nul = {'all': [], 'AFR-stratified': []}
        for kind in nul:
            for _ in range(NP):
                if kind == 'all': p = rng.permutation(n)
                else:
                    p = np.arange(n)
                    for s in (True, False):
                        w = np.where(afr == s)[0]; p[w] = w[rng.permutation(len(w))]
                mp = m[p]; nul[kind].append((within(C, mp), nn_agree(M, mp)))
        row = {'make_up': g, 'haplotypes': int(m.sum()), 'AFR': int((m & afr).sum()), 'tree': tree,
               'median_within_cophenetic': round(obs_w, 4), 'median_all_cophenetic': round(float(np.median(C[np.triu_indices(n, 1)])), 4),
               'nn_same_make_up': round(obs_nn, 3), 'leaf_runs': obs_runs}
        for kind, v in nul.items():
            v = np.array(v); row[f'nn_chance_{kind}'] = round(float(v[:, 1].mean()), 3)
            row[f'p_within_{kind}'] = (1 + int((v[:, 0] <= obs_w).sum())) / (NP + 1)
            row[f'p_nn_{kind}'] = (1 + int((v[:, 1] >= obs_nn).sum())) / (NP + 1)
        out.append(row)
T = pd.DataFrame(out); T.to_csv(f'{O}/lineage_tree_match.tsv', sep='\t', index=False)
# AB haplotypes: same cluster in both trees? cut each tree into k clusters and cross-tabulate
from scipy.cluster.hierarchy import fcluster
summ = {'haplotypes': n, 'make_up_counts': {k: int(v) for k, v in zip(*np.unique(lab, return_counts=True))}, 'permutations': NP}
ML = pd.read_csv(f'{O}/ml_hap_labels.tsv', sep='\t').set_index('src')   # step 19b: structure clusters (UPGMA cut at 16) and ML sequence clusters
assert set(rk) <= set(ML.index)
cs = fcluster(Zs, 16, 'maxclust'); assert all(ML.loc[r, 'struct_cluster'] == cs[i] for i, r in enumerate(rk)), 'structure clusters differ from step 19b'
cg = ML.loc[rk, 'seq_label'].to_numpy(dtype=object)
for g in ['AB', 'BB', 'B']:
    m = lab == g
    if not m.any(): continue
    ts = pd.Series(cs[m]).value_counts(); tg = pd.Series(cg[m]).value_counts()
    both = int(((cs == ts.index[0]) & (cg == tg.index[0]) & m).sum())
    summ[f'{g}_k16'] = {'largest_struct_cluster': int(ts.iloc[0]), 'struct_cluster_size': int((cs == ts.index[0]).sum()),
                         'largest_seq_cluster': int(tg.iloc[0]), 'seq_cluster': str(tg.index[0]), 'seq_cluster_size': int((cg == tg.index[0]).sum()), 'in_both': both}
summ['sequence_tree'] = 'maximum likelihood (step 04a)'
json.dump(summ, open(f'{O}/lineage_tree_match.json', 'w'), indent=1)
print(T.to_string()); print(json.dumps(summ, indent=1))
