# Pipeline step 21 (five_group_shuffles); ported from explore/structure_vs_sequence/18_five_group_shuffles.py (sha256 48de18d53e07, commit 5e2b3bb),
# paths rewritten to results/sv/ and figures/; run from the project root after the steps before it.
"""Step 18: the permutation tests of steps 02 and 14 rerun with labels shuffled only within each of the five
continental groups (AFR, AMR, EAS, EUR, SAS; ASL/unlabelled = 'other' and the two references kept as their own
strata), instead of within AFR / non-AFR. Same statistics, same inputs:
  Mantel (Spearman) of structure against SMN1-copy and SMN2-copy sequence distance (out/pair_distances.tsv.gz);
  copy nearest-neighbor label agreement (gene-body p-distance from out/copies.aln.fa, as step 02);
  lineage make-up nearest-neighbor agreement in both trees (as step 19c: structure UPGMA, sequence = the ML tree).
999 permutations each (1,999 for step 14's), seed 1. Writes out/five_group_shuffles.tsv."""
import numpy as np, pandas as pd
from scipy.stats import spearmanr
O = 'results/sv'; rng = np.random.default_rng(1)
L = pd.read_csv(f'{O}/lineage_haplotypes.tsv', sep='\t').set_index('src')
P = pd.read_csv(f'{O}/pair_distances.tsv.gz', sep='\t')
rows = sorted(set(P.a) | set(P.b)); ri = {r: i for i, r in enumerate(rows)}; n0 = len(rows)
def mat(col):
    M = np.full((n0, n0), np.nan); np.fill_diagonal(M, 0)
    ia, ib = P.a.map(ri).to_numpy(), P.b.map(ri).to_numpy(); M[ia, ib] = M[ib, ia] = P[col].to_numpy(); return M
S = mat('struct')
sp = np.array([L.loc[r, 'superpop'] if r in L.index else 'reference' for r in rows], dtype=object)
def strat_perm(st):
    p = np.arange(len(st))
    for g in sorted(set(st)):
        w = np.where(st == g)[0]; p[w] = w[rng.permutation(len(w))]
    return p
out = []
for col in ['seq_SMN1', 'seq_SMN2']:
    M = mat(col); k = np.where((~np.isnan(M)).sum(1) > 1)[0]
    Sk, Mk = S[np.ix_(k, k)], M[np.ix_(k, k)]; u = np.triu_indices(len(k), 1)
    obs = spearmanr(Sk[u], Mk[u]).statistic; null = []
    for _ in range(999):
        p = strat_perm(sp[k]); null.append(spearmanr(Sk[u], Mk[np.ix_(p, p)][u]).statistic)
    out.append({'test': 'mantel', 'what': col, 'n': len(k), 'observed': round(float(obs), 4), 'null_mean': round(float(np.mean(null)), 4),
                'p': (1 + sum(x >= obs for x in null)) / 1000})
    print(out[-1], flush=True)
# copy level: rebuild D exactly as step 02
names, seqs = [], []
for l in open(f'{O}/copies.aln.fa'):
    l = l.strip()
    if l.startswith('>'): names.append(l[1:]); seqs.append([])
    else: seqs[-1].append(l.upper())
A = np.array([np.frombuffer(''.join(s).encode(), dtype=np.uint8) for s in seqs])
_t = pd.read_csv(f'{O}/copy_tips.tsv', sep='\t').set_index('tip').loc[names]
keep_rows = (_t.copy_class != 'del7_8_SMN2like').to_numpy()
ref = np.where((_t.copy_class == 'SMN1').to_numpy() & _t.src.str.startswith('GCA_000001405').to_numpy())[0][0]
rpos = np.cumsum(A[ref] != ord('-')) - 1; rlen = rpos[-1] + 1
inside = (rpos >= 2000) & (rpos < rlen - 2000); colkeep = inside | ((A[ref] == ord('-')) & (rpos >= 1999) & (rpos < rlen - 2000))
A = A[keep_rows][:, colkeep]; names = [n for n, k in zip(names, keep_rows) if k]
assert colkeep.sum() == 28778
valid = np.isin(A, list(b'ACGT')); V = valid.astype(np.float32); nvalid = V @ V.T
var = np.array([len(set(A[valid[:, j], j])) > 1 for j in range(A.shape[1])])
same = np.zeros_like(nvalid)
for b in b'ACGT':
    X = (A[:, var] == b).astype(np.float32); same += X @ X.T
D = (V[:, var] @ V[:, var].T - same) / nvalid; np.fill_diagonal(D, 0)
lab = pd.read_csv(f'{O}/copy_labels.tsv', sep='\t').set_index('tip').loc[names].reset_index()
assert len(lab) == len(names)
for cls in ['SMN1', 'SMN2']:
    sel = np.where((lab.cls == cls).values)[0]
    Dc = D[np.ix_(sel, sel)].copy(); srcs = lab.src.to_numpy(dtype=object)[sel]
    Dc[srcs[:, None] == srcs[None, :]] = np.inf
    neigh = [np.where(Dc[i] == Dc[i].min())[0] for i in range(len(sel))]
    spc = lab.superpop.fillna('reference').to_numpy(dtype=object)[sel]
    for f in ['position', 'strand', 'naip3', 'env', 'arrangement']:
        y = lab[f].to_numpy(dtype=object)[sel]
        ok = np.array([v is not None and v == v for v in y])
        def agree(yy):
            return np.mean([np.mean(yy[ne[ok[ne]]] == yy[i]) for i, ne in enumerate(neigh) if ok[i] and ok[ne].any()])
        obs = agree(y); null = []
        for _ in range(999):
            yp = y.copy()
            for g in sorted(set(spc)):
                w = np.where(ok & (spc == g))[0]; yp[w] = y[w][rng.permutation(len(w))]
            null.append(agree(yp))
        out.append({'test': 'copy_neighbor', 'what': f'{cls}:{f}', 'n': int(ok.sum()), 'observed': round(float(obs), 3),
                    'null_mean': round(float(np.mean(null)), 3), 'p': (1 + sum(x >= obs for x in null)) / 1000})
        print(out[-1], flush=True)
# step 14: lineage make-up, nearest neighbor in each tree
G1 = mat('seq_SMN1'); k1 = np.where((~np.isnan(G1)).sum(1) > 1)[0]
rk = [rows[i] for i in k1]
mk = np.array([L.loc[r, 'composition'] if r in L.index else 'NA' for r in rk], dtype=object)
spk = sp[k1]
import sys; sys.path.insert(0, 'scripts'); from smnlib import ml_smn1_hap
Hml, _ = ml_smn1_hap(rk, O)   # sequence tree = the gene-body ML tree (step 04a), at the haplotype level, as step 19c
for tree, M in [('structure', S[np.ix_(k1, k1)]), ('SMN1 sequence', Hml)]:
    Mi = M.copy(); np.fill_diagonal(Mi, np.inf)
    neigh = [np.where(Mi[i] == Mi[i].min())[0] for i in range(len(k1))]
    for g in ['A', 'AA', 'AB', 'B', 'BB']:
        def nn(lb):
            m = lb == g; return np.mean([np.mean(lb[neigh[i]] == g) for i in np.where(m)[0]])
        obs = nn(mk); null = []
        for _ in range(1999):
            null.append(nn(mk[strat_perm(spk)]))
        out.append({'test': 'make_up_neighbor', 'what': f'{tree}:{g}', 'n': int((mk == g).sum()), 'observed': round(float(obs), 3),
                    'null_mean': round(float(np.mean(null)), 3), 'p': (1 + sum(x >= obs for x in null)) / 2000})
        print(out[-1], flush=True)
pd.DataFrame(out).to_csv(f'{O}/five_group_shuffles.tsv', sep='\t', index=False)
