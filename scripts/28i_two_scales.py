"""Pipeline step 28i (after 28b, whose arm types it reads): structure against sequence at two bundle scales.
Coarse = the study's decomposition (results/decomp/smn.bed, step 00; distance results/sv/smn_struct.dist, step 01).
Fine   = the clinical study's finer decomposition (results/bundles/smn.bed; distance results/sv/fine/smn_struct.dist,
         step 01b).
Sequence distances are step 07's haplotype copy-set distances (results/sv/pair_distances.tsv.gz: all copies,
SMN1 copies, SMN2 copies), so only the structural side differs between the scales.
Tests, each run identically at both scales (same rows, same permutations: one generator per test, seed 1):
  Mantel (Spearman, 999 permutations) for all copies, SMN1, SMN2, SMN1 within AFR / outside AFR;
  cophenetic correlation on the haplotypes with an SMN1 copy: the structure tree (UPGMA, optimal leaf order) against
  the gene-body maximum-likelihood tree (results/sv/mltree.treefile, step 28), whose haplotype distance is the
  patristic distance between the haplotypes' SMN1 copies, combined over copies as step 07 combines p-distances
  (mean of each copy's nearest copy, both ways); the UPGMA sequence tree is kept as 'upgma_seq_tree' in the json;
  SMN1 Mantel within one arm type (C1, C0, T0, T1; arm types from step 28b, results/sv/arm_element_haps.tsv),
  so that the arm-type split itself cannot drive the correlation.
Writes results/sv/two_scales.tsv and results/sv/two_scales.json."""
import numpy as np, pandas as pd, json
from scipy.cluster.hierarchy import average, cophenet, optimal_leaf_ordering
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr
O = 'results/sv'
AFR = {'ACB', 'ASW', 'ESN', 'GWD', 'LWK', 'MSL', 'YRI', 'MKK'}
pdist = pd.read_csv(f'{O}/pair_distances.tsv.gz', sep='\t')
rows = sorted(set(pdist.a) | set(pdist.b)); ri = {r: i for i, r in enumerate(rows)}; n = len(rows)
ia = pdist.a.map(ri).to_numpy(); ib = pdist.b.map(ri).to_numpy()
def mat(vals):
    M = np.full((n, n), np.nan); np.fill_diagonal(M, 0); M[ia, ib] = vals; M[ib, ia] = vals; return M
G, G1, G2 = mat(pdist.seq_all_copies.to_numpy()), mat(pdist.seq_SMN1.to_numpy()), mat(pdist.seq_SMN2.to_numpy())
scoarse = mat(pdist.struct.to_numpy())
def read_dist(path):
    st = pd.read_csv(path, sep=' ', header=None, usecols=[0, 1, 2], names=['a', 'b', 'd'])
    st = st[st.a.isin(ri) & st.b.isin(ri)]
    S = np.full((n, n), np.nan)
    S[st.a.map(ri).to_numpy(), st.b.map(ri).to_numpy()] = st.d.to_numpy()
    assert not np.isnan(S).any(); S = (S + S.T) / 2; np.fill_diagonal(S, 0); return S
S = {'coarse': read_dist(f'{O}/smn_struct.dist'), 'fine': read_dist(f'{O}/fine/smn_struct.dist')}
# the coarse matrix must be the one step 07 used
assert np.allclose(S['coarse'], scoarse, atol=1e-6), 'step 07 struct column differs from results/sv/smn_struct.dist'
hc = pd.read_csv('tables/haplotype_copy_counts.tsv', sep='\t'); pop = dict(zip(hc.src, hc.population))
afr = np.array([pop.get(r) in AFR for r in rows])
ah = pd.read_csv(f'{O}/arm_element_haps.tsv', sep='\t').set_index('src')
ctype = np.array([ah.cen_type.get(r, '') for r in rows]); ttype = np.array([ah.tel_type.get(r, '') for r in rows])
def mantel(Sx, M, keep, seed):
    rng = np.random.default_rng(seed)
    k = np.where(keep & ((~np.isnan(M)).sum(1) > 1))[0]
    Sk, Mk = Sx[np.ix_(k, k)], M[np.ix_(k, k)]; u = np.triu_indices(len(k), 1)
    obs = spearmanr(Sk[u], Mk[u]).statistic; ge = 0
    for _ in range(999):
        p = rng.permutation(len(k)); ge += spearmanr(Sk[u], Mk[np.ix_(p, p)][u]).statistic >= obs
    return int(len(k)), float(obs), (1 + ge) / 1000
all_ = np.ones(n, bool)
tests = [('all copies', G, all_), ('SMN1', G1, all_), ('SMN2', G2, all_), ('SMN1, AFR', G1, afr), ('SMN1, outside AFR', G1, ~afr),
         ('SMN1, C1 arms', G1, ctype == 'C1'), ('SMN1, C0 arms', G1, ctype == 'C0'),
         ('SMN1, T0 arms', G1, ttype == 'T0'), ('SMN1, T1 arms', G1, ttype == 'T1')]
out = []
for ti, (name, M, keep) in enumerate(tests):
    rec = {'test': name}
    for sc in ('coarse', 'fine'):
        h, r, p = mantel(S[sc], M, keep, seed=1 + ti)
        rec['haplotypes'] = h; rec[f'{sc}_spearman'] = round(r, 4); rec[f'{sc}_p'] = p
    out.append(rec)
k1 = np.where((~np.isnan(G1)).sum(1) > 1)[0]; G1k = G1[np.ix_(k1, k1)]
Zg = optimal_leaf_ordering(average(squareform(G1k)), squareform(G1k))
# the ML tree's patristic distances between SMN1 copies, combined per haplotype pair as step 07 does (smnlib.ml_smn1_hap)
import sys; sys.path.insert(0, 'scripts')
from smnlib import ml_smn1_hap
hs = [rows[i] for i in k1]; H, _ = ml_smn1_hap(hs, O); m = len(hs)
_tp = pd.read_csv(f'{O}/copy_tips.tsv', sep='\t'); ntips = int(((_tp.copy_class == 'SMN1') & _tp.src.isin(hs)).sum())
_u = np.triu_indices(m, 1)
ml_vs_p = round(float(spearmanr(H[_u], G1k[_u]).statistic), 4)
assert ml_vs_p > 0.9, f'ML patristic and p-distance disagree ({ml_vs_p})'
coph, coph_upgma = {}, {}
for sc in ('coarse', 'fine'):
    Sk = S[sc][np.ix_(k1, k1)]; Zs = optimal_leaf_ordering(average(squareform(Sk)), squareform(Sk))
    coph[sc] = round(float(np.corrcoef(cophenet(Zs), H[_u])[0, 1]), 4)
    coph_upgma[sc] = round(float(np.corrcoef(cophenet(Zs), cophenet(Zg))[0, 1]), 4)
out.append({'test': 'cophenetic correlation (SMN1 haplotypes)', 'haplotypes': int(len(k1)), 'coarse_spearman': coph['coarse'], 'coarse_p': None, 'fine_spearman': coph['fine'], 'fine_p': None})
tab = pd.DataFrame(out); tab.to_csv(f'{O}/two_scales.tsv', sep='\t', index=False)
# the scales' structural distances against each other
u = np.triu_indices(n, 1)
summ = {'haplotypes': n, 'cophenetic_seq_tree': 'maximum likelihood', 'upgma_seq_tree': coph_upgma, 'ml_patristic_vs_pdistance_spearman': ml_vs_p, 'ml_tips': ntips, 'coarse_vs_fine_struct_spearman': round(float(spearmanr(S['coarse'][u], S['fine'][u]).statistic), 4),
        'tests': out}
json.dump(summ, open(f'{O}/two_scales.json', 'w'), indent=1)
print(tab.to_string(index=False)); print('coarse vs fine structural distance, Spearman', summ['coarse_vs_fine_struct_spearman'])

# ---- the parameters of the two scales, read from the files themselves (researcher, 2026-11-15: document them)
import re, subprocess
help_ = subprocess.run(['pangea', 'cohort', 'bundle', 'decomp', '--help'], capture_output=True, text=True).stdout
defaults = {}; flag = None
for l in help_.splitlines():
    m = re.match(r'^\s+(?:-(\w),? ?)?(?:--([\w-]+))?\s*<\w+>\s*$', l)
    if m and (m.group(1) or m.group(2)): flag = m.group(2) or m.group(1); continue
    d = re.search(r'\[default: ([^\]]+)\]', l)
    if d and flag: defaults[flag] = d.group(1); flag = None
for key in ['w', 'k', 'r', 'min-span', 'min-cov', 'min-branch-size', 'bundle-length-cutoff', 'bundle-merge-distance']:
    assert key in defaults, f'no default parsed for {key}'
def cmdline(path):
    for l in open(path):
        if l.startswith('# cmd:'): return l[len('# cmd:'):].strip()
    raise SystemExit(f'no # cmd: line in {path}')
def given(cmd):
    t = cmd.split(); g = {}
    for i, x in enumerate(t):
        if x.startswith('-') and x not in ('-f', '-o', '-b') and i + 1 < len(t): g[x.lstrip('-')] = t[i + 1]
    return g
prow = []; ver = subprocess.run(['pangea', '--version'], capture_output=True, text=True).stdout.strip()
for sc, bed in (('coarse', 'results/decomp/smn.bed'), ('fine', 'results/bundles/smn.bed')):
    c = cmdline(bed); g = given(c)
    b = pd.read_csv(bed, sep='\t', comment='#', header=None, usecols=[0, 3])
    ids = b[3].astype(str).str.split(':').str[0]
    rec = {'scale': sc, 'bed': bed, 'decomp_cmd': c}
    for key in ['w', 'k', 'r', 'min-span', 'min-cov', 'min-branch-size', 'bundle-length-cutoff', 'bundle-merge-distance']:
        rec[key] = g.get(key, defaults[key]); rec[key + '_source'] = 'given' if key in g else 'default'
    rec['repeat_units_cmd'] = cmdline(bed.replace('.bed', '.repeat_units.tsv'))
    rec['dist_cmd'] = f"pangea cohort bundle dist -b {bed} --top-n 50"
    rec['rows'] = int(b[0].nunique()); rec['segments'] = int(len(b)); rec['bundle_ids'] = int(ids.nunique())
    rec['pangea_version_now'] = ver
    prow.append(rec)
P = pd.DataFrame(prow); P.to_csv(f'{O}/decomp_params.tsv', sep='\t', index=False)
summ['params'] = prow; json.dump(summ, open(f'{O}/two_scales.json', 'w'), indent=1)
print(P.drop(columns=['decomp_cmd', 'repeat_units_cmd', 'dist_cmd', 'pangea_version_now']).T.to_string())
