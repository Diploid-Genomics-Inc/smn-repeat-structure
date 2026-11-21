# Pipeline step 20 (s13_units); ported from explore/structure_vs_sequence/17_s13_units.py (sha256 1c606b06815c, commit 5e2b3bb),
# paths rewritten to results/sv/ and figures/; run from the project root after the steps before it.
"""Step 17: which large repeat-unit families set the S1-3 band of the Fig. 2 tanglegram apart from the central
"2- 1+" cluster, and does that unit go with SMN1 haplogroup S1-3 across the cohort?
Band and central cluster: the S1-3 band and the bulk of step 19b (ML-clade bands, results/sv/ml_cluster_pairs_members.tsv).
Unit families per haplotype: results/decomp/haplotype_units.tsv (families of pangea's repeat units, centromeric to telomeric).
Haplogroup per haplotype: out/chen_hg_copies.tsv (S1-3 if any SMN1 copy is assigned S1-3 alone).
Writes out/s13_unit_families.tsv (presence of each family, band vs central, Fisher) and out/s13_summary.json."""
import json, collections, pandas as pd
from scipy.stats import fisher_exact
O = 'results/sv'
# bands from step 19b (sequence side = supported ML clades; researcher's decision 2026-11-14); before that, step 19's
# UPGMA pairs (structure 3 x sequence 15, structure 6 x sequence 8), which hold the same 33 band haplotypes
M = pd.read_csv(f'{O}/ml_cluster_pairs_members.tsv', sep='\t')
band = M[M.band == 'Three']
cen = M[M.band == 'Bulk']
assert len(band) == 33 and len(cen) == 209
U = pd.read_csv('results/decomp/haplotype_units.tsv', sep='\t').set_index('src')
fams = sorted({t[0] for u in U.units.dropna() for t in u.split()})
def present(src, f):
    return f in {t[0] for t in str(U.loc[src, 'units']).split()}
rows = []
for label, sub in [('all', None), ('2- 1+ only', '2- 1+')]:
    b = band if sub is None else band[band.arrangement == sub]
    c = cen if sub is None else cen[cen.arrangement == sub]
    for f in fams:
        nb = sum(present(s, f) for s in b.src); nc = sum(present(s, f) for s in c.src)
        p = fisher_exact([[nb, len(b) - nb], [nc, len(c) - nc]]).pvalue
        rows.append({'subset': label, 'family': f, 'band_with': nb, 'band_n': len(b), 'central_with': nc, 'central_n': len(c), 'fisher_p': p})
T = pd.DataFrame(rows); T.to_csv(f'{O}/s13_unit_families.tsv', sep='\t', index=False)
# cohort-wide: family G against S1-3, over resolved haplotypes carrying an SMN1 copy
H = pd.read_csv(f'{O}/chen_hg_copies.tsv', sep='\t'); H['src'] = H.copy_id.str.split('|').str[0]
s13 = set(H[H.family == 'S1-3'].src)
L = pd.read_csv(f'{O}/lineage_haplotypes.tsv', sep='\t')
tree = sorted(L[(L.resolved == True) & (L.nA + L.nB > 0)].src)
g = {s: present(s, 'G') for s in tree if s in U.index}
a = sum(1 for s in g if g[s] and s in s13); b_ = sum(1 for s in g if g[s] and s not in s13)
c_ = sum(1 for s in g if not g[s] and s in s13); d = sum(1 for s in g if not g[s] and s not in s13)
R = pd.read_csv('results/decomp/smn.repeat_units.tsv', sep='\t', comment='#')
Gu = R[R.family == 'G']
summ = {'band': len(band), 'central': len(cen), 'band_G': int(sum(present(s, 'G') for s in band.src)),
        'central_G': int(sum(present(s, 'G') for s in cen.src)),
        'G_and_S13': a, 'G_not_S13': b_, 'S13_not_G': c_, 'neither': d, 'haplotypes': len(g),
        'G_S13_fisher_p': fisher_exact([[a, b_], [c_, d]]).pvalue,
        'G_units': int(len(Gu)), 'G_rows': int(Gu['row'].nunique()), 'G_span_median_bp': float(Gu.span_bp.median()),
        'band_S13': int(sum(s in s13 for s in band.src))}
json.dump(summ, open(f'{O}/s13_summary.json', 'w'), indent=1)
print(json.dumps(summ, indent=1)); print(T[T.fisher_p < 0.01].to_string())
# who else carries family G
Li = L.set_index('src')
Go = [s for s in g if g[s] and s not in s13]
B = pd.DataFrame({'src': Go, 'arrangement': [Li.loc[s, 'arrangement'] for s in Go], 'composition': [Li.loc[s, 'composition'] for s in Go], 'superpop': [Li.loc[s, 'superpop'] for s in Go],
                  'hg': ['; '.join(sorted(set(H[(H.src == s) & (H.lineage.isin(['A', 'B']))].family))) for s in Go]})
B.to_csv(f'{O}/s13_G_other_haplotypes.tsv', sep='\t', index=False)
for c in ['arrangement', 'composition', 'superpop', 'hg']: print(c, B[c].value_counts().head(8).to_dict())
# among G carriers: families C and D by haplogroup group
def grp(s):
    f = set(H[(H.src == s) & (H.lineage.isin(['A', 'B']))].family)
    if 'S1-3' in f: return 'S1-3'
    if f & {'S1-2', 'S1-2c'}: return 'S1-2'
    if f & {'S1-8', 'S1-9', 'S1-9c', 'S1-9d', 'S1-6'}: return 'S1-8/9/6'
    return 'other'
Gs = [s for s in g if g[s]]
C = pd.DataFrame({'src': Gs, 'group': [grp(s) for s in Gs], 'C': [present(s, 'C') for s in Gs], 'D': [present(s, 'D') for s in Gs]})
CT = C.groupby('group').agg(n=('src', 'size'), with_C=('C', 'sum'), with_D=('D', 'sum')).reset_index()
CT.to_csv(f'{O}/s13_G_carriers_CD.tsv', sep='\t', index=False); print(CT.to_string())
