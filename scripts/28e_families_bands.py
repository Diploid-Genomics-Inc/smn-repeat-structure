"""Pipeline step 28e: (1) what each repeat-unit family is made of, in arm elements; (2) the bands and structure
clusters against arm form.
(1) For every repeat-unit occurrence on the 270 haplotypes with the head-to-head A pair (results/decomp/smn.repeat_units.tsv),
    the bp share of each arm element it covers (results/sv/arm_element_segments.tsv, step 28b; 'Acore' = segments within
    the arm's A core, a_end_kb from arm_elements.tsv). Per family (>= 10 occurrences): mean share of each element and the
    commonest element combination (elements with share >= 0.1). -> results/sv/unit_family_elements.tsv
(2) Arm form (cenSMN C1/C0 | telSMN T0/T1, arm_element_haps.tsv) cross-tabulated against the structure clusters and
    sequence labels of step 19b (ml_hap_labels.tsv) and its bands (ml_cluster_pairs_members.tsv); family G carriers
    against telSMN type. -> results/sv/armform_struct.tsv, armform_bands.tsv, armform_summary.json"""
import json
import numpy as np, pandas as pd
from scipy.stats import fisher_exact
O = 'results/sv/'
SG = pd.read_csv(O + 'arm_element_segments.tsv', sep='\t')
A = pd.read_csv(O + 'arm_elements.tsv', sep='\t')[['src', 'side', 'a_end_kb']]
AX = pd.read_csv(O + 'palindrome_diagram_rows.tsv', sep='\t').set_index('src').axis
SG = SG.merge(A, on=['src', 'side'])
SG['lab'] = np.where(SG.a1 <= SG.a_end_kb + 0.001, 'Acore', SG.el.replace('.', 'other'))
ax = SG.src.map(AX)
SG['bgn'] = np.where(SG.side == 'cen', ax - SG.a1 * 1e3, ax + SG.a0 * 1e3)
SG['end'] = np.where(SG.side == 'cen', ax - SG.a0 * 1e3, ax + SG.a1 * 1e3)
RU = pd.read_csv('results/decomp/smn.repeat_units.tsv', sep='\t', comment='#')
RUp = RU[RU.row.isin(AX.index)]
by = {s: g for s, g in SG.groupby('src')}
rows = []
for r in RUp.itertuples():
    s = by.get(r.row)
    if s is None: continue
    s = s[(s.end > r.bgn) & (s.bgn < r.end)]
    ov = (np.minimum(s.end, r.end) - np.maximum(s.bgn, r.bgn)).clip(lower=0)
    bp = ov.groupby(s.lab).sum(); tot = bp.sum()
    if tot == 0: continue
    rows.append(dict(src=r.row, family=r.family, span_kb=r.span_bp / 1e3, **(bp / tot).to_dict()))
D = pd.DataFrame(rows).fillna(0)
EL = ['Acore', 'K', 'Xi', 'S', 'Xd', 'E', 'Fl', 'other']
for c in EL:
    if c not in D: D[c] = 0.0
NAME = {'Acore': 'A core', 'K': 'core end', 'Xi': 'Xi', 'S': 'S', 'Xd': 'Xd', 'E': 'E', 'Fl': 'flank', 'other': 'other'}
D['combo'] = D[EL].apply(lambda x: ' + '.join(NAME[c] for c in EL if x[c] >= 0.1), axis=1)
out = []
for f, g in D.groupby('family'):
    if len(g) < 10: continue
    vc = g.combo.value_counts()
    out.append(dict(family=f, occurrences=len(g), haplotypes=g.src.nunique(), **{f'share_{c}': round(g[c].mean(), 3) for c in EL},
                    top_combo=vc.index[0], top_combo_n=int(vc.iloc[0]), main_element=NAME[g[EL].mean().idxmax()]))
F = pd.DataFrame(out).sort_values('occurrences', ascending=False)
F.to_csv(O + 'unit_family_elements.tsv', sep='\t', index=False)
# (2) arm form against structure clusters and bands
H = pd.read_csv(O + 'arm_element_haps.tsv', sep='\t'); H['form'] = H.cen_type + '|' + H.tel_type
L = pd.read_csv(O + 'ml_hap_labels.tsv', sep='\t')
M = pd.read_csv(O + 'ml_cluster_pairs_members.tsv', sep='\t')
X = L.merge(H[['src', 'form', 'tel_type', 'cen_type']], on='src', how='left'); X['form'] = X.form.fillna('no pair')
st = pd.crosstab(X.struct_cluster.astype(int), X.form); st.to_csv(O + 'armform_struct.tsv', sep='\t')
B = M.merge(H[['src', 'form']], on='src', how='left'); B['form'] = B.form.fillna('no pair')
bt = pd.crosstab(B.band, B.form); bt.to_csv(O + 'armform_bands.tsv', sep='\t')
# the two big structure clusters: the bulk's cluster and the cluster shared by the S1-3 and B1+B2 bands
cb = lambda b: int(M[M.band == b].struct_cluster.mode()[0])
bulk, mixc = cb('Bulk'), cb('AB')   # the B1+B2 band's cluster (under the clinical decomposition the S1-3 band shared it)
if 'Three' in set(M.band): assert cb('Three') == mixc
P = X[X.form != 'no pair']
def tt(c):
    g = P[P.struct_cluster == c]; return int((g.tel_type == 'T1').sum()), len(g)
t2b, nb = tt(bulk); t2m, nm = tt(mixc)
p = fisher_exact([[t2m, nm - t2m], [t2b, nb - t2b]])[1]
T2all = int((H.tel_type == 'T1').sum()); T1all = int((H.tel_type == 'T0').sum())
T1bulk = int(((P.struct_cluster == bulk) & (P.tel_type == 'T0')).sum())
# the family that marks the T1 arm: for each family, its share of T1 and of T0 pair haplotypes; the largest gap wins
FT = {f: (len(set(g.row) & set(H[H.tel_type == 'T1'].src)) / max(1, int((H.tel_type == 'T1').sum())),
          len(set(g.row) & set(H[H.tel_type == 'T0'].src)) / max(1, int((H.tel_type == 'T0').sum()))) for f, g in RU.groupby('family')}
SIG = max(FT, key=lambda f: FT[f][0] - FT[f][1])
print('family vs telSMN arm type (share of T1, share of T0):', {f: tuple(round(x, 2) for x in v) for f, v in sorted(FT.items())}, 'signature', SIG)
g_src = set(RU[RU.family == SIG].row)
H['G'] = H.src.isin(g_src)
nopair_band = B[(B.form == 'no pair')]
S = dict(bulk_cluster=bulk, mixed_cluster=mixc, sig_family=SIG, family_T1_T0={f: [round(a, 3), round(b, 3)] for f, (a, b) in FT.items()},
         mixed_cluster_size=int((X.struct_cluster == mixc).sum()), bulk_cluster_size=int((X.struct_cluster == bulk).sum()),
         mixed_pair=nm, mixed_pair_T1=t2m, bulk_pair=nb, bulk_pair_T1=t2b, bulk_pair_T0=T1bulk, fisher_p_T1=p,
         T1_all=T2all, T0_all=T1all, T1_in_mixed=t2m,
         mixed_C0T1=int(((P.struct_cluster == mixc) & (P.form == 'C0|T1')).sum()),
         mixed_C1T1=int(((P.struct_cluster == mixc) & (P.form == 'C1|T1')).sum()),
         bands={b: {k: int(v) for k, v in r.items()} for b, r in bt.iterrows()},
         G_on_T1=[int((H.G & (H.tel_type == 'T1')).sum()), T2all], G_on_T0=[int((H.G & (H.tel_type == 'T0')).sum()), T1all],
         G_haps_all=len(g_src), G_haps_nopair=len(g_src - set(H.src)),
         S13_band_nopair=int((nopair_band.band == 'Three').sum()),
         S13_band_nopair_G=int(((nopair_band.band == 'Three') & nopair_band.src.isin(g_src)).sum()))
# the cenSMN arm type in the same two clusters (coarse scale: the clusters split by C, not by T)
_cb = P[P.struct_cluster == bulk]; _cm = P[P.struct_cluster == mixc]
_c0b, _c0m = int(_cb.form.str.startswith('C0').sum()), int(_cm.form.str.startswith('C0').sum())
S.update(bulk_pair_C0=_c0b, bulk_pair_C1=len(_cb) - _c0b, mixed_pair_C0=_c0m, mixed_pair_C1=len(_cm) - _c0m,
         fisher_p_C=float(fisher_exact([[_c0m, len(_cm) - _c0m], [_c0b, len(_cb) - _c0b]])[1]))
json.dump(S, open(O + 'armform_summary.json', 'w'), indent=1, sort_keys=True)
print(F.to_string()); print(st.to_string()); print(bt.to_string()); print(json.dumps(S, indent=1))
