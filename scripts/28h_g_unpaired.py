"""Pipeline step 28h (from exploration g01, 2026-11-14; researcher: 'G result: add it to S3.5? Ok'; researcher: 'yes' to asking whether the 79 haplotypes that carry family G but no
A+ A- pair carry a T1-like arm; and: the elements defined from the A+ A- axis may not apply to unpaired haplotypes).
Axis-free reading: every S site (S bundles joined within 30 kb) in results/decomp/smn.bed is read from its
nearest SMN copy (results/decomp/copy_units.tsv, mid_fwd) outward, by bundle ID only (no strand, no axis):
  E between copy and S; Xi-outer IDs (step 28b) within 80 kb before S (near) and after S (beyond);
  another SMN copy within 450 kb beyond S (the site lies between two copies: the main axis).
Class: axis | T1-like (E between, Xi-outer on both sides of S) | C1-like (no E between, both sides) |
       T0-like (Xi-outer near side only) | other.
Validated against step 28b on the 270 pair haplotypes (arm_element_spacers.tsv sites, arm_element_haps.tsv types).
Writes results/sv/g_s_sites.tsv, g_validation.tsv, g_unpaired_calls.tsv and g_unpaired.json."""
import pandas as pd, numpy as np
import json as _j
S_, XO, E_ = {20}, {13, 16, 17, 21, 24, 3}, {7, 14, 18, 19, 23, 4}   # step 28b's S, Xi-outer and E bundles (study decomposition)
SIG = _j.load(open('results/sv/armform_summary.json'))['sig_family']   # the family that marks the T1 arm (step 28e)
B = pd.read_csv('results/decomp/smn.bed', sep='\t', comment='#', header=None, names=['row', 'bgn', 'end', 'tag'])
B['bid'] = B.tag.str.split(':').str[0].astype(int); B['mid'] = (B.bgn + B.end) / 2
CU = pd.read_csv('results/decomp/copy_units.tsv', sep='\t')
rows = []
for row, b in B.groupby('row'):
    cps = np.sort(CU[CU.src == row].mid_fwd.values.astype(float))
    if not len(cps): continue
    s = b[b.bid.isin(S_)].sort_values('bgn'); sites = []
    for r in s.itertuples():
        if sites and r.bgn - sites[-1][1] < 30000: sites[-1][1] = max(sites[-1][1], r.end)
        else: sites.append([r.bgn, r.end])
    for s0, s1 in sites:
        sm = (s0 + s1) / 2; c = cps[np.argmin(np.abs(cps - sm))]; d = 1 if sm > c else -1
        lo, hi = sorted((c, sm))
        e_between = bool(((b.mid > lo) & (b.mid < hi) & b.bid.isin(E_)).any())
        near = (s0 - 80000, s0) if d == 1 else (s1, s1 + 80000); beyond = (s1, s1 + 80000) if d == 1 else (s0 - 80000, s0)
        xo_near = bool(((b.mid > near[0]) & (b.mid < near[1]) & b.bid.isin(XO)).any())
        xo_beyond = bool(((b.mid > beyond[0]) & (b.mid < beyond[1]) & b.bid.isin(XO)).any())
        cb = cps[(cps - sm) * d > 0]; copy_beyond = bool(len(cb) and np.min(np.abs(cb - sm)) < 450000)
        cls = ('axis' if copy_beyond else 'T1-like' if (xo_near and xo_beyond and e_between) else 'C1-like' if (xo_near and xo_beyond)
               else 'T0-like' if xo_near else 'other')
        rows.append(dict(row=row, s0=s0, s1=s1, copy_kb=round(abs(sm - c) / 1e3, 1), e_between=e_between, xo_near=xo_near, xo_beyond=xo_beyond, copy_beyond=copy_beyond, cls=cls))
T = pd.DataFrame(rows); T.to_csv('results/sv/g_s_sites.tsv', sep='\t', index=False)
# validation on pair haplotypes: match each S site to step 28b's site label, and the arm type
SP = pd.read_csv('results/sv/arm_element_spacers.tsv', sep='\t'); EH = pd.read_csv('results/sv/arm_element_haps.tsv', sep='\t').set_index('src')
v = []
for r in T[T.row.isin(EH.index)].itertuples():
    m = SP[(SP.src == r.row) & (SP.end > r.s0) & (SP.bgn < r.s1)]
    site = m.site.iloc[0] if len(m) else 'none'
    typ = 'axis' if site == 'axis' else (EH.cen_type[r.row] if site == 'cen' else EH.tel_type[r.row] if site == 'tel' else 'none')
    v.append(dict(row=r.row, site=site, truth=typ, cls=r.cls))
V = pd.DataFrame(v); V.to_csv('results/sv/g_validation.tsv', sep='\t', index=False)
print(pd.crosstab(V.truth, V.cls, margins=True))
H = V.groupby('row').agg(t2=('cls', lambda x: (x == 'T1-like').any())); H['tel'] = EH.tel_type.reindex(H.index)
print('pair haps: any T1-like site vs tel type'); print(pd.crosstab(H.tel, H.t2))
PU = pd.read_csv('results/sv/palindrome_units.tsv', sep='\t'); RU = pd.read_csv('results/decomp/smn.repeat_units.tsv', sep='\t', comment='#')
unp = {r for r in PU[~PU.has_pair].src if not r.startswith(('GCA_', 'chm13'))};  # donors only
g = set(RU[RU.family == SIG].row)
U = T[T.row.isin(unp)].groupby('row').cls.apply(lambda x: ' '.join(sorted(x))).rename('sites').to_frame()
U['G'] = U.index.isin(g); U['arr'] = PU.set_index('src').arr.reindex(U.index); U['t2'] = U.sites.str.contains('T1-like')
U.to_csv('results/sv/g_unpaired_calls.tsv', sep='\t')
print('unpaired haps with an S site:', len(U), 'of', len(unp)); print(pd.crosstab(U.G, U.t2, margins=True))
print(U[U.G].groupby('arr').t2.agg(['sum', 'size']).sort_values('size', ascending=False))

import json
L = pd.read_csv('results/sv/ml_hap_labels.tsv', sep='\t').set_index('src')
U['seq_label'] = L.seq_label.reindex(U.index)
UG = U[U.G]
bylab = {k: [int(v.t2.sum()), int(len(v))] for k, v in UG.groupby('seq_label')}
hv = pd.crosstab(H.tel, H.t2)
out = dict(sig_family=SIG, pair_T1_t2like=[int(hv.loc['T1', True]), int(hv.loc['T1'].sum())], pair_T0_t2like=[int(hv.loc['T0', True]), int(hv.loc['T0'].sum())],
           unpaired=len(unp), unpaired_with_S=len(U), G_unpaired=int(U.G.sum()), G_unpaired_t2like=int(UG.t2.sum()),
           noG_unpaired_t2like=[int(U[~U.G].t2.sum()), int((~U.G).sum())], G_unpaired_by_label=bylab,
           S12_G_sites=UG[UG.seq_label == 'S1-2'].sites.value_counts().head(5).to_dict())
assert out['G_unpaired'] == len(g & unp & set(U.index)), 'G haplotypes without the pair but with no S site'
json.dump(out, open('results/sv/g_unpaired.json', 'w'), indent=1); print(json.dumps(out, indent=1))
