"""Pipeline step 28b: the elements of each arm of the head-to-head family-A pair, read from the bundle decomposition.
Ported from explore/a_break/a29, a30, a34-a36 (commits fec1214-0af8629); replaces the long/short unit call.
For each of the haplotypes with an A+ A- pair (step 22's rows; axis = midpoint of the junction), every bundle
segment of results/decomp/smn.bed is put on the cenSMN arm (left of the axis) or the telSMN arm (right), by its
midpoint, measured in kb from the axis, and read OUTWARD: its token is bundle id + strand as seen walking away from
the axis. Tokens are named by element. The sets are the clinical decomposition's elements (exploration a25/a26/a30) carried over
to the study's decomposition (step 00) by base overlap on the same windows (explore/new_bundles/n02 and the inner/outer
split of Xi): every token below has >= 0.7 of its labelled bases in its element; all but 3 and 4 have 1.0.
  K   end of the shared A core       8- 9- 10-
  Xi  inverted copy of the core's inner part   1- 5- (inner half) 13- 16- 17- 21- 24- 3- (outer half)
  S   palindrome spacer              bundle 20, either strand
  Xd  Xi's outer half again, direct  13+ 16+ 17+ 21+ 24+ 3+
  E   GTF2H2-bearing stretch         7+ 14+ 18+ 19+ 23- 4+ (4+ also runs into the flank: 0.4-0.7 of its labelled bases are E)
  Fl  the window's unique flank      bundles 2 6 11, either strand
Bundle 15 lies in the ~20 kb tandem arrays (0.84 of its bases) and, like any other bundle, is not an element.
A = the arm from the axis to the outer end of the first K run beyond 200 kb (whatever bundles lie inside).
Layout = the elements of segments beyond A, in order, runs collapsed, other bundles skipped.
Arm types, by whether the arm carries Xd beyond A (checked: Xd is never a marginal stretch, >= 20 kb wherever present):
  cenSMN C1 = carries X (Xi, S and Xd; Xi and Xd are checked to come together); C0 = carries none of Xi, S, Xd.
  telSMN T0 = Xi and S but no Xd;  T1 = Xi, S and Xd (Xi's outer half, S, Xd, then Xi's inner half further out).
Each arm's layout is kept; the commonest layout of each type (CANON) is checked.
Also per arm: the SMN copy (step 22), tandem arrays (step 28, units >= 10 kb), S occurrences, NAIP/GTF2H2 matches
(results/copies/features.tsv, the window carrying the palindrome copies) and the element each lies in.
Out: results/sv/arm_elements.tsv (one row per arm), arm_element_segments.tsv, arm_element_genes.tsv,
arm_element_arrays.tsv, arm_elements_summary.json."""
import re, json, numpy as np, pandas as pd
O = 'results/sv/'
EL = {}
for t in '8- 9- 10-'.split(): EL[t] = 'K'
for t in '1- 5- 13- 16- 17- 21- 24- 3-'.split(): EL[t] = 'Xi'
XIN = set('1- 5-'.split())
for t in '13+ 16+ 17+ 21+ 24+ 3+'.split(): EL[t] = 'Xd'
for t in '20+ 20-'.split(): EL[t] = 'S'
for t in '7+ 14+ 18+ 19+ 23- 4+'.split(): EL[t] = 'E'
for b in '2 6 11'.split(): EL[b + '+'] = EL[b + '-'] = 'Fl'
SBID = [20]                         # the bundle(s) of S
NAMES = {'K': 'end of the core', 'Xi': 'Xi, inverted copy of the core\'s inner part', 'S': 'S, palindrome spacer', 'Xd': 'Xd, outer half of Xi again, opposite to Xi (direct, as the core)',
         'E': 'E, carries GTF2H2', 'Fl': 'unique flank'}
pd.DataFrame([dict(element=e, name=NAMES[e], tokens=' '.join(t for t in EL if EL[t] == e)) for e in ['K', 'Xi', 'S', 'Xd', 'E', 'Fl']]).to_csv(O + 'arm_element_sets.tsv', sep='\t', index=False)
AX = pd.read_csv(O + 'palindrome_diagram_rows.tsv', sep='\t').set_index('src')
bed = pd.read_csv('results/decomp/smn.bed', sep='\t', comment='#', header=None, names=['row', 'bgn', 'end', 'name'])
f = bed.name.str.split(':', expand=True); bed['bid'] = f[0].astype(int); bed['d'] = f[2].astype(int)
segs = []
for src, ax in AX.axis.items():
    g = bed[bed.row == src]
    mid = (g.bgn + g.end) / 2
    for b, d, bg, en, m in zip(g.bid, g.d, g.bgn, g.end, mid):
        if m < ax: side, a0, a1, tok = 'cen', (ax - en) / 1e3, (ax - bg) / 1e3, f'{b}{"+-"[d]}'
        else: side, a0, a1, tok = 'tel', (bg - ax) / 1e3, (en - ax) / 1e3, f'{b}{"+-"[1 - d]}'
        segs.append(dict(src=src, side=side, a0=round(a0, 3), a1=round(a1, 3), bid=b, tok=tok, el=EL.get(tok, '.')))
SG = pd.DataFrame(segs).sort_values(['src', 'side', 'a0', 'a1']).reset_index(drop=True)
arms = []
for (src, side), s in SG.groupby(['src', 'side'], sort=True):
    k = s[(s.el == 'K') & (s.a0 > 200)]
    aend = np.nan
    if len(k):
        aend = k.a1.iloc[0]
        for q in s[s.a0 >= k.a0.iloc[0]].itertuples():
            if q.el == 'K' and q.a0 - aend < 15: aend = max(aend, q.a1)
            elif q.el not in ('K', '.'): break
    beyond = s[(s.a0 >= aend) & (s.el != '.')] if aend == aend else s.iloc[0:0]
    lay = []
    for e in beyond.el:
        if not lay or lay[-1] != e: lay.append(e)
    arms.append(dict(src=src, side=side, a_end_kb=aend, layout='A ' + ' '.join(lay), reach_kb=s.a1.max()))
A = pd.DataFrame(arms)
# element kb beyond A, per arm
SA = SG.merge(A[['src', 'side', 'a_end_kb']], on=['src', 'side'])
SA = SA[(SA.a0 >= SA.a_end_kb) & (SA.el != '.')].assign(kb=lambda x: x.a1 - x.a0)
kb = SA.groupby(['src', 'side', 'el']).kb.sum().unstack(fill_value=0)
for e in ['Xi', 'S', 'Xd', 'E', 'K', 'Fl']: A[f'kb_{e}'] = [round(kb.loc[(a, b), e], 1) if (a, b) in kb.index else 0.0 for a, b in zip(A.src, A.side)]
A['xin_kb'] = [round(SA[(SA.src == a) & (SA.side == b) & SA.tok.isin(XIN)].kb.sum(), 1) for a, b in zip(A.src, A.side)]
has = lambda e: A[f'kb_{e}'] > 0
A['type'] = np.where(A.side == 'cen', np.where(has('Xd'), 'C1', np.where(~has('Xi') & ~has('S'), 'C0', 'other')),
                     np.where(has('Xd'), 'T1', np.where(has('Xi') & has('S'), 'T0', 'other')))
CANON = {'C1': 'A Xi S Xd K E Fl', 'C0': 'A E Fl', 'T0': 'A E K Xi S', 'T1': 'A E Xi S Xd K E K Xi S'}   # under the clinical decomposition a flank stretch sat between E and K; here bundle 4 joins it to E
A['canonical'] = A.layout == A.type.map(CANON)
# checks on the element rule
c = A[A.side == 'cen']
assert ((c.kb_Xi > 0) == (c.kb_Xd > 0)).all(), 'a cenSMN arm with Xi but no Xd, or the reverse'
assert (A.type != 'other').all(), A[A.type == 'other']
assert A[A.kb_Xd > 0].kb_Xd.min() >= 20, 'a marginal Xd'
assert A.a_end_kb.notna().all(), 'an arm without an end-of-core block'
for t, lay in CANON.items(): print(t, (A.type == t).sum(), A[A.type == t].layout.value_counts().head(4).to_dict())
for t, lay in CANON.items(): assert A[A.type == t].layout.value_counts().index[0] == lay, ('canonical layout is not the commonest', t)
# copies (step 22): the counted SMN copy on each arm
CP = pd.read_csv(O + 'palindrome_copies_lineage.tsv', sep='\t', keep_default_na=False)
CP['side'] = np.where(CP.d_kb < 0, 'cen', 'tel'); CP['d'] = CP.d_kb.abs()
cnt = CP[~CP.token.str[0].isin(['t', 'f'])]
assert (cnt.groupby(['src', 'side']).size() == 1).all() and len(cnt) == 2 * len(AX)
cnt = cnt.set_index(['src', 'side'])
A['copy_d_kb'] = [cnt.loc[(a, b), 'd'] for a, b in zip(A.src, A.side)]; A['copy_mark'] = [cnt.loc[(a, b), 'mark'] for a, b in zip(A.src, A.side)]
assert (A.copy_d_kb < A.a_end_kb).all(), 'an SMN copy outside A'
# tandem arrays (step 28), units >= 10 kb, placed on arms; inner = starts inside A, outer = beyond A
T = pd.read_csv(O + 'smn.tandem_repeats.tsv', sep='\t', comment='#'); T = T[T.unit_size >= 10000]
ar = []
for r in T[T.row.isin(AX.index)].itertuples():
    ax = AX.axis[r.row]; m = (r.start + r.end) / 2
    side = 'cen' if m < ax else 'tel'
    d0, d1 = ((ax - r.end) / 1e3, (ax - r.start) / 1e3) if side == 'cen' else ((r.start - ax) / 1e3, (r.end - ax) / 1e3)
    ae = A[(A.src == r.row) & (A.side == side)].a_end_kb.iloc[0]
    ar.append(dict(src=r.row, side=side, d0=round(d0, 1), d1=round(d1, 1), unit_kb=round(r.unit_size / 1e3, 2), copies=r.copies,
                   where='inner' if d0 < ae else 'outer'))
AR = pd.DataFrame(ar).sort_values(['src', 'side', 'd0']); AR.to_csv(O + 'arm_element_arrays.tsv', sep='\t', index=False)
g = AR.groupby(['src', 'side', 'where']).copies.sum().unstack(fill_value=0)
for w in ('inner', 'outer'): A[f'{w}_copies'] = [round(g.loc[(a, b), w], 2) if (a, b) in g.index else 0.0 for a, b in zip(A.src, A.side)]
u = AR.groupby(['src', 'side', 'where']).unit_kb.first().unstack()
for w in ('inner', 'outer'): A[f'{w}_unit_kb'] = [u.loc[(a, b), w] if (a, b) in u.index else np.nan for a, b in zip(A.src, A.side)]
# genes in the palindrome window, placed on arms and elements
F = pd.read_csv('results/copies/features.tsv', sep='\t'); F = F[F.kind == 'gene']
win = CP.assign(w=CP.copy_id.str.extract(r'\|w(\d+)\|')[0].astype(int)).groupby('src').w.agg(lambda s: s.mode()[0])
gl = []
for src, ax in AX.axis.items():
    for r in F[(F.src == src) & (F.window == win[src])].itertuples():
        m = (r.win_start + r.win_end) / 2; side = 'cen' if m < ax else 'tel'; d = abs(m - ax) / 1e3
        ae = A[(A.src == src) & (A.side == side)].a_end_kb.iloc[0]
        s = SG[(SG.src == src) & (SG.side == side) & (SG.a0 >= ae)]
        if d <= ae: el = 'A'
        else:
            dist = np.maximum(s.a0 - d, d - s.a1).clip(lower=0)
            el = s.el[dist.idxmin()] if len(s) and dist.min() <= 10 else 'gap'
            el = 'other' if el == '.' else el
        gl.append(dict(src=src, side=side, gene=r.gene, d_kb=round(d, 1), tx_frac=r.tx_frac, whole=r.tx_frac >= 0.9, tx_from=r.tx_from, tx_to=r.tx_to, element=el))
GL = pd.DataFrame(gl).sort_values(['src', 'side', 'd_kb']); GL.to_csv(O + 'arm_element_genes.tsv', sep='\t', index=False)
for gn in ('NAIP', 'GTF2H2'):
    for wh in (True, False):
        q = GL[(GL.gene == gn) & (GL.whole == wh)].groupby(['src', 'side']).size()
        A[f'{gn}_{"whole" if wh else "part"}'] = [int(q.get((a, b), 0)) for a, b in zip(A.src, A.side)]
# S occurrences: runs of S-bundle segments less than 30 kb apart, anywhere in the window (window coordinates)
SS = []
for src, ax in AX.axis.items():
    q = bed[(bed.row == src) & bed.bid.isin(SBID)].sort_values('bgn')
    run = []
    for r in q.itertuples():
        if run and r.bgn - max(x.end for x in run) > 30000: SS.append((src, ax, run)); run = []
        run.append(r)
    if run: SS.append((src, ax, run))
srows = []
for src, ax, run in SS:
    b, e = run[0].bgn, max(x.end for x in run); m = (b + e) / 2
    order = '+' if run[0].d == 0 else '-'   # strand of S on the window
    site = 'axis' if (b <= ax <= e or abs(m - ax) < 30000) else ('cen' if m < ax else 'tel')
    srows.append(dict(src=src, site=site, bgn=b, end=e, span_kb=round((e - b) / 1e3, 1), d_kb=round((m - ax) / 1e3, 1),
                      bundle_kb=round(sum(x.end - x.bgn for x in run) / 1e3, 1), toks=' '.join(f'{x.bid}{"+-"[x.d]}' for x in run), window_order=order))
SO = pd.DataFrame(srows); SO.to_csv(O + 'arm_element_spacers.tsv', sep='\t', index=False)
# haplotype table
CL = pd.read_csv(O + 'copy_labels.tsv', sep='\t').dropna(subset=['superpop']).groupby('src').superpop.first()
H = A.pivot(index='src', columns='side', values='type').rename(columns={'cen': 'cen_type', 'tel': 'tel_type'})
H['cen_mark'] = A[A.side == 'cen'].set_index('src').copy_mark; H['tel_mark'] = A[A.side == 'tel'].set_index('src').copy_mark
H['superpop'] = H.index.map(CL); assert H.superpop.notna().all()
H['afr'] = H.superpop == 'AFR'
H['b1b2'] = ((H.cen_mark == 'SMN1 B1') & (H.tel_mark == 'SMN1 B2')) | ((H.cen_mark == 'SMN1 B2') & (H.tel_mark == 'SMN1 B1'))
H['axis_S'] = [int(((SO.src == s_) & (SO.site == 'axis')).sum()) for s_ in H.index]
H.to_csv(O + 'arm_element_haps.tsv', sep='\t')
A.to_csv(O + 'arm_elements.tsv', sep='\t', index=False)
SG.to_csv(O + 'arm_element_segments.tsv', sep='\t', index=False)
print(pd.crosstab(A.side, A.type)); print(A.groupby('type').canonical.agg(['sum', 'size']))
print(SO.site.value_counts().to_dict(), 'axis S per hap', H.axis_S.value_counts().to_dict())
# ---- NAIP pieces against A's end and X (spans, not midpoints)
NP = []
for src, ax in AX.axis.items():
    for r in F[(F.src == src) & (F.window == win[src]) & (F.gene == 'NAIP')].itertuples():
        m = (r.win_start + r.win_end) / 2; side = 'cen' if m < ax else 'tel'
        d0, d1 = ((ax - r.win_end) / 1e3, (ax - r.win_start) / 1e3) if side == 'cen' else ((r.win_start - ax) / 1e3, (r.win_end - ax) / 1e3)
        a = A[(A.src == src) & (A.side == side)].iloc[0]
        x = SG[(SG.src == src) & (SG.side == side) & SG.el.isin(['Xi', 'S', 'Xd']) & (SG.a0 >= a.a_end_kb)]
        piece = 'whole' if r.tx_from == 0 and r.tx_to >= 7000 else '5p' if r.tx_from == 0 else '3p' if r.tx_to >= 7000 and r.tx_from < 3000 else 'fragment'
        NP.append(dict(src=src, side=side, type=a.type, tx_from=r.tx_from, tx_to=r.tx_to, piece=piece, d0=round(d0, 2), d1=round(d1, 2), a_end_kb=a.a_end_kb,
                       next_x_kb=round((x.a0 - d1).abs().min(), 2) if len(x) else np.nan, spans_a_end=bool(d0 < a.a_end_kb < d1), ends_at_a_end=bool(abs(d1 - a.a_end_kb) < 5)))
NP = pd.DataFrame(NP).sort_values(['src', 'side', 'd0']); NP.to_csv(O + 'arm_naip_pieces.tsv', sep='\t', index=False)
# ---- summary
from scipy.stats import fisher_exact, spearmanr, mannwhitneyu
J = {}
J['n_haps'] = len(H); J['types'] = {t: int(n) for t, n in A.type.value_counts().items()}
J['canonical'] = {t: int(n) for t, n in A[A.canonical].type.value_counts().items()}
J['canon_layout'] = CANON
ct = pd.crosstab(H.cen_type, H.tel_type); J['cen_x_tel'] = {f'{a}{b}': int(ct.loc[a, b]) for a in ct.index for b in ct.columns}
J['cen_x_tel_p'] = fisher_exact(ct.values)[1]
for nm, q in (('afr', H[H.afr]), ('nonafr', H[~H.afr])):
    c2 = pd.crosstab(q.cen_type, q.tel_type).reindex(index=['C0', 'C1'], columns=['T0', 'T1'], fill_value=0)
    J[f'cen_x_tel_{nm}'] = {f'{a}{b}': int(c2.loc[a, b]) for a in c2.index for b in c2.columns}; J[f'cen_x_tel_{nm}_p'] = fisher_exact(c2.values)[1]
J['C0_by_superpop'] = {g: [int((q.cen_type == 'C0').sum()), len(q)] for g, q in H.groupby('superpop')}
J['T1_by_superpop'] = {g: [int((q.tel_type == 'T1').sum()), len(q)] for g, q in H.groupby('superpop')}
b1 = H.cen_mark == 'SMN1 B1'; a_ = H[H.afr]
J['B1cen_C0'] = [int((b1 & (H.cen_type == 'C0')).sum()), int(b1.sum())]
t = pd.crosstab(a_.cen_mark == 'SMN1 B1', a_.cen_type == 'C0').reindex(index=[True, False], columns=[True, False], fill_value=0)
J['afr_B1cen_C0'] = [int(t.loc[True, True]), int(t.loc[True].sum())]; J['afr_otherCen_C0'] = [int(t.loc[False, True]), int(t.loc[False].sum())]; J['afr_B1cen_C0_p'] = fisher_exact(t.values)[1]
b2 = H.tel_mark == 'SMN1 B2'
J['B2tel_T1'] = [int((b2 & (H.tel_type == 'T1')).sum()), int(b2.sum())]
t = pd.crosstab(a_.tel_mark == 'SMN1 B2', a_.tel_type == 'T1').reindex(index=[True, False], columns=[True, False], fill_value=0)
J['afr_B2tel_T1'] = [int(t.loc[True, True]), int(t.loc[True].sum())]; J['afr_otherTel_T1'] = [int(t.loc[False, True]), int(t.loc[False].sum())]; J['afr_B2tel_T1_p'] = fisher_exact(t.values)[1]
bb = H[H.b1b2]; J['b1b2'] = dict(n=len(bb), C0=int((bb.cen_type == 'C0').sum()), T1=int((bb.tel_type == 'T1').sum()), C0T1=int(((bb.cen_type == 'C0') & (bb.tel_type == 'T1')).sum()))
nb = H[~H.b1b2]; J['C0_not_b1b2'] = [int((nb.cen_type == 'C0').sum()), int((nb[nb.cen_type == 'C0'].afr).sum())]
for side in ('cen', 'tel'):
    q = A[A.side == side]
    J[f'{side}_copy_d_by_type'] = {t: float(g.copy_d_kb.median()) for t, g in q.groupby('type')}
    J[f'{side}_copy_d_p'] = mannwhitneyu(*[g.copy_d_kb for _, g in q.groupby('type')])[1]
    J[f'{side}_aend_by_type'] = {t: float(g.a_end_kb.median()) for t, g in q.groupby('type')}
    J[f'{side}_inner_by_type'] = {t: float(g.inner_copies.median()) for t, g in q.groupby('type')}
    ok = q.inner_copies > 0
    J[f'{side}_rho_copy_inner'] = float(spearmanr(q.copy_d_kb[ok], q.inner_copies[ok])[0]); J[f'{side}_n_inner'] = int(ok.sum())
    J[f'{side}_copy_d_by_inner'] = {int(k): float(v) for k, v in q[ok].groupby(q[ok].inner_copies.round(0)).copy_d_kb.median().items()}
    J[f'{side}_inner_unit_kb'] = float(q.inner_unit_kb.median())
    J[f'{side}_outer_by_type'] = {t: dict(with_array=int((g.outer_copies > 0).sum()), n=len(g), copies_median=float(g[g.outer_copies > 0].outer_copies.median()) if (g.outer_copies > 0).any() else 0.0,
                                          unit_kb=float(g.outer_unit_kb.median()) if g.outer_unit_kb.notna().any() else 0.0) for t, g in q.groupby('type')}
    J[f'{side}_element_kb'] = {t: {e: float(g[g[f'kb_{e}'] > 0][f'kb_{e}'].median()) if (g[f'kb_{e}'] > 0).any() else 0.0 for e in ('Xi', 'S', 'Xd', 'E')} for t, g in q.groupby('type')}
    J[f'{side}_reach_min'] = float(q.reach_kb.min())
t2 = A[A.type == 'T1'].outer_copies; t1 = A[A.type == 'T0'].outer_copies; J['tel_outer_T0_T1_p'] = mannwhitneyu(t1[t1 > 0], t2[t2 > 0])[1]
# NAIP
w = NP[NP.piece == 'whole']
J['naip_whole_spans_aend'] = {t: [int(g.spans_a_end.sum()), int(A[A.type == t].shape[0])] for t, g in w.groupby('type')}
J['naip_whole_past_aend_kb'] = float((w[w.spans_a_end].d1 - w[w.spans_a_end].a_end_kb).median())
c1 = NP[(NP.type == 'C1') & (NP.piece == '3p') & NP.ends_at_a_end]
J['C1_3p_at_aend'] = [int(c1.src.nunique()), int((A.type == 'C1').sum())]; J['C1_3p_txfrom'] = {int(k): int(v) for k, v in c1.tx_from.value_counts().items()}
J['C1_3p_next_x_kb'] = {int(k): float(g.next_x_kb.median()) for k, g in c1.groupby('tx_from')}
c5 = NP[(NP.type == 'C1') & (NP.piece == '5p')]
J['C1_5p'] = [int(c5.src.nunique()), int(c5.tx_to.median())]
J['C1_whole_cen'] = int(NP[(NP.type == 'C1') & (NP.piece == 'whole')].src.nunique())
t13 = NP[(NP.type == 'T0') & (NP.piece == '3p')]
J['T0_3p_before_x'] = [int((t13.next_x_kb < 5).sum()), int((A.type == 'T0').sum())]
t23 = NP[(NP.type == 'T1') & (NP.piece == '3p')]; t25 = NP[(NP.type == 'T1') & (NP.piece == '5p')]
J['T1_3p_before_x'] = [int((t23.next_x_kb < 5).sum()), int((A.type == 'T1').sum())]; J['T1_5p'] = int(t25.src.nunique())
J['naip_whole_per_hap'] = {k: int(v) for k, v in NP[NP.piece == 'whole'].groupby('src').size().reindex(H.index, fill_value=0).groupby(H.cen_type).apply(lambda s: (s == 2).sum()).items()}
J['naip_whole2_C0'] = [int((NP[NP.piece == 'whole'].groupby('src').size().reindex(H.index, fill_value=0)[H.cen_type == 'C0'] == 2).sum()), int((H.cen_type == 'C0').sum())]
J['naip_whole1_C1'] = [int((NP[NP.piece == 'whole'].groupby('src').size().reindex(H.index, fill_value=0)[H.cen_type == 'C1'] == 1).sum()), int((H.cen_type == 'C1').sum())]
# E carries GTF2H2
gt = GL[GL.gene == 'GTF2H2']
J['gtf2h2_in_E'] = [int((gt.element == 'E').sum()), len(gt)]; J['gtf2h2_whole_in_E'] = [int(((gt.element == 'E') & gt.whole).sum()), int((gt.element == 'E').sum())]
J['E_arms'] = {s_: int((A[A.side == s_].kb_E > 0).sum()) for s_ in ('cen', 'tel')}
J['axis_S'] = int((H.axis_S == 1).sum())
J['S_per_window'] = {int(k): int(v) for k, v in SO.groupby('src').size().value_counts().sort_index().items()}
# the earlier long/short call (step 27: the family-A unit of the arm is 500 kb or more) against the arm types
U27 = pd.read_csv(O + 'arm_units.tsv', sep='\t').assign(side=lambda d: d.side.map({'left': 'cen', 'right': 'tel'}))
UL = A.merge(U27[['src', 'side', 'span_kb']], on=['src', 'side'], validate='one_to_one')
J['unit_long_by_type'] = {t: [int((g.span_kb >= 500).sum()), len(g)] for t, g in UL.groupby('type')}
json.dump(J, open(O + 'arm_elements_summary.json', 'w'), indent=1, default=float)
print(json.dumps(J, indent=1, default=float))
# checks for the statements the note makes
assert J['b1b2']['C0'] == J['b1b2']['n'], 'a B1+B2 haplotype carries X on its cenSMN arm'
_c1 = set(c1.src) | set(NP[(NP.type == 'C1') & (NP.piece == 'whole') & NP.spans_a_end].src)
assert len(_c1) == J['C1_3p_at_aend'][1] and not (set(c1.src) & set(NP[(NP.type == 'C1') & (NP.piece == 'whole') & NP.spans_a_end].src)), 'a C1 arm whose NAIP is neither split at A end nor whole across it'
assert J['C1_3p_at_aend'][0] >= 0.95 * J['C1_3p_at_aend'][1], 'text: NAIP split at A end on nearly every C1 arm'
assert all(v < 5 for v in J['C1_3p_next_x_kb'].values()), 'X does not start right after the NAIP 3-prime piece'
assert J['naip_whole_spans_aend']['C0'][0] == J['naip_whole_spans_aend']['C0'][1], 'a C0 arm without a whole NAIP across A end'
assert J['cen_outer_by_type']['C0']['with_array'] == 0, 'an outer array on an arm without X'
assert J['gtf2h2_whole_in_E'][0] == J['gtf2h2_whole_in_E'][1]
