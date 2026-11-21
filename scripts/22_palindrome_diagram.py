# Pipeline step 22 (palindrome_diagram); ported from explore/structure_vs_sequence/19_palindrome_diagram.py (sha256 141fa992f412, commit 5e2b3bb),
# paths rewritten to results/sv/ and figures/; run from the project root after the steps before it.
"""Exploration step 19: diagram of the family-A palindrome (Section 4.11).
(a) one representative haplotype drawn to scale: its repeat units, the adjacent head-to-head A pair, the axis, the SMN copies;
(b) all haplotypes with an adjacent A+ A- pair, aligned on the axis, arms as bars and copies as marks;
(c) each haplotype's nearest left copy against nearest right copy, distance from the axis (mirror = on the diagonal).
Reads results/decomp/smn.repeat_units.tsv, results/decomp/copy_units.tsv, tables/haplotype_copy_counts.tsv,
out/palindrome_axis_copies.tsv (step 10, used as a check), doc/marker_lineage.tsv (SMN1 lineage A/B, by copy_id).
Also reads results/decomp/smn.bed (nested repeats, step 24's rule).
Writes out/palindrome_diagram_rows.tsv, out/palindrome_copies_lineage.tsv, out/palindrome_arm_marks.tsv (the figure moved to step 28d)."""
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, Rectangle
# (exploration folder variable removed)
RU = pd.read_csv('results/decomp/smn.repeat_units.tsv', sep='\t', comment='#')
CU = pd.read_csv('results/decomp/copy_units.tsv', sep='\t')
hc = pd.read_csv('tables/haplotype_copy_counts.tsv', sep='\t')
res = set(hc[hc.resolved].src)
OA10 = pd.read_csv('results/sv/palindrome_axis_copies.tsv', sep='\t')
rows, cps = [], []
for src, g in RU.groupby('row'):
    if src not in res: continue
    g = g.sort_values('bgn').reset_index(drop=True)
    for i in range(len(g) - 1):
        a, b = g.iloc[i], g.iloc[i + 1]
        if a.family == 'A' and b.family == 'A' and a.strand == '+' and b.strand == '-':
            axis = (a.end + b.bgn) / 2
            rows.append(dict(src=src, axis=axis, lb=(a.bgn - axis) / 1e3, le=(a.end - axis) / 1e3, rb=(b.bgn - axis) / 1e3, re=(b.end - axis) / 1e3))
            for c in CU[CU.src == src].itertuples():
                if a.bgn <= c.mid_fwd <= b.end: cps.append(dict(src=src, copy_id=c.copy_id, token=c.token, d_kb=(c.mid_fwd - axis) / 1e3))
R = pd.DataFrame(rows); C = pd.DataFrame(cps)
# check: same copies and distances as step 10
assert len(C) == len(OA10) and np.allclose(np.sort(C.d_kb.values), np.sort(OA10.d_kb.values)), 'differs from step 10'
assert R.src.is_unique
# a head-to-head pair counts as the palindrome only when each arm carries an SMN copy; the others are listed
_side = C.assign(side=np.sign(C.d_kb)).groupby('src').side.agg(lambda x: set(x))
_keep = {s for s, v in _side.items() if {-1, 1} <= v}
R[~R.src.isin(_keep)].assign(copies_kb=lambda x: [' '.join(f'{d:.1f}' for d in sorted(C[C.src == s].d_kb)) for s in x.src]).to_csv('results/sv/palindrome_one_arm.tsv', sep='\t', index=False)
print('head-to-head pairs', len(R), 'with a copy on each arm', len(_keep), 'left out:', sorted(set(R.src) - _keep))
R = R[R.src.isin(_keep)].reset_index(drop=True); C = C[C.src.isin(_keep)].reset_index(drop=True)
def gene(t):
    return 'SMN1' if t[0] == '1' else 'SMN2' if t[0] == '2' else 'hybrid' if t[0] == 'h' else 'other'
import sys; sys.path.insert(0, 'scripts'); from smnlib import COPYCOL, NESTCOL
COL = {'SMN1 A': COPYCOL['A'], 'SMN1 B1': COPYCOL['B1'], 'SMN1 B2': COPYCOL['B2'], 'SMN1 B-': COPYCOL['B-'], 'SMN1 B': COPYCOL['B1'], 'SMN1 not in tree': '#fb9a99', 'SMN1': COPYCOL['A'], 'SMN2': COPYCOL['SMN2'], 'hybrid': COPYCOL['hybrid'], 'other': COPYCOL['other']}
C['gene'] = C.token.map(gene)
# SMN1 lineage A/B (Section 4.4): B = in_marker_subclade among SMN1 tips of the copy tree
ML = pd.read_csv('results/sv/marker_lineage.tsv', sep='\t'); ML = ML[ML['class'] == 'SMN1']
lin = dict(zip(ML.copy_id, np.where(ML.in_marker_subclade, 'B', 'A'))); sublin = dict(zip(ML.copy_id, ML.sublineage))
C['lineage'] = [lin.get(i, 'not in tree') if g_ == 'SMN1' else '' for i, g_ in zip(C.copy_id, C.gene)]
C['arm'] = np.where(C.d_kb < 0, 'left (A+, centromeric)', 'right (A-, telomeric)')
C['sublineage'] = [sublin.get(i, '') if g_ == 'SMN1' else '' for i, g_ in zip(C.copy_id, C.gene)]
C['mark'] = np.where(C.gene == 'SMN1', 'SMN1 ' + np.where(C.sublineage != '', C.sublineage, C.lineage), C.gene)
C.to_csv('results/sv/palindrome_copies_lineage.tsv', sep='\t', index=False)
print(pd.crosstab(C[C.gene == 'SMN1'].lineage, C[C.gene == 'SMN1'].arm))
# nearest copy on each arm
L = C[C.d_kb < 0].assign(a=lambda x: -x.d_kb).groupby('src').a.min(); Rr = C[C.d_kb > 0].groupby('src').d_kb.min()
R = R.set_index('src').join(L.rename('near_left')).join(Rr.rename('near_right')).reset_index()
assert R.near_left.notna().all() and R.near_right.notna().all(), 'a pair haplotype lacks a copy on one arm'
R['mirror_kb'] = (R.near_left - R.near_right).abs()
R = R.sort_values(['lb', 'src']).reset_index(drop=True)
# checks guarding the Fig. caption: 'other' = 7 delta7-8 copies + 1 truncated piece; arms vary more than twofold while
# every copy stays within the axis-distance range the text states (palDmin-palDmax)
assert C[C.gene == 'other'].token.value_counts().to_dict() == {'d2-': 7, 't-': 1}
for a_ in [R['le'] - R['lb'], R['re'] - R['rb']]: assert a_.max() / a_.min() > 2
print('copy distance from the axis, kb: min', round(C.d_kb.abs().min(), 1), 'max', round(C.d_kb.abs().max(), 1))   # the range the text states is checked in step 40
R.to_csv('results/sv/palindrome_diagram_rows.tsv', sep='\t', index=False)
print('pair haplotypes', len(R), 'copies', len(C), 'median mirror kb', R.mirror_kb.median(), 'genes', C.gene.value_counts().to_dict())

import sys; sys.path.insert(0, 'scripts'); from smnlib import FAMCOL as fam_col   # same family colours as Figs. 1 and 2
# continental group of each haplotype (copy_labels.tsv, step 04); AFR against the rest, as in step 27's arm_slots.tsv
SP = pd.read_csv('results/sv/copy_labels.tsv', sep='\t').dropna(subset=['superpop']).groupby('src').superpop.first()
assert R.src.map(SP).notna().all() and C.src.map(SP).notna().all()
R['grp'] = np.where(R.src.map(SP) == 'AFR', 'AFR', 'other'); C['grp'] = np.where(C.src.map(SP) == 'AFR', 'AFR', 'other')
C['side'] = np.where(C.d_kb < 0, 'left', 'right')
ORDER = ['SMN1 B1', 'SMN1 B2', 'SMN1 B-', 'SMN1 A', 'SMN2', 'hybrid', 'other']
assert set(C['mark']) <= set(ORDER), set(C['mark'])
T = pd.crosstab([C.side, C.grp], C['mark']).reindex(columns=ORDER, fill_value=0)
T.to_csv('results/sv/palindrome_arm_marks.tsv', sep='\t')
print(T)
# The figure (Fig. 4) and the nested-repeat tables this step used to write are superseded by steps 28b and 28d (2026-11-14):
# the arms are now read as named elements (scripts/28b_arm_elements.py) and drawn by scripts/28d_fig4.py.
