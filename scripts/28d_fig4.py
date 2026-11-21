"""Pipeline step 28d: Fig. 4, the family-A palindrome read as named elements (replaces step 22's drawing).
(a) what each arm holds, AFR against the other groups (results/sv/palindrome_arm_marks.tsv, step 22);
(b) the four haplotype forms (cenSMN type C1/C0 x telSMN type T0/T1, step 28b), each drawn from the median block
    positions of the arms with that type's canonical layout;
(c) every head-to-head haplotype's whole window aligned on its axis, coloured by element (step 28b segments: A from the
    axis to the end of the core, then Xi, S, Xd, E, flank; pieces of one element joined across bundle-free gaps,
    'other' pieces under 5 kb skipped), S also at the main axis, the ~20 kb tandem arrays (step 28), SMN copies
    (step 22) and NAIP matches (step 28b; filled = whole, open = partial);
(d) nearest copy on each arm, kb from the axis (mirror positions).
Writes figures/palindrome.png/.pdf and results/sv/arm_element_blocks.tsv."""
import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Patch; from matplotlib.lines import Line2D
sys.path.insert(0, 'scripts'); from smnlib import COPYCOL, ELCOL, ELLAB, elcol, eltxt, xi_part_labels, XI_PRIME
O = 'results/sv/'
A = pd.read_csv(O + 'arm_elements.tsv', sep='\t'); SG = pd.read_csv(O + 'arm_element_segments.tsv', sep='\t')
H = pd.read_csv(O + 'arm_element_haps.tsv', sep='\t').set_index('src'); AR = pd.read_csv(O + 'arm_element_arrays.tsv', sep='\t')
NP = pd.read_csv(O + 'arm_naip_pieces.tsv', sep='\t'); SO = pd.read_csv(O + 'arm_element_spacers.tsv', sep='\t')
CP = pd.read_csv(O + 'palindrome_copies_lineage.tsv', sep='\t', keep_default_na=False)
D = pd.read_csv(O + 'palindrome_diagram_rows.tsv', sep='\t').set_index('src')
T = pd.read_csv(O + 'palindrome_arm_marks.tsv', sep='\t', index_col=[0, 1])
CC = pd.read_csv('results/copies/copy_calls.tsv', sep='\t', keep_default_na=False).set_index('copy_id').copy_class
# ---- blocks per arm
BL = []
for (src, side), sg in SG.groupby(['src', 'side'], sort=True):
    ae = A[(A.src == src) & (A.side == side)].a_end_kb.iloc[0]
    blk = [[0.0, ae, 'A']]
    for q in sg[sg.a1 > ae].sort_values('a0').itertuples():
        e = {'.': 'other', 'K': 'A'}.get(q.el, q.el)
        if e == 'other' and q.a1 - q.a0 < 5: continue
        if blk[-1][2] == e: blk[-1][1] = max(blk[-1][1], q.a1)
        else: blk.append([max(q.a0, ae), q.a1, e])
    for i, (a0, a1, e) in enumerate(blk): BL.append(dict(src=src, side=side, i=i, a0=round(a0, 2), a1=round(a1, 2), el=e))
BL = pd.DataFrame(BL); BL.to_csv(O + 'arm_element_blocks.tsv', sep='\t', index=False)
BL['lab'] = xi_part_labels(BL)   # display only: Xi' (inner part) or Xi (outer part)
# ---- copy marks
GC = {'SMN1 A': COPYCOL['A'], 'SMN1 B1': COPYCOL['B1'], 'SMN1 B2': COPYCOL['B2'], 'SMN1 B-': COPYCOL['B-'], 'SMN2': COPYCOL['SMN2'],
      'hybrid (h12 / h21)': COPYCOL['hybrid'], 'Δ7–8 copy': '#ff7f00', 'truncated piece': '#bdbdbd'}
def cmark(c):
    if c.mark in GC: return c.mark
    k = CC.get(c.copy_id, '')
    if k.startswith('hyb'): return 'hybrid (h12 / h21)'
    if k.startswith('del7_8'): return 'Δ7–8 copy'
    if k in ('truncated', 'fragment'): return 'truncated piece'
    raise SystemExit(f'unclassified copy {c.copy_id} {c.mark} {k}')
CP['cm'] = [cmark(c) for c in CP.itertuples()]
CP.cm.value_counts().rename_axis('mark').rename('items').to_csv('results/sv/palindrome_dot_classes.tsv', sep='\t')   # Fig. 4c dot classes (caption)
# ---- rows: grouped by the haplotype form, AFR first within a group, then by cenSMN reach
FORMS = [('C1', 'T0'), ('C1', 'T1'), ('C0', 'T0'), ('C0', 'T1')]
reach = A[A.side == 'cen'].set_index('src').reach_kb
H['form'] = list(zip(H.cen_type, H.tel_type)); H['fi'] = H.form.map({f: i for i, f in enumerate(FORMS)})
assert H.fi.notna().all()
R = H.assign(rc=reach, na=~H.afr).sort_values(['fi', 'na', 'rc'], ascending=[False, True, True]).reset_index()
fig = plt.figure(figsize=(12, 15))
gs = fig.add_gridspec(3, 2, height_ratios=[0.8, 0.9, 3.3], width_ratios=[2.5, 1], hspace=0.28, wspace=0.32)
# (a)
ax = fig.add_subplot(gs[0, :])
ORDER = ['SMN1 B1', 'SMN1 B2', 'SMN1 B-', 'SMN1 A', 'SMN2', 'hybrid', 'other']
COLA = {'SMN1 A': COPYCOL['A'], 'SMN1 B1': COPYCOL['B1'], 'SMN1 B2': COPYCOL['B2'], 'SMN1 B-': COPYCOL['B-'], 'SMN2': COPYCOL['SMN2'], 'hybrid': COPYCOL['hybrid'], 'other': COPYCOL['other']}
NAME = {'left': 'cenSMN arm', 'right': 'telSMN arm'}
for side, grp, y in [('left', 'AFR', 3.6), ('left', 'other', 2.6), ('right', 'AFR', 1.0), ('right', 'other', 0.0)]:
    row = T.loc[(side, grp)]; n = int(row.sum()); x = 0.0
    for m in ORDER:
        v = int(row.get(m, 0)); w = v / n
        if v == 0: continue
        ax.barh(y, w, left=x, height=0.8, color=COLA[m], edgecolor='white', lw=0.6)
        if w >= 0.045: ax.text(x + w / 2, y, str(v), ha='center', va='center', fontsize=8, color='white' if m != 'other' else 'black')
        x += w
    nh = int((H.afr == (grp == 'AFR')).sum())
    ax.text(-0.01, y, f'{NAME[side]}, {"AFR" if grp == "AFR" else "other groups"}\n{n} items, {nh} haplotypes', ha='right', va='center', fontsize=8.5)
ax.set_xlim(0, 1); ax.set_ylim(-0.6, 4.2); ax.set_yticks([]); ax.set_xlabel('share of the items on the arm', fontsize=9)
for s_ in ['top', 'right', 'left']: ax.spines[s_].set_visible(False)
ax.legend(handles=[Patch(color=COLA[m], label=m.replace('SMN1 B-', 'SMN1 B, neither clade')) for m in ORDER if T[m].sum()], loc='upper center', bbox_to_anchor=(0.5, -0.3), ncol=7, fontsize=8, frameon=False)
ax.set_title('a  What each arm holds, AFR haplotypes against the other groups', loc='left', fontsize=10.5)
# (b) schematic of the four forms
ax = fig.add_subplot(gs[1, :])
CAN = A[A.canonical]
sch = {}
for t, g in CAN.groupby('type'):
    b = BL.merge(g[['src', 'side']], on=['src', 'side'])
    sch[t] = b.groupby('i').agg(a0=('a0', 'median'), a1=('a1', 'median'), el=('el', 'first'), lab=('lab', lambda x: x.mode().iloc[0])).reset_index()
for k, (ct, tt) in enumerate(FORMS):
    y = 3 - k; n = int(((H.cen_type == ct) & (H.tel_type == tt)).sum()); na = int(((H.cen_type == ct) & (H.tel_type == tt) & H.afr).sum())
    for t, sgn in ((ct, -1), (tt, 1)):
        for q in sch[t].itertuples():
            ax.barh(y, sgn * (q.a1 - q.a0), left=sgn * q.a0, height=0.55, color=elcol(q.el, q.lab), lw=0)
            if q.a1 - q.a0 >= 30 and q.el != 'A': ax.text(sgn * (q.a0 + q.a1) / 2, y, q.lab, ha='center', va='center', fontsize=7.5, color=eltxt(q.lab))
    ax.barh(y, 21.3, left=-10.65, height=0.55, color=ELCOL['S'], lw=0)
    for t, sgn in ((ct, -1), (tt, 1)):   # bracket over X = Xi S Xd (C1 cenSMN arms, T1 telSMN arms)
        if t not in ('C1', 'T1'): continue
        q = sch[t]; s0 = q[q.el == 'S'].a0.min(); xi = q[(q.lab == 'Xi') & (q.a1 <= s0 + 1)]; xd = q[q.el == 'Xd']
        assert len(xi) and len(xd), t
        x0, x1 = sgn * xi.a0.max(), sgn * xd.a1.max(); yb = y + 0.36   # X = Xi (the outer part, next to S), S, Xd
        ax.plot([x0, x0, x1, x1], [yb - 0.06, yb, yb, yb - 0.06], color='k', lw=0.8)
        ax.text((x0 + x1) / 2, yb + 0.03, 'X (Xi–S–Xd)' if sgn < 0 else 'X, inverted copy (Xi–S–Xd)', ha='center', va='bottom', fontsize=7.5)
    ax.text(-1480, y + 0.12, f'{ct} | {tt}', ha='left', va='center', fontsize=10, weight='bold')
    ax.text(-1480, y - 0.2, f'{n} haplotypes ({na} AFR)', ha='left', va='center', fontsize=7.5)
ax.axvline(0, color='k', ls='--', lw=0.6); ax.set_xlim(-1500, 1250); ax.set_ylim(-0.6, 4.05); ax.set_yticks([])
ax.set_xlabel('kb from the axis (median block positions of the arms with the type\'s commonest layout)', fontsize=9)
for s_ in ['top', 'right', 'left']: ax.spines[s_].set_visible(False)
ax.text(-1000, 3.75, '← cenSMN arm', ha='center', fontsize=9); ax.text(900, 3.75, 'telSMN arm →', ha='center', fontsize=9)
ax.set_title('b  The four forms of the palindrome: cenSMN arm C1 (Xi′, then X = Xi S Xd) or C0 (without), telSMN arm T0 (Xi′ Xi S) or T1 (Xi S Xd … Xi′ Xi S)', loc='left', fontsize=10.5)
# (c) every haplotype
ax = fig.add_subplot(gs[2, 0]); h = 0.9
idx = {s: i for i, s in enumerate(R.src)}
for r in BL.itertuples():
    y = idx[r.src]; sgn = -1 if r.side == 'cen' else 1
    ax.barh(y, sgn * (r.a1 - r.a0), left=sgn * r.a0, height=h, color=elcol(r.el, r.lab), lw=0, zorder=1)
for r in SO[SO.site == 'axis'].itertuples():
    ax.barh(idx[r.src], r.span_kb, left=r.d_kb - r.span_kb / 2, height=h, color=ELCOL['S'], lw=0, zorder=1.5)
for r in AR.itertuples():
    sgn = -1 if r.side == 'cen' else 1; ax.plot([sgn * r.d0, sgn * r.d1], [idx[r.src]] * 2, color='#333', lw=0.6, zorder=2)
for r in NP[NP.piece.isin(['whole', '3p', '5p'])].itertuples():
    x = (-1 if r.side == 'cen' else 1) * (r.d0 + r.d1) / 2
    ax.plot(x, idx[r.src], '^', ms=3.0, color='#2ca25f' if r.piece == 'whole' else 'none', mec='#2ca25f', mew=0.5, zorder=3)
for c in CP.itertuples():
    ax.plot(c.d_kb, idx[c.src], 'o', ms=3.2, color=GC[c.cm], mec='white', mew=0.3, zorder=4)
ax.axvline(0, color='k', ls='--', lw=0.6)
for k, f in enumerate(FORMS):
    ys = [idx[s] for s in R[R.form == f].src]
    ax.axhline(min(ys) - 0.5, color='k', lw=0.7)
    ax.text(1.005, (np.mean(ys) + 0.5) / len(R), f'{f[0]} | {f[1]}\n({len(ys)})', transform=ax.transAxes, fontsize=9, va='center', weight='bold')
for s, i in idx.items():
    if H.afr[s]: ax.plot(-0.008, i, 's', ms=1.6, color='#7f3b08', transform=ax.get_yaxis_transform(), clip_on=False)
ax.set_ylim(-0.5, len(R) - 0.5); ax.set_xlim(-A[A.side == 'cen'].reach_kb.max() - 20, A[A.side == 'tel'].reach_kb.max() + 20)
ax.set_xlabel('kb from the axis (cenSMN arm left, telSMN arm right), whole window', fontsize=9)
ax.set_ylabel(f'{len(R)} haplotypes with the head-to-head R$_A$ pair (brown square = AFR)', fontsize=9)
ax.set_title('c  Every pair haplotype, by element', loc='left', fontsize=10.5)
hd = [Patch(color=ELCOL[e], label='S, spacer at the axis and at X\'s center' if e == 'S' else ELLAB[e]) for e in ['A', XI_PRIME, 'Xi', 'S', 'Xd', 'E', 'Fl', 'other']] + [Line2D([], [], color='#333', lw=1.5, label='~20 kb tandem array')]
hd += [Line2D([], [], marker='o', ls='', color=c, label=g) for g, c in GC.items() if (CP.cm == g).any()]
hd += [Line2D([], [], marker='^', ls='', color='#2ca25f', label='NAIP, whole'), Line2D([], [], marker='^', ls='', mfc='none', mec='#2ca25f', label='NAIP, 5′ or 3′ piece')]
ax.legend(handles=hd, loc='upper center', bbox_to_anchor=(0.62, -0.06), ncol=3, fontsize=7.6, frameon=False)
# (d) mirror
ax = fig.add_subplot(gs[2, 1])
P = D.assign(xl=D.near_left.round(), xr=D.near_right.round()).groupby(['xl', 'xr']).size().reset_index(name='n')
assert P.n.sum() == len(D)
ax.scatter(P.xl, P.xr, s=6 + 4 * P.n, color='#444', alpha=0.55, lw=0)
for n_ in [1, 10, 30]: ax.scatter([], [], s=6 + 4 * n_, color='#444', alpha=0.55, lw=0, label=f'{n_} haplotype' + ('s' if n_ > 1 else ''))
ax.legend(loc='upper left', fontsize=7, frameon=False, labelspacing=1.2)
m0, m1 = min(D.near_left.min(), D.near_right.min()) - 10, max(D.near_left.max(), D.near_right.max()) + 10
ax.plot([m0, m1], [m0, m1], color='grey', lw=0.8, ls=':'); ax.set_xlim(m0, m1); ax.set_ylim(m0, m1); ax.set_aspect('equal')
ax.set_xlabel('nearest copy, cenSMN arm (kb from axis)', fontsize=9); ax.set_ylabel('nearest copy, telSMN arm (kb from axis)', fontsize=9)
ax.set_title(f'd  Mirror positions\n(median difference {D.mirror_kb.median():.1f} kb)', loc='left', fontsize=10.5)
fig.savefig('figures/palindrome.png', dpi=200, bbox_inches='tight')
fig.savefig('figures/palindrome.pdf', bbox_inches='tight', metadata={'CreationDate': None})
print('rows', len(R), {f'{a}{b}': int((R.form == (a, b)).sum()) for a, b in FORMS})
