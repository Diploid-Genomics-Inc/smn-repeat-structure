"""Pipeline step 31: the locus as the study defines it (Fig. 1). Two real haplotypes drawn to scale around the axis of
their head-to-head pair of family-A repeat units (R_A+ then R_A-), with the end of the core marked on each arm, one of each of the two commonest forms of step 28b:
(a) cenSMN arm C1 and telSMN arm T0 (X on the cenSMN arm), the two-copy '2- 1+' haplotype with both arms in their
type's commonest layout whose arm lengths and mirror difference are nearest the medians (step 22's rule);
(b) cenSMN arm C0 and telSMN arm T1, the same rule among '1- 1+' haplotypes carrying one B1 and one B2 copy.
On each: repeat units (results/decomp/smn.repeat_units.tsv), genes found in the window (results/copies/features.tsv,
gene calls only), SMN copies (results/decomp/copy_units.tsv) and, under the units, the arm elements of step 28b
(results/sv/arm_element_blocks.tsv). (c) The copy classes: each copy is named by its own bases at c.840 (exon 7) and
c.*239 (exon 8). Writes figures/locus_diagram.png/.pdf and results/sv/locus_diagram_examples.tsv."""
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, Polygon, Rectangle
O = 'results/sv'
RU = pd.read_csv('results/decomp/smn.repeat_units.tsv', sep='\t', comment='#')
CU = pd.read_csv('results/decomp/copy_units.tsv', sep='\t')
F = pd.read_csv('results/copies/features.tsv', sep='\t'); F = F[F.kind == 'gene']
W = pd.read_csv('results/copies/windows.tsv', sep='\t'); W = W[W.n_copies > 0].set_index('src')
D = pd.read_csv(f'{O}/palindrome_diagram_rows.tsv', sep='\t')
BLK = pd.read_csv(f'{O}/arm_element_blocks.tsv', sep='\t'); EH = pd.read_csv(f'{O}/arm_element_haps.tsv', sep='\t').set_index('src')
EA = pd.read_csv(f'{O}/arm_elements.tsv', sep='\t'); canon = EA.groupby('src').canonical.all()
SPC = pd.read_csv(f'{O}/arm_element_spacers.tsv', sep='\t')
HU = pd.read_csv('results/decomp/haplotype_units.tsv', sep='\t').set_index('src')
ml = pd.read_csv(f'{O}/marker_lineage.tsv', sep='\t').set_index('copy_id')
D['arr'] = D.src.map(HU.arrangement); D['ncopy'] = D.src.map(lambda s: int(CU[(CU.src == s) & CU.token.str[:-1].isin(['1', '2', 'h12', 'h21', 'd2'])].shape[0]))
D['left'] = -D.lb; D['right'] = D.re
def nearest(sub):
    m = sub[['left', 'right', 'mirror_kb']].median()
    sc = ((sub[['left', 'right', 'mirror_kb']] - m).abs() / sub[['left', 'right', 'mirror_kb']].std()).sum(axis=1)
    return sub.loc[sc.sort_values(kind='stable').index[0]]
D['form'] = D.src.map(EH.cen_type + EH.tel_type); D['canon'] = D.src.map(canon); D['b1b2'] = D.src.map(EH.b1b2)
short = nearest(D[(D.form == 'C1T0') & D.canon & (D.arr == '2- 1+') & (D.ncopy == 2)])
long_ = nearest(D[(D.form == 'C0T1') & D.canon & (D.arr == '1- 1+') & (D.ncopy == 2) & D.b1b2])
import sys; sys.path.insert(0, 'scripts'); from smnlib import FAMCOL, famlab   # same family colours as Figs. 2 and 4
from smnlib import ELCOL, elcol, eltxt, xi_part_labels   # element colors shared with Fig. 4
BLK['lab'] = xi_part_labels(BLK)   # Xi' (inner part) or Xi (outer part), display only
PARTIAL = 0.9   # a gene match covering less of the transcript than this is drawn pale and labelled 'partial'
GCOL = {'SERF1B': '#2c7fb8', 'NAIP': '#8c510a', 'GTF2H2': '#636363', 'GTF2H2C': '#636363', 'OCLN': '#bdbdbd'}
USED = set()   # tick colours drawn in a and b: the legend lists only these
def ccol(c):
    USED.add(_ccol(c)); return _ccol(c)
def _ccol(c):
    k = c.token[:-1]
    if k == '1': return {'A': '#e41a1c', 'B1': '#000000', 'B2': '#a65628', 'B-': '#666666'}[ml.sublineage[c.copy_id]] if c.copy_id in ml.index else '#e41a1c'
    return {'2': '#377eb8', 'h12': '#984ea3', 'h21': '#984ea3'}.get(k, '#999999')
def draw(ax, r, title):
    s = r.src; ax0 = r.axis
    for u in RU[RU.row == s].itertuples():
        x0u, x1u = (u.bgn - ax0) / 1e3, (u.end - ax0) / 1e3; hl = min(28, (x1u - x0u) * 0.35)   # block arrow: points along the unit's strand
        pts = [(x0u, -0.25), (x1u - hl, -0.25), (x1u, 0), (x1u - hl, 0.25), (x0u, 0.25)] if u.strand == '+' else [(x1u, -0.25), (x0u + hl, -0.25), (x0u, 0), (x0u + hl, 0.25), (x1u, 0.25)]
        ax.add_patch(Polygon(pts, closed=True, facecolor=FAMCOL.get(u.family, '#d9d9d9'), alpha=1 if u.strand == '+' else .45, edgecolor='#333333', lw=0.6))
        xm = (u.bgn + u.end) / 2
        if any(abs(c - xm) < 30000 for c in CU[CU.src == s].mid_fwd): xm += 60000 if u.end - xm > 80000 else -60000
        ax.text((xm - ax0) / 1e3, 0, famlab(u.family, u.strand), ha='center', va='center', fontsize=11, color='w' if u.family in ('F', 'J') and u.strand == '+' else 'k', clip_on=True)
    w = W.loc[s, 'window']
    lane_end = [-1e9] * 6   # right edge (kb) of the last label in each lane
    for g in F[(F.src == s) & (F.window == w)].sort_values('win_start').itertuples():
        x0, x1 = (g.win_start - ax0) / 1e3, (g.win_end - ax0) / 1e3
        lw_ = 13 * len(g.gene.replace('SERF1B', 'SERF1').replace('GTF2H2C', 'GTF2H2'))   # label width, ~13 kb per character at 8.5 pt
        cx = (x0 + x1) / 2; lo_, hi_ = min(x0, cx - lw_ / 2), max(x1, cx + lw_ / 2)   # arrow and label together
        L = next(k for k in range(6) if lane_end[k] < lo_ - 15); lane_end[L] = hi_
        y = 0.45 + 0.3 * L
        name = g.gene.replace('SERF1B', 'SERF1').replace('GTF2H2C', 'GTF2H2')
        if x1 - x0 < 8: xm_ = (x0 + x1) / 2; x0, x1 = xm_ - 4, xm_ + 4   # draw short matches at 8 kb so their arrow shows
        dx = (x1 - x0) if g.strand == '+' else -(x1 - x0)
        part = g.tx_frac < PARTIAL
        ax.add_patch(FancyArrow(x0 if g.strand == '+' else x1, y, dx, 0, width=0.06, head_width=0.12, head_length=min(6, abs(dx) / 2), length_includes_head=True, color=GCOL.get(g.gene, '#999'), alpha=.4 if part else 1, lw=0))
        ax.text((x0 + x1) / 2, y + 0.06, name, ha='center', va='bottom', fontsize=8.5, alpha=.6 if part else 1, color=GCOL.get(g.gene, '#555'), clip_on=True)
    for c in CU[CU.src == s].itertuples():
        x = (c.mid_fwd - ax0) / 1e3
        ax.plot([x, x], [-0.45, 0.45], color=ccol(c), lw=3, solid_capstyle='butt', zorder=5)
        lab = {'1': 'SMN1', '2': 'SMN2'}.get(c.token[:-1], c.token[:-1])
        ax.text(x, -0.55, f'{lab} ({"+" if c.token[-1] == "+" else "−"})', ha='center', va='top', fontsize=10, color=ccol(c))
    ax.axvline(0, color='k', lw=.6, ls=':')
    ax.text(0, 2.0, 'axis', ha='center', fontsize=9)
    ax.text(-r.near_left, -1.0, 'cenSMN position', ha='center', fontsize=9)
    ax.text(r.near_right, -1.0, 'telSMN position', ha='center', fontsize=9)
    for q in SPC[(SPC.src == s) & (SPC.site == 'axis')].itertuples():   # S at the main axis
        ax.add_patch(Rectangle((q.d_kb - q.span_kb / 2, -0.40), q.span_kb, 0.11, color=ELCOL['S'], lw=0, zorder=4.5))
    for b in BLK[BLK.src == s].itertuples():   # element strip under the units (arm elements, step 28b)
        x0, x1 = (-b.a1, -b.a0) if b.side == 'cen' else (b.a0, b.a1)
        ax.add_patch(Rectangle((x0, -0.40), x1 - x0, 0.11, color=elcol(b.el, b.lab), lw=0, zorder=4))
        if b.el in ('Xi', 'S', 'Xd', 'E') and x1 - x0 >= 18: ax.text((x0 + x1) / 2, -0.345, b.lab, ha='center', va='center', fontsize=7, color=eltxt(b.lab), zorder=5)
    for q in EA[EA.src == s].itertuples():   # end of the core (step 28b a_end_kb): the R_A unit ends before or beyond it
        x = -q.a_end_kb if q.side == 'cen' else q.a_end_kb
        ax.plot([x, x], [-0.47, 0.30], color='#5a4a2a', lw=1.3, ls='--', zorder=6)
    ax.set_ylim(-1.25, 2.15); ax.set_yticks([])
    for sp in ('left', 'right', 'top'): ax.spines[sp].set_visible(False)
    ax.set_title(title, fontsize=11, loc='left')
fig = plt.figure(figsize=(11, 9.8))
axa = fig.add_axes([0.05, 0.715, 0.92, 0.245]); axb = fig.add_axes([0.05, 0.395, 0.92, 0.245]); axc = fig.add_axes([0.05, 0.01, 0.92, 0.31])
# the whole window of both haplotypes: no unit or gene is cut at the edge
edge = lambda r: (min(RU[RU.row == r.src].bgn.min(), F[(F.src == r.src) & (F.window == W.loc[r.src, 'window'])].win_start.min()) - r.axis, max(RU[RU.row == r.src].end.max(), F[(F.src == r.src) & (F.window == W.loc[r.src, 'window'])].win_end.max()) - r.axis)
xl = (min(edge(short)[0], edge(long_)[0]) / 1e3 - 60, max(edge(short)[1], edge(long_)[1]) / 1e3 + 15)
draw(axa, short, f'a  {short.src.split("_hprc")[0].replace("_", " ")} ("2- 1+"; cenSMN arm C1, telSMN arm T0): one SMN copy on each arm of an inverted R$_A$ pair, X on the cenSMN arm')
draw(axb, long_, f'b  {long_.src.split("_hprc")[0].replace("_", " ")} ("1- 1+", B1 + B2; C0, T1): no X on the cenSMN arm; Xd on the telSMN arm')
for a in (axa, axb): a.set_xlim(*xl); a.set_xlabel('distance from the axis of the R$_A$+ R$_A$− pair (kb)', fontsize=10); a.tick_params(labelsize=10)
# c: copy classes
axc.axis('off'); axc.set_xlim(0, 10); axc.set_ylim(-1.6, 4)
axc.text(0, 3.7, 'c  Copy classes, read from each copy\'s own bases at two sites', fontsize=11)
hdr = ['class', 'c.840 (exon 7)', 'c.*239 (exon 8)', 'exons 7–8']
rowsC = [('SMN1', 'C', 'G', 'present', '#e41a1c'), ('SMN2', 'T', 'A', 'present', '#377eb8'), ('h12 (hybrid)', 'C', 'A', 'present', '#984ea3'),
         ('h21 (hybrid)', 'T', 'G', 'present', '#984ea3'), ('Δ7–8', '–', '–', 'absent', '#999999')]
for j, h in enumerate(hdr): axc.text(0.3 + 1.6 * j, 3.1, h, fontsize=10, weight='bold')
for i, (n, a, b, e, c) in enumerate(rowsC):
    y = 2.6 - 0.5 * i
    for j, v in enumerate([n, a, b, e]): axc.text(0.3 + 1.6 * j, y, v, fontsize=10)
axc.text(6.9, 3.1, 'SMN copy ticks', fontsize=10, weight='bold')
for i, (c, l) in enumerate([x for x in [('#e41a1c', 'SMN1, lineage A'), ('#000000', 'SMN1, lineage B1'), ('#a65628', 'SMN1, lineage B2'), ('#377eb8', 'SMN2'), ('#984ea3', 'hybrid')] if x[0] in USED]):
    axc.plot([6.95, 6.95], [2.55 - .45 * i, 2.8 - .45 * i], color=c, lw=3); axc.text(7.1, 2.6 - .45 * i, l, fontsize=10)
axc.text(6.9, 1.0, 'Unit arrows: repeat-unit family (R$_A$, R$_B$, ...),\npointing along its strand; light = reverse strand. Families\nother than R$_A$ vary between haplotypes (Fig. 2).\nGenes: arrows in their direction of\ntranscription; pale = partial match\n(<90% of the transcript).', fontsize=8.5, va='top')
axc.text(0.3, 0.05, 'Strip under the units: arm elements (core tan, Xi′ pale blue and Xi blue, S magenta, Xd green, E gold, flank gray;\nS also at the axis; Fig. 4). Dashed brown line: the end of the core; an R$_A$ unit ends before it or beyond it.', fontsize=8.5, va='top')
fig.savefig('figures/locus_diagram.png', dpi=220); fig.savefig('figures/locus_diagram.pdf', metadata={'CreationDate': None})
pd.DataFrame([dict(panel='a', **short.to_dict()), dict(panel='b', **long_.to_dict())]).to_csv(f'{O}/locus_diagram_examples.tsv', sep='\t', index=False)
print(short.src, long_.src)
