"""Pipeline step 39 (supplement figure): how a window is taken and how PANGEA repeat-units turns its bundle
decomposition into repeat units and families, drawn on real data.
Inputs: results/sv/window_flanks.tsv and .json (step 36), results/decomp/smn.bed and smn.repeat_units.tsv (step 00).
Outputs: figures/repeat_units_algo.png/.pdf, results/sv/repeat_units_algo.json (the numbers the figure prints)."""
import json, collections
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrow
import sys; sys.path.insert(0, 'scripts')
from smnlib import FAMCOL, FAM_OTHER
plt.rcParams['svg.hashsalt'] = 'smn'
EX = 'HG01099_mat_hprc_r2_v1.0.1'
WF = pd.read_csv('results/sv/window_flanks.tsv', sep='\t'); WJ = json.load(open('results/sv/window_flanks.json'))
seg = collections.defaultdict(list)
for l in open('results/decomp/smn.bed'):
    if l.startswith('#'): continue
    s, b, e, n = l.rstrip('\n').split('\t')[:4]
    f = n.split(':'); seg[s].append((int(b), int(e), int(f[0]), '+' if f[2] == '0' else '-'))
for v in seg.values(): v.sort()
RU = pd.read_csv('results/decomp/smn.repeat_units.tsv', sep='\t', comment='#')
ex = seg[EX]; exu = RU[RU.row == EX].sort_values('bgn')
assert len(exu) == 6 and list(exu.family) == ['C', 'B', 'A', 'A', 'D', 'B'], 'example units changed'
c3, c4 = exu.iloc[2], exu.iloc[3]
assert c3.strand == '+' and c4.strand == '-', 'example: R_A+ then R_A- at the axis'
axis = (c3.end + c4.bgn) / 2
wr = WF.set_index('src').loc[EX]
short = WF.loc[WF.length.idxmin()]; long_ = WF.loc[WF.length.idxmax()]
ids = sorted({x[2] for v in seg.values() for x in v})
cm = plt.get_cmap('tab20'); cm2 = plt.get_cmap('tab20b')
BC = {b: (cm(i) if i < 20 else cm2(i - 20)) for i, b in enumerate(ids)}

fig = plt.figure(figsize=(13, 12))
gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 2.4, 1.05], width_ratios=[1.35, 1], hspace=0.42, wspace=0.18)
# a: how the window is taken
ax = fig.add_subplot(gs[0, :]); kb = 1e3
Q = WJ['query_len'] / kb
ax.add_patch(Rectangle((0, 3), Q, 0.5, color='#555'))
ax.text(Q / 2, 3.75, f"GRCh38 query chr5:{WJ['grch38_start1']:,}–{WJ['grch38_end']:,} ({Q:,.0f} kb)", ha='center', fontsize=9)
rowsA = [(EX.split('_hprc')[0].replace('_', ' '), wr.length, wr.cen_unique_bp, wr.tel_unique_bp, 'example below'),
         (short.src.split('_hprc')[0].replace('_', ' '), short.length, short.cen_unique_bp, short.tel_unique_bp, 'shortest'),
         (long_.src.split('_hprc')[0].replace('_', ' '), long_.length, long_.cen_unique_bp, long_.tel_unique_bp, 'longest')]
for k, (nm, L, cu, tu, why) in enumerate(rowsA):
    y = 2 - k * 0.9; L /= kb; cu /= kb; tu /= kb
    ax.add_patch(Rectangle((0, y), cu, 0.45, color='#bdbdbd')); ax.add_patch(Rectangle((cu, y), L - cu - tu, 0.45, color='#e6550d', alpha=.55))
    ax.add_patch(Rectangle((L - tu, y), tu, 0.45, color='#bdbdbd'))
    ax.text(L + 25, y + 0.22, f'{nm}: {L:,.0f} kb window ({why})', va='center', fontsize=8.5)
ax.text(-15, 1.1, 'each haplotype:\nthe contig stretch\nthe query hits', ha='right', va='center', fontsize=8.5)
ax.set_xlim(-330, max(Q, long_.length / kb) + 520); ax.set_ylim(-0.1, 4.2); ax.axis('off')
ax.add_patch(Rectangle((max(Q, long_.length / kb) + 40, 3.0), 30, 0.3, color='#bdbdbd')); ax.text(max(Q, long_.length / kb) + 80, 3.15, 'single-copy bundles (unique flank)', fontsize=8, va='center')
ax.add_patch(Rectangle((max(Q, long_.length / kb) + 40, 2.5), 30, 0.3, color='#e6550d', alpha=.55)); ax.text(max(Q, long_.length / kb) + 80, 2.65, 'first to last bundle that recurs in the window', fontsize=8, va='center')
ax.set_title(f"a  The window: all {WJ['n_windows']} start within 2 kb of the query start and begin and end in single-copy bundles "
             f"(≥{WJ['cen_unique_kb_min']:.1f} kb centromeric, ≥{WJ['tel_unique_kb_min']:.1f} kb telomeric)", loc='left', fontsize=10)
# b: self-match dot plot in bundles
ax = fig.add_subplot(gs[1, 0]); L = wr.length / kb
nd = ni = 0
for i, (b1, e1, i1, s1) in enumerate(ex):
    for j, (b2, e2, i2, s2) in enumerate(ex):
        if i1 != i2: continue   # i == j: a segment with itself, the main diagonal of a self dot plot
        x, y = (b1 + e1) / 2 / kb, (b2 + e2) / 2 / kb
        d = s1 == s2; ax.plot(x, y, 'o', ms=2.6, color='#2b6cb0' if d else '#d6604d', zorder=3)
        if i < j: nd += d; ni += (not d)
ax.plot([0, L], [0, L], color='#ccc', lw=.6)
for _, u in exu.iterrows():
    ax.axvspan(u.bgn / kb, u.end / kb, color=FAMCOL.get(u.family, FAM_OTHER), alpha=.13, lw=0)
ax.axvline(axis / kb, ls=':', color='#555', lw=.7, zorder=1); ax.axhline(axis / kb, ls=':', color='#555', lw=.7, zorder=1)
ax.text(axis / kb + 12, L * 0.03, 'palindrome\ncenter', fontsize=8)
ax.set_xlim(0, L); ax.set_ylim(0, L); ax.set_aspect('equal')
ax.set_xlabel(f'{rowsA[0][0]} window (kb)'); ax.set_ylabel('same window (kb)')
ax.plot([], [], 'o', color='#2b6cb0', label='direct match (same bundle, same strand;\nincl. each segment with itself, the diagonal)')
ax.plot([], [], 'o', color='#d6604d', label='inverted match (same bundle, opposite strand)')
ax.legend(loc='upper left', fontsize=7.5, framealpha=.9)
ax.set_title('b  Steps 1–3: every two segments of one bundle match (a dot);\nmatches chain along diagonals; chains are cut at palindrome centers', loc='left', fontsize=10)
# c: the steps in words, next to b
ax = fig.add_subplot(gs[1, 1]); ax.axis('off')
steps = [('1  Self-match', 'two segments with the same bundle id match:\ndirect on the same strand, inverted on opposite strands'),
         ('2  Chains', 'matches join along a diagonal; a step may skip ≤3 segments\non both copies or ≤12 on one (an indel); score = matched bp\nminus ¼ of the difference in bp skipped on the two copies'),
         ('3  Palindrome centers', 'where segments mirror outward from a point for at least\nthe minimum unit span, chains crossing that point are cut'),
         ('4  Tandem arrays', 'a chain whose copies overlap at one offset is an array\nof ≥2.5 units; its copies are snapped to whole units'),
         ('5  Copies', 'chains taken best first; each copy keeps only segments no\nearlier copy holds; copies found by one chain are linked\ninto one repeat (here one repeat with 6 copies)'),
         ('6  Families', 'over all windows, copies longest first: each joins the family\nwhose mean bundle content (bp per bundle id) is most similar,\nif bp-weighted Jaccard ≥0.5; otherwise it founds a family')]
y = 0.98
for h, t in steps:
    ax.text(0, y, h, fontsize=9.5, weight='bold', va='top'); ax.text(0.04, y - 0.045, t, fontsize=8.3, va='top')
    y -= 0.045 + 0.04 * (t.count('\n') + 1) + 0.03
ax.set_title('The steps of repeat-units', loc='left', fontsize=10)
# d: segments and the units it calls
ax = fig.add_subplot(gs[2, :])
for b, e, i, s in ex:
    ax.add_patch(Rectangle((b / kb, 1.0), (e - b) / kb, 0.5, color=BC[i], alpha=1 if s == '+' else .45, lw=0))
for _, u in exu.iterrows():
    x0, x1 = (u.bgn, u.end) if u.strand == '+' else (u.end, u.bgn)
    ax.add_patch(FancyArrow(x0 / kb, 0.3, (x1 - x0) / kb, 0, width=0.36, head_width=0.5, head_length=18, length_includes_head=True,
                            color=FAMCOL.get(u.family, FAM_OTHER)))
    ax.text((u.bgn + u.end) / 2 / kb, -0.25, f'R$_{{{u.family}}}${u.strand}', ha='center', fontsize=9)
ax.axvline(axis / kb, ls='--', color='k', lw=1)
ax.text(-10, 1.25, 'bundle\nsegments', ha='right', va='center', fontsize=8.5); ax.text(-10, 0.3, 'repeat\nunits', ha='right', va='center', fontsize=8.5)
ax.set_xlim(-120, L + 10); ax.set_ylim(-0.6, 1.7); ax.set_yticks([]); ax.set_xlabel('kb along the window')
for sp in ('left', 'right', 'top'): ax.spines[sp].set_visible(False)
ax.set_title(f'c  Input and output for {rowsA[0][0]}: {len(ex)} bundle segments (one color per bundle id; pale, reverse strand) '
             f'and the {len(exu)} copies repeat-units calls, by family', loc='left', fontsize=10)
fig.savefig('figures/repeat_units_algo.png', dpi=150, bbox_inches='tight')
fig.savefig('figures/repeat_units_algo.pdf', bbox_inches='tight', metadata={'CreationDate': None, 'ModDate': None})
json.dump(dict(example=EX, n_segments=len(ex), n_units=len(exu), families=list(exu.family), axis_kb=round(axis / kb, 1),
               n_direct_pairs=nd, n_inverted_pairs=ni, short_src=short.src, short_kb=round(short.length / kb, 1),
               long_src=long_.src, long_kb=round(long_.length / kb, 1)), open('results/sv/repeat_units_algo.json', 'w'), indent=1)
print(open('results/sv/repeat_units_algo.json').read())
