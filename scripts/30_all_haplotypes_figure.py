"""Pipeline step 30: every resolved haplotype in one figure. Each row is one haplotype's copy-bearing window (the
window pangea decomposed, results/decomp/smn.bed), drawn to scale from its centromeric end: its large repeat units
(pangea cohort bundle repeat-units, results/decomp/smn.repeat_units.tsv) as bars coloured by family, light for the
reverse strand, and its SMN copies (results/decomp/copy_units.tsv) as ticks coloured by class and SMN1 lineage
(results/sv/marker_lineage.tsv). Each row is aligned on its most telomeric counted SMN copy (x = 0). Rows are grouped by arrangement type (most common first; types on fewer than 5
haplotypes pooled as 'other types') and ordered within a group by UPGMA on the structural distance (step 01).
A strip gives each haplotype's continental group. Writes figures/all_haplotypes.png/.pdf and
results/sv/all_haplotypes_rows.tsv (row order and what each row holds)."""
import numpy as np, pandas as pd, collections, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import average, optimal_leaf_ordering, leaves_list
from scipy.spatial.distance import squareform
O = 'results/sv'
RU = pd.read_csv('results/decomp/smn.repeat_units.tsv', sep='\t', comment='#')
CU = pd.read_csv('results/decomp/copy_units.tsv', sep='\t')
hc = pd.read_csv('tables/haplotype_copy_counts.tsv', sep='\t')
W = pd.read_csv('results/copies/windows.tsv', sep='\t'); W = W[W.n_copies > 0].set_index('src')
ml = pd.read_csv(f'{O}/marker_lineage.tsv', sep='\t').set_index('copy_id')
NOT_IN_TREE = 0   # SMN1 copies drawn that are not in the copy tree (legend entry only if any)
AFR = {'ACB', 'ASW', 'ESN', 'GWD', 'LWK', 'MSL', 'YRI', 'MKK'}
SUPER = {**{p: 'AFR' for p in AFR}, **{p: 'AMR' for p in ['PUR', 'PEL', 'CLM', 'MXL']},
         **{p: 'EAS' for p in ['JPT', 'KHV', 'CHS', 'CDX', 'CHB']}, **{p: 'EUR' for p in ['GBR', 'FIN', 'TSI', 'IBS']},
         **{p: 'SAS' for p in ['BEB', 'PJL', 'GIH', 'STU', 'ITU']}}
SUPER_COL = {'AFR': '#8c510a', 'AMR': '#d73027', 'EAS': '#4575b4', 'EUR': '#1a9850', 'SAS': '#984ea3', 'other': '#bdbdbd'}
res = hc[hc.resolved & hc.src.isin(set(RU.row))].copy()
res['superpop'] = res.population.map(lambda p: SUPER.get(p, 'other'))
tc = res.arrangement.value_counts()
res['group'] = res.arrangement.where(res.arrangement.map(tc) >= 5, 'other types')
order_groups = [g for g in tc.index if tc[g] >= 5] + ['other types']
# structural distance within each group
st = pd.read_csv(f'{O}/smn_struct.dist', sep=' ', header=None, usecols=[0, 1, 2], names=['a', 'b', 'd'])
st = st[st.a.isin(set(res.src)) & st.b.isin(set(res.src))]
SD = {(a, b): d for a, b, d in st.itertuples(index=False)}
rows = []
for g in order_groups:
    s = sorted(res.src[res.group == g])
    if len(s) > 2:
        M = np.array([[0 if a == b else (SD.get((a, b), SD.get((b, a))) ) for b in s] for a in s], float)
        M = (M + M.T) / 2; np.fill_diagonal(M, 0); y = squareform(M, checks=False)
        s = [s[i] for i in leaves_list(optimal_leaf_ordering(average(y), y))]
    rows += [(g, x) for x in s]
assert len(rows) == len(res)
fam_tot = RU[RU.row.isin(set(res.src))].family.value_counts()
FAMS = list(fam_tot.index[:8])
import sys; sys.path.insert(0, 'scripts'); from smnlib import FAMCOL, famlab
assert set(FAMS) == set(FAMCOL), (FAMS, FAMCOL)   # the eight most frequent families are the eight with a colour
FCOL = {f: FAMCOL[f] for f in FAMS}
def copycol(r):
    k = r.token[:-1]
    if k == '1':
        if r.copy_id not in ml.index:
            global NOT_IN_TREE; NOT_IN_TREE += 1; return '#fb9a99'   # SMN1 copy not in the copy tree
        return {'A': '#e41a1c', 'B1': '#000000', 'B2': '#a65628', 'B-': '#666666'}[ml.sublineage[r.copy_id]]
    return {'2': '#377eb8', 'h12': '#984ea3', 'h21': '#984ea3'}.get(k, '#999999')
XL, XR = [], []
x0s = {}
for _, s in rows:
    cs_ = CU[CU.src == s]; cn_ = cs_[cs_.token.str[:-1].isin(['1', '2', 'h12', 'h21', 'd2'])]
    x0s[s] = cn_.mid_fwd.max() if len(cn_) else cs_.mid_fwd.max()
XMIN = min(-x0s[s] / 1e3 for _, s in rows); XMAX = max((W.length[s] - x0s[s]) / 1e3 for _, s in rows)
n = len(rows); fig = plt.figure(figsize=(11, 14)); ax = fig.add_axes([0.2, 0.08, 0.77, 0.88])
out = []; ygroup = collections.defaultdict(list)
for i, (g, s) in enumerate(rows):
    y = n - 1 - i; ygroup[g].append(y)
    cs = CU[CU.src == s]
    cnt = cs[cs.token.str[:-1].isin(['1', '2', 'h12', 'h21', 'd2'])]
    x0 = cnt.mid_fwd.max() if len(cnt) else cs.mid_fwd.max()
    for u in RU[RU.row == s].itertuples():
        c = FCOL.get(u.family, '#d9d9d9')
        ax.add_patch(plt.Rectangle(((u.bgn - x0) / 1e3, y - .42), u.span_bp / 1e3, .84, color=c, alpha=1 if u.strand == '+' else .45, lw=0))
    for c in cs.itertuples():
        ax.plot([(c.mid_fwd - x0) / 1e3] * 2, [y - .7, y + .7], color=copycol(c), lw=2.2, solid_capstyle='butt', zorder=5)
    XL.append((-x0) / 1e3); XR.append((W.length[s] - x0) / 1e3)
    ax.add_patch(plt.Rectangle((XMIN - 80, y - .5), 60, 1, color=SUPER_COL[res.set_index('src').superpop[s]], lw=0))
    out.append(dict(order=i, src=s, group=g, window_kb=round(W.length[s] / 1e3, 1), aligned_at_kb=round(x0 / 1e3, 1),
                    units=' '.join(f'{u.family}{u.strand}' for u in RU[RU.row == s].sort_values('bgn').itertuples()),
                    copies=' '.join(c.token for c in cs.sort_values('mid_fwd').itertuples())))
for g, ys in ygroup.items():
    ax.plot([XMIN - 10, XMIN - 10], [min(ys) - .4, max(ys) + .4], color='k', lw=1)
    ax.text(XMIN - 100, (min(ys) + max(ys)) / 2, f'{g}  ({len(ys)})', ha='right', va='center', fontsize=7)
    ax.axhline(min(ys) - .5, color='#bbbbbb', lw=.3)
ax.set_xlim(XMIN - 100, XMAX + 20); ax.axvline(0, color='k', lw=.4, ls=':', zorder=0); ax.set_ylim(-1, n)
ax.set_yticks([]); ax.set_xlabel('position (kb) relative to the most telomeric SMN copy; centromeric to the left')
for sp in ('left', 'right', 'top'): ax.spines[sp].set_visible(False)
h1 = [plt.Rectangle((0, 0), 1, 1, color=FCOL[f], label=f'unit family {famlab(f)}') for f in FAMS] + [plt.Rectangle((0, 0), 1, 1, color='#d9d9d9', label='other families')]
h2 = [plt.Line2D([], [], color=c, lw=2, label=l) for c, l in [('#e41a1c', 'SMN1, lineage A'), ('#000000', 'SMN1, lineage B1 (S1-9/S1-6)'), ('#a65628', 'SMN1, lineage B2 (S1-8)'), ('#666666', 'SMN1, lineage B, in neither'), ('#377eb8', 'SMN2'), ('#984ea3', 'hybrid (h12, h21)')] + ([('#fb9a99', 'SMN1, not in the copy tree')] if NOT_IN_TREE else []) + [('#999999', 'other copy or piece')]]
h3 = [plt.Rectangle((0, 0), 1, 1, color=c, label=('ASL or no label' if k == 'other' else k)) for k, c in SUPER_COL.items()]
l1 = fig.legend(handles=h1, loc='lower left', bbox_to_anchor=(0.03, 0.0), ncol=3, fontsize=7, frameon=False, title='repeat-unit family (bars)', title_fontsize=7.5)
l2 = fig.legend(handles=h2, loc='lower center', bbox_to_anchor=(0.575, 0.0), ncol=3, columnspacing=0.8, handlelength=1.5, fontsize=7, frameon=False, title='SMN copy (ticks)', title_fontsize=7.5)
l3 = fig.legend(handles=h3, loc='lower right', bbox_to_anchor=(0.99, 0.0), ncol=2, fontsize=7, frameon=False, title='continental group (left strip)', title_fontsize=7.5)
fig.text(0.2, 0.965, 'unit bars: dark = forward strand, light = reverse strand; strip at left: continental group', fontsize=7)
fig.savefig('figures/all_haplotypes.png', dpi=220); fig.savefig('figures/all_haplotypes.pdf', metadata={'CreationDate': None})
pd.DataFrame(out).to_csv(f'{O}/all_haplotypes_rows.tsv', sep='\t', index=False)
shown = CU[CU.src.isin([s for _, s in rows])]
print('tick colors:', collections.Counter(copycol(c) for c in shown.itertuples()))
print(n, 'haplotypes;', len(order_groups), 'groups;', 'families colored', FAMS)
