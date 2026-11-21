# Pipeline step 06 (lineage_detail); ported from explore/structure_vs_sequence/06_lineage_detail.py (sha256 fb5c5f9866e6, commit 5e2b3bb),
# paths rewritten to results/sv/ and figures/; run from the project root after the steps before it.
"""Exploration step 6: the two SMN1 lineages of the copy tree in detail (A = other SMN1, B = the SMN1-only clade that
holds most c.*3+80T>G copies, doc/marker_lineage.tsv in_marker_subclade). Cross-tabulations by position x strand,
arrangement type, 1000 Genomes population, downstream context, and per-haplotype lineage composition; one figure.
Reads out/copy_labels.tsv, doc/marker_lineage.tsv, tables/haplotype_copy_counts.tsv. Nothing here feeds the note."""
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.stats import fisher_exact
O = 'results/sv'; F = 'figures'
lab = pd.read_csv(f'{O}/copy_labels.tsv', sep='\t'); mk = pd.read_csv(f'{O}/marker_lineage.tsv', sep='\t')
hc = pd.read_csv('tables/haplotype_copy_counts.tsv', sep='\t')[['sample', 'hap', 'src', 'population', 'resolved']]
d = lab.merge(mk[['tip', 'c3p80_alt', 'in_marker_subclade']], on='tip', how='left')
assert d.in_marker_subclade.notna().all()
d = d[d.cls == 'SMN1'].copy(); d['lineage'] = np.where(d.in_marker_subclade, 'B', 'A')
d = d.merge(hc, on='src', how='left')
d.loc[d.superpop.isna(), 'population'] = 'reference'
d['position'] = d.position.fillna('unresolved'); d['superpop'] = d.superpop.fillna('reference')
d['group'] = np.where(d.superpop == 'AFR', 'AFR', 'non-AFR')
d.to_csv(f'{O}/lineage_copies.tsv', sep='\t', index=False)

rows = []
for g, sub in [('all', d), ('AFR', d[d.group == 'AFR']), ('non-AFR', d[d.group == 'non-AFR'])]:
    for (lin, pos, st), n in sub.groupby(['lineage', 'position', 'strand']).size().items():
        rows.append(dict(group=g, lineage=lin, position=pos, strand=st, copies=n))
pd.DataFrame(rows).to_csv(f'{O}/lineage_position_strand.tsv', sep='\t', index=False)
ar = pd.crosstab(d.arrangement.fillna('unresolved'), d.lineage); ar['B_share'] = ar.B / (ar.A + ar.B)
ar = ar.sort_values(['B', 'A'], ascending=False); ar.to_csv(f'{O}/lineage_arrangement.tsv', sep='\t')
po = d.groupby(['superpop', 'population']).agg(A=('lineage', lambda s: (s == 'A').sum()), B=('lineage', lambda s: (s == 'B').sum()),
                                                haplotypes=('src', 'nunique')).reset_index()
po['B_share'] = po.B / (po.A + po.B); po.to_csv(f'{O}/lineage_population.tsv', sep='\t', index=False)
ctx = []
for c_ in ['env', 'naip3', 'c3p80_alt']:
    for (lin, v), n in d.groupby(['lineage', d[c_].astype(str)]).size().items():
        ctx.append(dict(label=c_, value=v, lineage=lin, copies=n))
pd.DataFrame(ctx).to_csv(f'{O}/lineage_context.tsv', sep='\t', index=False)
h = d.groupby('src').agg(sample=('sample', 'first'), hap=('hap', 'first'), population=('population', 'first'),
                         superpop=('superpop', 'first'), resolved=('resolved', 'first'), arrangement=('arrangement', 'first'),
                         nA=('lineage', lambda s: (s == 'A').sum()), nB=('lineage', lambda s: (s == 'B').sum()),
                         marker_copies=('c3p80_alt', 'sum')).reset_index()
h['composition'] = h.apply(lambda r: 'A' * r.nA + 'B' * r.nB, axis=1)
h['B_copies'] = h.src.map(lambda s_: ';'.join(f'{p}{st}' for p, st in zip(*d[(d.src == s_) & (d.lineage == 'B')][['position', 'strand']].values.T)) if ((d.src == s_) & (d.lineage == 'B')).any() else '')
h.to_csv(f'{O}/lineage_haplotypes.tsv', sep='\t', index=False)
comp = pd.crosstab(h.composition, h.superpop); comp.to_csv(f'{O}/lineage_composition.tsv', sep='\t')
ha = h[(h.superpop == 'AFR') & (h.resolved == True)]
t = [[int(((ha.nB > 0) & (ha.nA + ha.nB >= 2)).sum()), int(((ha.nB > 0) & (ha.nA + ha.nB == 1)).sum())],
     [int(((ha.nB == 0) & (ha.nA + ha.nB >= 2)).sum()), int(((ha.nB == 0) & (ha.nA + ha.nB == 1)).sum())]]
pd.DataFrame([dict(test='AFR resolved haplotypes with SMN1: carries a B copy vs carries >=2 SMN1 copies', B_multi=t[0][0], B_single=t[0][1],
                   noB_multi=t[1][0], noB_single=t[1][1], p=fisher_exact(t)[1])]).to_csv(f'{O}/lineage_tests.tsv', sep='\t', index=False)

# figure: lineage colours of the shared scheme (scripts/smnlib.py COPYCOL); B split into B1, B2 and B- (step 05 sublineage)
import sys; sys.path.insert(0, 'scripts'); from smnlib import COPYCOL
d['sub'] = d.tip.map(mk.set_index('tip').sublineage).where(d.lineage == 'B', 'A')
assert set(d['sub']) <= {'A', 'B1', 'B2', 'B-'}, set(d['sub'])
SUBS = ['A', 'B1', 'B2', 'B-']; col = {k: COPYCOL[k] for k in SUBS}; col['A'] = '#d98880'  # softer red for this dense figure (researcher: too red)
SUBLAB = {'A': 'A', 'B1': 'B1', 'B2': 'B2', 'B-': 'B, neither clade'}
fig, ax = plt.subplots(2, 2, figsize=(12, 9.5))
a = ax[0, 0]; cats = [(p, s) for p in ['centromeric', 'middle', 'telomeric', 'single', 'unresolved'] for s in ['+', '-']]
x = np.arange(len(cats)); w = 0.21
for i, lin in enumerate(SUBS):
    v = [int(((d['sub'] == lin) & (d.position == p) & (d.strand == s)).sum()) for p, s in cats]
    a.bar(x + (i - 1.5) * w, v, w, color=col[lin], label=f'{SUBLAB[lin]} ({(d["sub"] == lin).sum()} copies)')
    for xi, vi in zip(x, v):
        if vi: a.text(xi + (i - 1.5) * w, vi + 3, str(vi), ha='center', fontsize=6)
a.set_xticks(x); a.set_xticklabels([f'{p} {s}' for p, s in cats], rotation=45, ha='right', fontsize=8)
a.set_ylabel('SMN1 copies'); a.set_title('a  Position and strand of each SMN1 copy', loc='left', fontsize=10); a.legend(fontsize=8, frameon=False)
hs = d.groupby('src').agg(superpop=('superpop', 'first'), subs=('sub', lambda s: '+'.join(sorted(s))))
def hcat(c):
    p = c.split('+')
    if set(p) == {'A'}: return 'A' if len(p) == 1 else 'A, two or more'
    if 'B1' in p and 'B2' in p: return 'B1 + B2'
    if 'A' in p: return 'A + B'
    return 'B only, not B1 + B2'
hs['cat'] = hs.subs.map(hcat)
order = ['A', 'A, two or more', 'A + B', 'B1 + B2', 'B only, not B1 + B2']; sps = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']
ccol = {'A': '#f2d0cc', 'A, two or more': '#c0605a', 'A + B': '#7f6a6a', 'B1 + B2': COPYCOL['B1'], 'B only, not B1 + B2': COPYCOL['B2']}
bot = np.zeros(len(sps)); a = ax[0, 1]
for c_ in order:
    v = np.array([int(((hs.superpop == s) & (hs['cat'] == c_)).sum()) for s in sps]); a.bar(sps, v, bottom=bot, color=ccol[c_], label=c_)
    for xi, (vi, b0) in enumerate(zip(v, bot)):
        if vi >= 6: a.text(xi, b0 + vi / 2, str(vi), ha='center', va='center', fontsize=7, color='w' if c_ in ('A, two or more', 'A + B', 'B1 + B2', 'B only, not B1 + B2') else 'k')
    bot += v
a.set_ylabel('haplotypes with an SMN1 copy in the tree'); a.legend(title='SMN1 copies by lineage', fontsize=8, title_fontsize=8, frameon=False, ncol=3, loc='upper right')
a.set_ylim(0, bot.max() * 1.2); a.set_title('b  Lineages of the SMN1 copies on each haplotype', loc='left', fontsize=10)
def stack(a, df, idx):
    left = np.zeros(len(df))
    for lin in SUBS:
        v = df[lin].values; a.barh(idx, v, left=left, color=col[lin], label=SUBLAB[lin]); left += v
pop_sub = pd.crosstab(d.population, d['sub']).reindex(columns=SUBS, fill_value=0)
a = ax[1, 0]; pp = po[po.superpop == 'AFR'].sort_values('B_share', ascending=False)
stack(a, pop_sub.loc[pp.population], pp.population)
for yi, r in enumerate(pp.itertuples()): a.text(r.A + r.B + 1, yi, f'B {r.B} of {r.A + r.B}', va='center', fontsize=7)
a.invert_yaxis(); a.set_xlabel('SMN1 copies'); a.set_title('c  AFR populations: SMN1 copies by lineage', loc='left', fontsize=10)
a.set_xlim(0, (pp.A + pp.B).max() * 1.3); a.legend(fontsize=8, frameon=False, loc='lower right')
a = ax[1, 1]; top = ar[ar.index != 'unresolved'].copy(); top['n'] = top.A + top.B; top = top[(top.B > 0) | (top.n >= 10)].sort_values('n')
ars = pd.crosstab(d.arrangement.fillna('unresolved'), d['sub']).reindex(columns=SUBS, fill_value=0).loc[top.index]
stack(a, ars, top.index)
for yi, r in enumerate(top.itertuples()): a.text(r.A + r.B + 2, yi, f'B {r.B} of {r.A + r.B}', va='center', fontsize=7)
a.set_xlabel('SMN1 copies'); a.set_title('d  Arrangement types: SMN1 copies by lineage', loc='left', fontsize=10)
a.set_xlim(0, top.n.max() * 1.2); a.tick_params(axis='y', labelsize=8); a.legend(fontsize=8, frameon=False, loc='lower right')
fig.tight_layout(); fig.savefig(f'{F}/lineage_detail.png', dpi=200); fig.savefig(f'{F}/lineage_detail.pdf', metadata={'CreationDate': None})
print(ar.head(12)); print(comp); print(po.to_string()); print(pd.read_csv(f'{O}/lineage_tests.tsv', sep='\t').T)
print(h[h.composition.isin(['B', 'BB'])].sort_values(['superpop', 'population'])[['sample', 'hap', 'population', 'resolved', 'arrangement', 'composition', 'B_copies', 'marker_copies']].to_string())
print(pd.read_csv(f'{O}/lineage_context.tsv', sep='\t').to_string())
