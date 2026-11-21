# Pipeline step 17 (distance_context); ported from explore/structure_vs_sequence/13_distance_context.py (sha256 a893c19ae84c, commit 5e2b3bb),
# paths rewritten to results/sv/ and figures/; run from the project root after the steps before it.
"""Exploration step 13: within SMN1 lineage A and within lineage B, does a copy's mean distance to the SMN2
copies (out/group_copy_means.tsv, gene body, per kb) go with where the copy sits: its position (centromeric /
telomeric / single), whether it is on an arm of the family-A palindrome (out/palindrome_axis_copies.tsv), its
haplotype's arrangement type, its downstream environment and its 3' NAIP form (out/copy_labels.tsv)?
Tests: Mann-Whitney (two levels) or Kruskal-Wallis (more), within each lineage, all copies and AFR copies only.
Copies are not independent (two per haplotype, related haplotypes), so p-values are descriptive.
Writes out/distance_context.tsv and out/distance_context_tests.tsv."""
import pandas as pd
from scipy.stats import mannwhitneyu, kruskal
O = 'results/sv'
g = pd.read_csv(f'{O}/group_copy_means.tsv', sep='\t')
l = pd.read_csv(f'{O}/copy_labels.tsv', sep='\t')
pal = pd.read_csv(f'{O}/palindrome_axis_copies.tsv', sep='\t')
l['token'] = l.copy_id.map(lambda c: None)
tips = pd.read_csv(f'{O}/copy_tips.tsv', sep='\t')[['tip', 'token']]
d = g[g.group.isin(['A', 'B'])].merge(l, on='tip').drop(columns='token').merge(tips, on='tip')
d = d[~d.src.str.lower().str.startswith(('gca_000001405', 'gca_009914755', 'chm13'))]   # GRCh38 and CHM13 (source name 'chm13v2.0...')
keyp = set(zip(pal.src, pal.token))
d['palindrome_arm'] = [('yes' if (s, t) in keyp else 'no') for s, t in zip(d.src, d.token)]
vc = d.arrangement.value_counts()
d['arrangement_top'] = d.arrangement.where(d.arrangement.map(vc) >= 10, 'other')
# arm of the head-to-head R_A pair, as in Fig. 1 and Fig. 3: cenSMN (left of the axis) or telSMN (right), only on haplotypes
# with a copy on both arms (the palindrome set); every other copy 'no R_A pair'
both = {s for s, x in pal.groupby('src') if (x.d_kb < 0).any() and (x.d_kb > 0).any()}
side = {(s, t): ('cenSMN' if k < 0 else 'telSMN') for s, t, k in zip(pal.src, pal.token, pal.d_kb) if s in both}
d['arm'] = [side.get((s, t), 'no pair') for s, t in zip(d.src, d.token)]
d['position3'] = d.position.where(d.position.isin(['centromeric', 'telomeric', 'single']), 'other')
ctx = ['arm', 'position3', 'palindrome_arm', 'arrangement_top', 'env', 'naip3']
rows, tests = [], []
for subset, dd in [('all', d), ('AFR', d[d.afr.astype(bool)])]:
    for lin in ['A', 'B']:
        x = dd[dd.group == lin]
        for c in ctx:
            lv = [v for v, n in x[c].value_counts().items() if n >= 5 and v != 'other' and pd.notna(v)]
            for v in lv:
                y = x[x[c] == v]
                rows.append(dict(subset=subset, lineage=lin, context=c, level=v, copies=len(y),
                                 median_to_SMN2=round(y.mean_to_SMN2.median(), 3), median_to_A=round(y.mean_to_A.median(), 3),
                                 median_to_B=round(y.mean_to_B.median(), 3)))
            if len(lv) < 2: continue
            gr = [x[x[c] == v].mean_to_SMN2 for v in lv]
            p = mannwhitneyu(*gr).pvalue if len(gr) == 2 else kruskal(*gr).pvalue
            tests.append(dict(subset=subset, lineage=lin, context=c, levels=len(lv), copies=sum(map(len, gr)),
                              test='Mann-Whitney' if len(gr) == 2 else 'Kruskal-Wallis', p=float('%.2g' % p),
                              spread=round(max(v.median() for v in gr) - min(v.median() for v in gr), 3)))
R = pd.DataFrame(rows); T = pd.DataFrame(tests)
R.to_csv(f'{O}/distance_context.tsv', sep='\t', index=False); T.to_csv(f'{O}/distance_context_tests.tsv', sep='\t', index=False)
print(R.to_string()); print(T.to_string())
# A copies: by the SMN1 lineages on their haplotype (A only / two A / one A + one B) and by Chen 2023 haplogroup
h = pd.read_csv(f'{O}/lineage_haplotypes.tsv', sep='\t')[['src', 'composition']]
c = pd.read_csv(f'{O}/chen_hg_copies.tsv', sep='\t')[['tip', 'hg']]
A = d[d.group == 'A'].merge(h, on='src').merge(c, on='tip', how='left')
A['hg1'] = A.hg.where(~A.hg.fillna('').str.contains(r'\|') | A.hg.isin(['S1-2|S1-2c']), 'tied')
rows = []
for col in ['composition', 'hg1']:
    for v, y in A.groupby(col):
        rows.append(dict(by=col, level=v, copies=len(y), afr=int(y.afr.astype(bool).sum()),
                         median_to_SMN2=round(y.mean_to_SMN2.median(), 3), median_to_B=round(y.mean_to_B.median(), 3)))
HA = pd.DataFrame(rows); HA.to_csv(f'{O}/distance_A_by_haplotype.tsv', sep='\t', index=False)
X = pd.crosstab(A.hg1, A.composition); X.to_csv(f'{O}/distance_A_hg_by_composition.tsv', sep='\t')
print(HA.to_string()); print(X.to_string())
