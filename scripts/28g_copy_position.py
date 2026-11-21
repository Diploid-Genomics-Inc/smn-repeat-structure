"""Pipeline step 28g: why the SMN copy sits further from the axis on C0 (cenSMN) and T1 (telSMN) arms.

For each arm of the 270 A+ A- pair haplotypes (step 28b: results/sv/arm_element_segments.tsv, arm_elements.tsv),
splits the axis-to-copy distance into the stretch covered by bundle 15 inside the A core (bundle 15 lies in the
~20-kb tandem arrays, step 28b; inside the core it is the core's array) and the rest. Reads the core only (to the end of
the core). The stretch is measured in bp, not from tandem-repeats' copy counts (which are fractional or 0 on some
arms). Writes results/sv/copy_position.tsv and copy_position.json."""
import json
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

s = pd.read_csv('results/sv/arm_element_segments.tsv', sep='\t')
e = pd.read_csv('results/sv/arm_elements.tsv', sep='\t')
rows = []
for (src, side), g in s.groupby(['src', 'side']):
    a = e[(e.src == src) & (e.side == side)].iloc[0]
    g = g[(g.a0 >= -1) & (g.a1 <= a.a_end_kb + 1)]
    q = g[g.bid == 15]
    ok = len(q) > 0 and bool((q.a1 < a.copy_d_kb).all())
    rows.append(dict(src=src, side=side, type=a.type, copy_d_kb=a.copy_d_kb, measured=ok,
                     arr_kb=round(q.a1.max() - q.a0.min(), 2) if ok else np.nan,
                     arr_start_kb=round(q.a0.min(), 2) if ok else np.nan))
d = pd.DataFrame(rows)
d['rest_kb'] = (d.copy_d_kb - d.arr_kb).round(2)
d.to_csv('results/sv/copy_position.tsv', sep='\t', index=False)
out = {}
for side, (t0, t1) in (('cen', ('C0', 'C1')), ('tel', ('T1', 'T0'))):
    x = d[(d.side == side) & d.measured]
    A, B = x[x.type == t0], x[x.type == t1]
    o = dict(n_arms=int((d.side == side).sum()), n_measured=len(x), n=[len(A), len(B)], types=[t0, t1])
    for col in ('copy_d_kb', 'arr_kb', 'rest_kb'):
        o[col] = dict(median=[round(float(A[col].median()), 1), round(float(B[col].median()), 1)],
                      iqr=[[round(float(v), 1) for v in A[col].quantile([.25, .75])], [round(float(v), 1) for v in B[col].quantile([.25, .75])]],
                      p=float(mannwhitneyu(A[col], B[col]).pvalue))
    o['diff_copy'] = round(o['copy_d_kb']['median'][0] - o['copy_d_kb']['median'][1], 1)
    o['diff_arr'] = round(o['arr_kb']['median'][0] - o['arr_kb']['median'][1], 1)
    o['diff_rest'] = round(o['rest_kb']['median'][0] - o['rest_kb']['median'][1], 1)
    o['rho_copy_arr'] = round(float(spearmanr(x.copy_d_kb, x.arr_kb)[0]), 2)
    # the stretch comes in steps: commonest lengths (rounded to 1 kb)
    o['arr_modes'] = {str(int(k)): int(v) for k, v in x.arr_kb.round(0).value_counts().head(4).items()}
    out[side] = o
    print(side, json.dumps(o))
c = out['cen']
assert c['n_measured'] >= 0.95 * c['n_arms']
print('C0 minus C1, median kb: copy distance', c['diff_copy'], 'array stretch', c['diff_arr'], 'rest', c['diff_rest'])   # the text's reading is checked in step 40
json.dump(out, open('results/sv/copy_position.json', 'w'), indent=1)
