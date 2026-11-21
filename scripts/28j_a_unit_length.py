"""Pipeline step 28j: why the family-A repeat units differ in length when the core does not.

For each arm of the 270 A+ A- pair haplotypes: the length of its family-A unit (step 13, palindrome_arms.tsv),
the end of the core (element A of step 28b, a_end_kb: from the axis to the end of the first copy of the core's end),
the element that follows the core in the arm's layout, and the difference unit - core. A unit 'runs past the core'
when it ends beyond the core's end. Summarised by arm form (cenSMN type | telSMN type).
Writes results/sv/a_unit_length.tsv and results/sv/a_unit_length.json."""
import json
import pandas as pd

O = 'results/sv'
arms = pd.read_csv(f'{O}/palindrome_arms.tsv', sep='\t').set_index('src')
el = pd.read_csv(f'{O}/arm_elements.tsv', sep='\t')
hap = pd.read_csv(f'{O}/arm_element_haps.tsv', sep='\t').set_index('src')
ex = pd.read_csv(f'{O}/locus_diagram_examples.tsv', sep='\t')
el = el[el.src.isin(arms.index) & el.src.isin(hap.index)].copy()
el['unit_kb'] = [arms.loc[s, 'left_kb' if d == 'cen' else 'right_kb'] for s, d in zip(el.src, el.side)]
el['form'] = [hap.loc[s, 'cen_type'] + '|' + hap.loc[s, 'tel_type'] for s in el.src]
el['next_el'] = el.layout.str.split().str[1]
el['over_kb'] = (el.unit_kb - el.a_end_kb).round(2)
el['past_core'] = el.over_kb > 0
el[['src', 'side', 'form', 'type', 'unit_kb', 'a_end_kb', 'over_kb', 'past_core', 'next_el', 'layout']].to_csv(
    f'{O}/a_unit_length.tsv', sep='\t', index=False)
out = dict(n_haps=int(el.src.nunique()), core_all=dict(median=round(float(el.a_end_kb.median()), 1),
           q10=round(float(el.a_end_kb.quantile(.1)), 1), q90=round(float(el.a_end_kb.quantile(.9)), 1)), by_form={})
for (f, d), x in el.groupby(['form', 'side']):
    out['by_form'].setdefault(f, {})[d] = dict(
        n=len(x), unit_med=round(float(x.unit_kb.median()), 1), core_med=round(float(x.a_end_kb.median()), 1),
        over_med=round(float(x.over_kb.median()), 1), past=int(x.past_core.sum()),
        next_E=int((x.next_el == 'E').sum()), next_Xi=int((x.next_el == 'Xi').sum()))
# the unit runs past the core on both arms, on neither, or on one
w = el.pivot(index='src', columns='side', values='past_core')
w['form'] = [hap.loc[s, 'cen_type'] + '|' + hap.loc[s, 'tel_type'] for s in w.index]
out['both_neither'] = {f: dict(both=int((x.cen & x.tel).sum()), neither=int((~x.cen & ~x.tel).sum()),
                               tel_only=int((~x.cen & x.tel).sum()), cen_only=int((x.cen & ~x.tel).sum()), n=len(x))
                       for f, x in w.groupby('form')}
exs = {}
for r in ex.itertuples():
    s = r.src if hasattr(r, 'src') else r[1]
    x = el[el.src == s]
    exs[s] = {d: dict(unit_kb=round(float(y.unit_kb.iloc[0]), 1), core_kb=round(float(y.a_end_kb.iloc[0]), 1),
                      next_el=y.next_el.iloc[0]) for d, y in x.groupby('side')}
out['examples'] = exs
# pooled over both arms: C0 haplotypes (either telSMN type), C1|T0, and the telSMN / cenSMN arms of C1|T1
el['ctype'] = el.form.str[:2]
def pool(x):
    return dict(n_arms=len(x), unit_med=round(float(x.unit_kb.median()), 1), over_med=round(float(x.over_kb.median()), 1),
                over_q10=round(float(x.over_kb.quantile(.1)), 1), over_q90=round(float(x.over_kb.quantile(.9)), 1), past=int(x.past_core.sum()))
out['pooled'] = dict(C0=pool(el[el.ctype == 'C0']), C1T0=pool(el[el.form == 'C1|T0']),
                     C1T1_tel=pool(el[(el.form == 'C1|T1') & (el.side == 'tel')]), C1T1_cen=pool(el[(el.form == 'C1|T1') & (el.side == 'cen')]))
bn = out['both_neither']
out['pooled']['C0_haps'] = bn['C0|T0']['n'] + bn['C0|T1']['n']; out['pooled']['C0_both'] = bn['C0|T0']['both'] + bn['C0|T1']['both']
# the next element after the core: E on C0 arms and on the telSMN arm, Xi on the cenSMN arm of C1 haplotypes
out['next_after_core'] = {f'{k[0]}|{k[1]}': {e: int(v) for e, v in g.next_el.value_counts().items()} for k, g in el.groupby(['ctype', 'side'])}
json.dump(out, open(f'{O}/a_unit_length.json', 'w'), indent=1, sort_keys=True)
print(json.dumps(out, indent=1, sort_keys=True))
