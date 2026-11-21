"""Step 19d: the supported pure-SMN2 clade of 40+ copies in the gene-body ML tree, and the S1-2 SMN1 link.

Reads results/sv/mltree.treefile, copy_labels.tsv, ml_hap_labels.tsv.
Writes results/sv/smn2_s12_clade.json (counts) and smn2_s12_clade.tsv (one row per copy of the clade).
"""
import json, sys
import pandas as pd
sys.path.insert(0, 'scripts')
from smnlib import parse_newick
O = 'results/sv'
T = parse_newick(open(f'{O}/mltree.treefile').read().strip().rstrip(';'))
CL = pd.read_csv(f'{O}/copy_labels.tsv', sep='\t').set_index('tip')
HL = pd.read_csv(f'{O}/ml_hap_labels.tsv', sep='\t')
def tips(n): return [n['name']] if not n['kids'] else [x for k in n['kids'] for x in tips(k)]
cands = []
def walk(n):
    if not n['kids']: return
    tp = [x for x in tips(n) if x in CL.index]
    if n['sup'] is not None and n['sup'] >= 95 and len(tp) >= 40 and set(CL.loc[tp, 'cls']) == {'SMN2'}: cands.append((n['sup'], tp))
    for k in n['kids']: walk(k)
walk(T)
assert len(cands) == 1, 'exactly one supported pure-SMN2 clade of 40+ copies'
sup, cl = cands[0]
h_all = set(CL.loc[cl, 'src'])
h12 = set(HL[HL.seq_label.str.split('+').apply(lambda v: 'S1-2' in v)].src)
h = h_all & set(HL.src)   # haplotypes in the structure tree (resolved, with an SMN1 copy)
out = {'n_copies': len(cl), 'support': int(sup), 'haps_all': len(h_all), 'haps_in_tree': len(h),
       'haps_with_s1_2': len(h & h12), 'afr_copies': int((CL.loc[cl, 'superpop'] == 'AFR').sum()),
       'copies_on_2m_1m': int((CL.loc[cl, 'arrangement'] == '2- 1-').sum())}
json.dump(out, open(f'{O}/smn2_s12_clade.json', 'w'), indent=1)
t = CL.loc[cl, ['copy_id', 'src', 'superpop', 'arrangement']].copy(); t['src_has_s1_2'] = t.src.isin(h12)
t.sort_values('copy_id').to_csv(f'{O}/smn2_s12_clade.tsv', sep='\t', index=False)
print(out)
