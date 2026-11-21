"""Pipeline step 04: each aligned copy's structural labels (position among the SMN copies of its window, strand on the
GRCh38-oriented window, NAIP form at its 3' end, downstream environment, arrangement type, ancestry group).
The label code is that of explore/structure_vs_sequence/02_compare.py, moved here so that the lineage steps can run
before the comparison step (07), which checks that it computes the same table. Writes results/sv/copy_labels.tsv."""
import pandas as pd
O = 'results/sv'
AFR = {'ACB', 'ASW', 'ESN', 'GWD', 'LWK', 'MSL', 'YRI', 'MKK'}
SUPER = {**{p: 'AFR' for p in AFR}, **{p: 'AMR' for p in ['PUR', 'PEL', 'CLM', 'MXL']},
         **{p: 'EAS' for p in ['JPT', 'KHV', 'CHS', 'CDX', 'CHB']}, **{p: 'EUR' for p in ['GBR', 'FIN', 'TSI', 'IBS']},
         **{p: 'SAS' for p in ['BEB', 'PJL', 'GIH', 'STU', 'ITU']}}
names = [l[1:].strip() for l in open(f'{O}/copies.aln.fa') if l.startswith('>')]
_t = pd.read_csv(f'{O}/copy_tips.tsv', sep='\t').set_index('tip').loc[names]
names = [n for n, c in zip(names, _t.copy_class) if c != 'del7_8_SMN2like']   # del7-8 rows left out (step 07)
tips = pd.read_csv(f'{O}/copy_tips.tsv', sep='\t').set_index('tip').loc[names].reset_index()
nb = pd.read_csv('tables/neighbour_copies.tsv', sep='\t').set_index('copy_id')
hc = pd.read_csv('tables/haplotype_copy_counts.tsv', sep='\t')
res = set(hc[hc.resolved].src) | {s for s in nb.index.str.split('|').str[0] if 'PanSN' in s}
pop = dict(zip(hc.src, hc.population))
hu = pd.read_csv('results/decomp/haplotype_units.tsv', sep='\t'); arr = dict(zip(hu.src, hu.arrangement))
labels = pd.DataFrame({'tip': names, 'copy_id': tips.copy_id, 'src': tips.src, 'cls': tips.copy_class.replace({'hyb_1_2': 'h12', 'hyb_2_1': 'h21', 'del7_8_SMN2like': 'd78'})})
def position(cid):
    if cid not in nb.index: return None
    r = nb.loc[cid]
    if not r.resolved: return None
    if r.n_in_window == 1: return 'single'
    return 'centromeric' if r.rank_from_cen == 1 else ('telomeric' if r.rank_from_cen == r.n_in_window else 'middle')
labels['position'] = [position(c) for c in labels.copy_id]
labels['strand'] = [nb.strand_grch38.get(c) for c in labels.copy_id]
labels['naip3'] = [nb.naip_3p.get(c) for c in labels.copy_id]
labels['env'] = [nb.environment.get(c) for c in labels.copy_id]
labels['arrangement'] = [arr.get(s) if s in res else None for s in labels.src]
labels['afr'] = [pop.get(s) in AFR for s in labels.src]
labels['superpop'] = [SUPER.get(pop.get(s), 'other') if s in pop else None for s in labels.src]
labels.to_csv(f'{O}/copy_labels.tsv', sep='\t', index=False)
print(len(labels), 'copies labelled')
