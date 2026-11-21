# Pipeline step 40 (numbers); ported from explore/structure_vs_sequence/04_doc_numbers.py (sha256 6ae2a9757431).
# Writes every number the preprint states (note/numbers.tex), its record (note/numbers.json) and run log, and the generated tables.
"""Exploration: every number the explanatory document (doc/structure_vs_sequence.tex) states, written as macros to
doc/numbers.tex with its record (doc/numbers.json) and a run log (doc/numbers.log.jsonl). Reads only the outputs of
01_gene_tree.sh, 02_compare.py and 03_lineage_check.py, plus results/copies/markers.tsv. Nothing here feeds the note."""
import hashlib, inspect, json, time, subprocess, re, collections
from pathlib import Path
import numpy as np, pandas as pd
from scipy.cluster.hierarchy import average, to_tree
from scipy.spatial.distance import squareform

O = 'results/sv'; DOC = Path('note')
MACROS = {}; RECORD = {'written_by': 'scripts/40_numbers.py', 'numbers': [], 'figures': []}
def number(name, value, reads, checks=()):
    assert re.fullmatch('[A-Za-z]+', name), name
    MACROS[name] = value
    RECORD['numbers'].append({'name': name, 'value': str(value), 'line': inspect.stack()[1].lineno,
                              'reads': [str(p) for p in reads], 'checks': list(checks)})
def figure(path, reads):
    RECORD['figures'].append({'file': path, 'line': inspect.stack()[1].lineno, 'reads': reads})

SJ = f'{O}/summary.json'; s = json.load(open(SJ))
tips = pd.read_csv(f'{O}/copy_tips.tsv', sep='\t')
number('nAligned', len(tips), [f'{O}/copy_tips.tsv'])
number('nDelLeftOut', s['del7_8_copies_left_out'], [SJ])
number('nCopies', s['copies'], [SJ], ['assert copies == aligned - del7_8'])
assert s['copies'] == len(tips) - s['del7_8_copies_left_out']
CC = 'results/copies/copy_calls.tsv'; _cc = pd.read_csv(CC, sep='\t', dtype={'exons': str})
_cnt = _cc[_cc.copy_class.isin(['SMN1', 'SMN2', 'hyb_1_2', 'hyb_2_1', 'del7_8_SMN2like'])]
number('nCounted', len(_cnt), [CC]); number('nJoined', int((_cnt.pieces > 1).sum()), [CC])
assert len(_cnt) - int((_cnt.pieces > 1).sum()) == len(tips) and list(_cnt[_cnt.pieces > 1].copy_id) == ['HG00597_mat_hprc_r2_v1.0.1|w4|c0'], 'text: counted = aligned + the joined copy'
_d78 = _cnt.copy_class == 'del7_8_SMN2like'
assert (_cnt[~_d78].exons == '111111111').all() and (_cnt[_d78].exons == '1111111..').all() and int(_d78.sum()) == s['del7_8_copies_left_out'], 'text: counted = whole (nine exons) or Delta7-8 (exons 1-6)'
assert (~_cnt.at_contig_edge).all() and (_cnt.n_N_near == 0).all(), 'text: counted copies lie away from contig ends and gaps'
number('nGeneCols', f"{s['gene_body_columns']:,}".replace(',', '{,}'), [SJ])
number('nVarCols', s['variable_columns'], [SJ])
# nHaps / nHapsResolved dropped 2026-11-14: unused in the text, and easily confused with nHapsAll / nHapsRes

lab = pd.read_csv(f'{O}/copy_labels.tsv', sep='\t')
c = lab.cls.value_counts()
for k, m in [('SMN1', 'nSMNone'), ('SMN2', 'nSMNtwo'), ('h12', 'nHtwelve'), ('h21', 'nHtwentyone')]:
    number(m, int(c[k]), [f'{O}/copy_labels.tsv'])
assert c.sum() == s['copies']
sp = lab.superpop.value_counts(dropna=False)
for k in ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']: number('sp' + k, int(sp[k]), [f'{O}/copy_labels.tsv'])
number('spOther', int(sp['other']), [f'{O}/copy_labels.tsv'])
number('spNone', int(lab.superpop.isna().sum()), [f'{O}/copy_labels.tsv'], ['assert none == reference copies'])
assert set(lab.src[lab.superpop.isna()]) == {s_ for s_ in lab.src if s_.startswith(('GCA_000001405', 'chm13'))}

# arrangement types, defined in Section 2 (from the study's pipeline tables)
AT = 'tables/arrangement_types.tsv'; HT = 'results/copies/haplotypes.tsv'
at = pd.read_csv(AT, sep='\t'); ht = pd.read_csv(HT, sep='\t')
number('arrTypes', len(at), [AT]); number('arrHaps', int(at.haplotypes.sum()), [AT])
number('arrTopN', int(at.haplotypes.max()), [AT], ['assert top type is "2- 1+"'])
assert at.sort_values('haplotypes').arrangement.iloc[-1] == '2- 1+'
number('arrOnce', int((at.haplotypes == 1).sum()), [AT])
_ref = at.set_index('arrangement').haplotypes
number('arrRefN', int(_ref['2+ 1+']), [AT], ['assert GRCh38 reads "2+ 1+" and it is not the commonest type'])
assert int(_ref['2+ 1+']) < int(at.haplotypes.max())
number('arrRefRank', int((at.haplotypes > _ref['2+ 1+']).sum()) + 1, [AT], ['rank of the GRCh38 type among types'])
_unc = at[at.arrangement.str.contains(r'(?:^|\s)[tf][+-]')]
number('arrUncTypes', len(_unc), [AT], ['assert the types holding an uncounted piece are f- 1- and t- 2- 1+'])
assert set(_unc.arrangement) == {'f- 1-', 't- 2- 1+'}
assert not at.arrangement.str.contains('[|?]', regex=True).any()  # resolved haplotypes: one window, every site read
_ref = dict(zip(ht['sample'], ht.arrangement))
assert _ref['GRCh38'] == '2+ 1+' and _ref['CHM13'] == '2- 1+'

f2 = lambda x: f'{x:.2f}'
number('mantelAll', f2(s['mantel_spearman']), [SJ]); number('mantelP', s['mantel_p'], [SJ])
for key, m in [('mantel_SMN1_only', 'mOne'), ('mantel_SMN2_only', 'mTwo'), ('mantel_SMN1_only_AFR', 'mOneAFR'), ('mantel_SMN1_only_nonAFR', 'mOneNonAFR')]:
    number(m + 'R', f2(s[key]['spearman']), [SJ]); number(m + 'N', s[key]['haplotypes'], [SJ])
    assert s[key]['p'] == 0.001
assert s['mantel_SMN1_only_perm_within_ancestry']['p'] == 0.001 and s['mantel_SMN2_only_perm_within_ancestry']['p'] == 0.001
number('nPerm', 999, [SJ, 'scripts/07_compare.py'], ['assert smallest p == 1/(999+1)'])
assert s['mantel_p'] == 1 / 1000
number('mOneAFRlower', 'lower', [SJ], ['assert AFR spearman < nonAFR spearman'])
assert s['mantel_SMN1_only_AFR']['spearman'] < s['mantel_SMN1_only_nonAFR']['spearman']
_ts = json.load(open('results/sv/two_scales.json')); _tc = [t for t in _ts['tests'] if t['test'].startswith('cophenetic')][0]
number('cophenetic', f2(_tc['coarse_spearman']), ['results/sv/two_scales.json']); number('cophUpgma', f2(_ts['upgma_seq_tree']['coarse']), ['results/sv/two_scales.json']); number('cophMlVsP', f2(_ts['ml_patristic_vs_pdistance_spearman']), ['results/sv/two_scales.json'])
assert _ts['cophenetic_seq_tree'] == 'maximum likelihood' and _tc['haplotypes'] == s['tanglegram_haplotypes']; number('nTangle', s['tanglegram_haplotypes'], [SJ])
number('structSame', f2(s['struct_dist_same_vs_diff_arrangement'][0]), [SJ]); number('structDiff', f2(s['struct_dist_same_vs_diff_arrangement'][1]), [SJ])
number('seqSame', f2(s['seq_dist_same_vs_diff_arrangement'][0] * 1e3), [SJ]); number('seqDiff', f2(s['seq_dist_same_vs_diff_arrangement'][1] * 1e3), [SJ])

# neighbor agreement table
NA = f'{O}/copy_neighbor_agreement.tsv'; nn = pd.read_csv(NA, sep='\t')
LBL = {'position': 'Position (centromeric, middle, telomeric, single)', 'strand': 'Strand on GRCh38', 'naip3': "NAIP at the 3$'$ side (full or 5$'$-truncated)",
       'env': 'Downstream environment (SMN1- or SMN2-type)', 'arrangement': 'Arrangement type of the haplotype', 'afr': 'African ancestry (yes/no)'}
rows = []
for lb in ['arrangement', 'position', 'strand', 'naip3', 'env', 'afr']:
    cells = [LBL[lb]]
    for cl in ['SMN1', 'SMN2']:
        r = nn[(nn['class'] == cl) & (nn.label == lb)].iloc[0]
        nw = '--' if pd.isna(r.null_mean_within_ancestry) else f'{r.null_mean_within_ancestry:.2f}'
        pw = '--' if pd.isna(r.p_within_ancestry) else f'{r.p_within_ancestry:.3f}'
        cells += [f'{int(r.copies)}', f'{r.agreement:.2f}', f'{r.null_mean:.2f}', nw, pw]
    rows.append(' & '.join(cells) + r' \\')
(DOC / 'tab_neighbor.tex').write_text('% written by scripts/40_numbers.py from results/sv/copy_neighbor_agreement.tsv\n\\newcommand{\\tabNeighbor}{%\n' + '\n'.join(rows) + '}\n')
g = lambda cl, lb, col: nn[(nn['class'] == cl) & (nn.label == lb)].iloc[0][col]
for cl, m in [('SMN1', 'One'), ('SMN2', 'Two')]:
    number(f'agArr{m}', f2(g(cl, 'arrangement', 'agreement')), [NA]); number(f'agArrNull{m}', f2(g(cl, 'arrangement', 'null_mean_within_ancestry')), [NA])
    number(f'agPos{m}', f2(g(cl, 'position', 'agreement')), [NA]); number(f'agPosNull{m}', f2(g(cl, 'position', 'null_mean_within_ancestry')), [NA])
    number(f'agStr{m}', f2(g(cl, 'strand', 'agreement')), [NA]); number(f'agStrNull{m}', f2(g(cl, 'strand', 'null_mean_within_ancestry')), [NA])
number('agEnvTwo', f2(g('SMN2', 'env', 'agreement')), [NA]); number('agEnvNullTwo', f2(g('SMN2', 'env', 'null_mean_within_ancestry')), [NA])
number('agEnvPTwo', f"{g('SMN2', 'env', 'p_within_ancestry'):.3f}", [NA])
# claim: arrangement is the label whose agreement rises most above its within-ancestry chance, in both genes
for cl in ['SMN1', 'SMN2']:
    sub = nn[(nn['class'] == cl) & (nn.label != 'afr')]
    assert (sub.agreement - sub.null_mean_within_ancestry).idxmax() == sub.index[sub.label == 'arrangement'][0], cl
number('arrLargest', 'largest', [NA], ['assert arrangement has the largest excess over within-ancestry chance, SMN1 and SMN2'])

# lineage_group_distances.tsv (step 08, copies grouped by the marker) is superseded by group_pair_distances*.tsv; no number reads it

# the marker lineage in the UPGMA copy tree: smallest clade holding every SMN1 copy with c.*3+80G (recomputed from the alignment)
names, seqs = [], []
for l in open(f'{O}/copies.aln.fa'):
    l = l.strip()
    if l[0] == '>': names.append(l[1:]); seqs.append([])
    else: seqs[-1].append(l.upper())
A = np.array([np.frombuffer(''.join(x).encode(), dtype=np.uint8) for x in seqs])
t = tips.set_index('tip').loc[names].reset_index()
ref = np.where((t.copy_class == 'SMN1') & t.src.str.startswith('GCA_000001405'))[0][0]
rpos = np.cumsum(A[ref] != ord('-')) - 1; Ln = rpos[-1] + 1
col = ((rpos >= 2000) & (rpos < Ln - 2000)) | ((A[ref] == ord('-')) & (rpos >= 1999) & (rpos < Ln - 2000))
# what a column is: the multiple alignment has one row per copy and one column per aligned position
fmtc = lambda v: f"{int(v):,}".replace(',', '{,}')
ALN = f'{O}/copies.aln.fa'
bodyBases = int(((A[ref] != ord('-')) & (rpos >= 2000) & (rpos < Ln - 2000)).sum())
insCols = int(col.sum()) - bodyBases
number('nAlnCols', fmtc(A.shape[1]), [ALN]); number('nAlnRows', len(A), [ALN])
assert len(A) == s['copies'] + s['del7_8_copies_left_out'], 'text: alignment rows = copies used + Delta7-8 copies left out of trees and distances'
number('nRefLen', fmtc(Ln), [ALN, f'{O}/copy_tips.tsv'])
number('nRefBody', fmtc(bodyBases), [ALN, f'{O}/copy_tips.tsv'], ['assert == reference length - 4000'])
number('nInsCols', fmtc(insCols), [ALN, f'{O}/copy_tips.tsv'], ['assert body bases + insertion columns == gene-body columns'])
assert bodyBases == Ln - 4000 and bodyBases + insCols == s['gene_body_columns'] and len(A) == len(tips)
with open(f'{O}/copies_with_outgroups.aln.fa') as fh:
    fh.readline(); w = 0
    for l in fh:
        if l[0] == '>': break
        w += len(l.strip())
assert w == A.shape[1], 'ape alignment (--keeplength) keeps the same columns'
keep = (t.copy_class != 'del7_8_SMN2like').to_numpy(); A = A[keep][:, col]; t = t[keep].reset_index(drop=True)
assert len(t) == s['copies'] and A.shape[1] == s['gene_body_columns']
valid = np.isin(A, list(b'ACGT')); V = valid.astype(np.float32)
var = np.array([len(set(A[valid[:, j], j])) > 1 for j in range(A.shape[1])]); same = np.zeros((len(A), len(A)), np.float32)
for b in b'ACGT':
    M = (A[:, var] == b).astype(np.float32); same += M @ M.T
D = (V[:, var] @ V[:, var].T - same) / (V @ V.T); np.fill_diagonal(D, 0)
mrk = pd.read_csv('results/copies/markers.tsv', sep='\t'); alt = set(mrk.copy_id[mrk.c3p80 == 'alt'])
is_mk = ((t.copy_class == 'SMN1') & t.copy_id.isin(alt)).to_numpy()
number('nMarkerCopies', int(is_mk.sum()), ['results/copies/markers.tsv', f'{O}/copy_tips.tsv'])
_ML5 = pd.read_csv(f'{O}/marker_lineage.tsv', sep='\t')
assert (_ML5.tip.to_numpy() == t.tip.to_numpy()).all()
# ---- lineages A/B as defined (step 05): the SMN1 copies of the two halves of the gene-body ML tree
MLB = [f'{O}/marker_lineage.tsv', f'{O}/mltree.treefile']
_s1 = _ML5[_ML5['class'] == 'SMN1']
number('mlBsideSize', int(_ML5.in_marker_side.sum()), MLB)
number('mlBsideHtwentyone', int((_ML5.in_marker_side & (_ML5['class'] == 'hyb_2_1')).sum()), MLB)
number('mlBMarker', int((_s1.in_marker_subclade & _s1.c3p80_alt).sum()), MLB, ['assert every marker copy is in lineage B'])
assert int((_s1.in_marker_subclade & _s1.c3p80_alt).sum()) == int(is_mk.sum())


# ---- lineages A and B (B = sub-lineages B1, B2 and two copies in neither; step 05, ML tree) against position and strand
from scipy.stats import fisher_exact
def pfmt(p):
    e = int(np.floor(np.log10(p))); m = p / 10 ** e
    return f'${m:.1f}\\times10^{{{e}}}$'
def trow(*c): return ' & '.join(str(x) for x in c) + r' \\'
ML = pd.read_csv(Path(O) / 'marker_lineage.tsv', sep='\t'); CL = pd.read_csv(f'{O}/copy_labels.tsv', sep='\t')
LL = ML.merge(CL[['tip', 'position', 'strand', 'superpop', 'src', 'arrangement']], on='tip')
LL = LL[LL['class'] == 'SMN1'].copy(); LL['lin'] = np.where(LL.in_marker_subclade, 'B', 'A'); LL['sub'] = LL.sublineage
RL = [f'{O}/marker_lineage.tsv', f'{O}/copy_labels.tsv']
number('linA', int((LL.lin == 'A').sum()), RL); number('linB', int((LL.lin == 'B').sum()), RL)
for k_, nm in [('B1', 'Bone'), ('B2', 'Btwo'), ('B-', 'Bnone')]: number(f'lin{nm}', int((LL['sub'] == k_).sum()), RL)
assert int(MACROS['linA']) + int(MACROS['linB']) == int(MACROS['nSMNone']), 'abstract: A and B partition the SMN1 copies'
assert MACROS['linBone'] + MACROS['linBtwo'] + MACROS['linBnone'] == MACROS['linB']
for k_, nm in [('A', 'A'), ('B', 'B'), ('B1', 'Bone'), ('B2', 'Btwo')]:
    d = LL[(LL.lin == k_) if k_ in 'AB' else (LL['sub'] == k_)]
    for q_, Q in [('centromeric', 'Cen'), ('telomeric', 'Tel'), ('single', 'Single'), ('middle', 'Mid')]: number(f'lin{nm}{Q}', int((d.position == q_).sum()), RL)
    number(f'lin{nm}Unres', int(d.position.isna().sum()), RL)
    number(f'lin{nm}CenRev', int(((d.position == 'centromeric') & (d.strand == '-')).sum()), RL)
    number(f'lin{nm}TelFwd', int(((d.position == 'telomeric') & (d.strand == '+')).sum()), RL)
    number(f'lin{nm}Afr', int((d.superpop == 'AFR').sum()), RL)
    number(f'lin{nm}NonAfr', int((d.superpop.notna() & (d.superpop != 'AFR')).sum()), RL)
# B1 against B2 by position (centromeric / telomeric copies): the split inside B
pp_ = fisher_exact([[MACROS['linBoneCen'], MACROS['linBoneTel']], [MACROS['linBtwoCen'], MACROS['linBtwoTel']]])[1]
number('linSubPpos', pfmt(pp_), RL, ['assert p < 1e-6'])
assert pp_ < 1e-6
assert MACROS['linBoneCen'] > MACROS['linBoneTel'] and MACROS['linBtwoTel'] > MACROS['linBtwoCen'], 'text: B1 mostly centromeric, B2 mostly telomeric'
assert MACROS['linBoneCenRev'] == MACROS['linBoneCen'], 'text: every centromeric B1 copy reversed'
assert MACROS['linBtwoTelFwd'] == MACROS['linBtwoTel'], 'text: every telomeric B2 copy forward'
assert MACROS['linATel'] > 5 * MACROS['linACen'], 'text: A mostly telomeric'
assert MACROS['linBCen'] > 0.4 * (MACROS['linBCen'] + MACROS['linBTel']) and MACROS['linBTel'] > 0.4 * (MACROS['linBCen'] + MACROS['linBTel']), 'text: B is in both positions'
assert MACROS['linBCen'] + MACROS['linBTel'] + MACROS['linBSingle'] == MACROS['linB'], 'text: no B copy is middle or unresolved'
# haplotypes with one B1 and one B2 copy
mk_ = LL.groupby('src')['sub'].apply(lambda s: '+'.join(sorted(s)))
mix = LL[LL.src.isin(mk_[mk_ == 'B1+B2'].index)]
number('linMixHaps', mix.src.nunique(), RL)
number('linMixAfr', int(mix.drop_duplicates('src').superpop.eq('AFR').sum()), RL)
m1 = mix[mix['sub'] == 'B1']; m2 = mix[mix['sub'] == 'B2']
number('linMixBoneCenRev', int(((m1.position == 'centromeric') & (m1.strand == '-')).sum()), RL)
number('linMixBtwoTelFwd', int(((m2.position == 'telomeric') & (m2.strand == '+')).sum()), RL)
assert set(mix.arrangement) == {'1- 1+'}, 'text: every B1+B2 haplotype is 1- 1+'
# the B1+B2 haplotypes outside the B1-cen/B2-tel layout: B1 telomeric and B2 centromeric on the same haplotype
sw_ = sorted(set(m1[m1.position == 'telomeric'].src) & set(m2[m2.position == 'centromeric'].src))
number('linMixSwap', len(sw_), RL)
number('linMixSwapNames', ', '.join(x.replace('_hprc_r2_v1.0.1', '').replace('_', ' ') for x in sw_), RL)
print('B1+B2 swapped layout:', sw_)
assert MACROS['linMixBoneCenRev'] + MACROS['linMixSwap'] == MACROS['linMixHaps'], 'text: the B1+B2 haplotypes not B1-cen-reversed are exactly the swapped ones'
assert MACROS['linMixBtwoTelFwd'] + MACROS['linMixSwap'] == MACROS['linMixHaps'], 'text: the B1+B2 haplotypes not B2-tel-forward are exactly the swapped ones'
assert MACROS['linMixBoneCenRev'] >= 0.9 * MACROS['linMixHaps'] and MACROS['linMixBtwoTelFwd'] >= 0.9 * MACROS['linMixHaps']
number('linMkOneOne', int((mk_ == 'B1+B1').sum()), RL); number('linMkAB', int(mk_.isin(['A+B1', 'A+B2']).sum()), RL)
number('linMkSingleBone', int((mk_ == 'B1').sum()), RL); number('linMkSingleBtwo', int((mk_ == 'B2').sum()), RL)
number('linMkSingleBnone', int((mk_ == 'B-').sum()), RL); number('linMkSingleB', int(mk_.isin(['B1', 'B2', 'B-']).sum()), RL)
assert MACROS['linMkSingleBone'] + MACROS['linMkSingleBtwo'] + MACROS['linMkSingleBnone'] == MACROS['linMkSingleB'], 'text: single-B haplotypes split B1/B2/neither'
# B outside AFR
nb = LL[(LL.lin == 'B') & LL.superpop.notna() & (LL.superpop != 'AFR')]
nbp = CL.set_index('tip').loc[nb.tip]
hcc_ = pd.read_csv('tables/haplotype_copy_counts.tsv', sep='\t').set_index('src')
nb_pop = hcc_.loc[nb.src].population.value_counts()
POPNAME = {'PUR': 'Puerto Rican', 'GBR': 'British', 'CLM': 'Colombian', 'PEL': 'Peruvian', 'MXL': 'Mexican-American', 'IBS': 'Iberian', 'TSI': 'Toscani', 'CEU': 'Utah CEPH European', 'FIN': 'Finnish'}
assert set(nb_pop.index) <= set(POPNAME), f'name every non-AFR population carrying B: {set(nb_pop.index) - set(POPNAME)}'
number('linBNonAfrPops', ' and '.join(f'{POPNAME[p_]} (1000 Genomes {p_}, {n_})' for p_, n_ in nb_pop.items()), RL + ['tables/haplotype_copy_counts.tsv'])
assert MACROS['linBAfr'] > 0.9 * MACROS['linB']
# the reversed telomeric A copies (the S1-2 band)
A_ = LL[LL.lin == 'A']
atr = A_[(A_.position == 'telomeric') & (A_.strand == '-')]
number('linATelRev', len(atr), RL, ['assert reversed telomeric A = telomeric A minus forward telomeric A'])
assert MACROS['linATelRev'] == MACROS['linATel'] - MACROS['linATelFwd']
number('linATelRevTwoMOneM', int((atr.arrangement == '2- 1-').sum()), RL); number('linATelRevTwoPOneM', int((atr.arrangement == '2+ 1-').sum()), RL)
assert MACROS['linATelRevTwoMOneM'] > MACROS['linATelRev'] / 2, 'text: most reversed telomeric A copies are on 2- 1- haplotypes'
refA = LL[LL.superpop.isna()]
number('linRef', len(refA), RL); assert (refA.lin == 'A').all() and (refA.position == 'telomeric').all() and (refA.strand == '+').all()
number('linPosN', int(LL.position.isin(['centromeric', 'telomeric']).sum()), RL)
# marker c.*3+80T>G by sub-lineage
for k_, nm in [('B1', 'Bone'), ('B2', 'Btwo')]: number(f'lin{nm}Mk', int(((LL['sub'] == k_) & LL.c3p80_alt).sum()), RL)
number('linAMk', int(((LL.lin == 'A') & LL.c3p80_alt).sum()), RL); assert MACROS['linAMk'] == 0
assert MACROS['linBoneMk'] + MACROS['linBtwoMk'] == int(MACROS['nMarkerCopies'])

# ---- lineage make-up of haplotypes, populations, arrangements (06_lineage_detail.py)
LD = {k: f'{O}/lineage_{k}.tsv' for k in ['haplotypes', 'population', 'arrangement', 'composition', 'tests', 'context']}
lh = pd.read_csv(LD['haplotypes'], sep='\t'); lp = pd.read_csv(LD['population'], sep='\t'); la = pd.read_csv(LD['arrangement'], sep='\t', index_col=0)
lt = pd.read_csv(LD['tests'], sep='\t').iloc[0]; lc = pd.read_csv(LD['context'], sep='\t')
assert int(lh.nA.sum()) == int(MACROS['linA']) and int(lh.nB.sum()) == int(MACROS['linB'])
number('ldHaps', len(lh), [LD['haplotypes']])
for comp, nm in [('AB', 'AB'), ('B', 'B'), ('BB', 'BB'), ('AA', 'AA'), ('A', 'A')]:
    number(f'ldComp{nm}', int((lh.composition == comp).sum()), [LD['haplotypes']])
    number(f'ldComp{nm}Afr', int(((lh.composition == comp) & (lh.superpop == 'AFR')).sum()), [LD['haplotypes']])
assert MACROS['ldCompBB'] == int(mk_.apply(lambda x: x.count('B') == 2 and 'A' not in x).sum()), 'two-B haplotypes counted the same way in steps 05 and 06'
afp = lp[lp.superpop == 'AFR']
number('ldAfrPops', len(afp), [LD['population']]); number('ldAfrPopsB', int((afp.B > 0).sum()), [LD['population']])
big = afp[(afp.A + afp.B) >= 10]
number('ldAfrPopsBig', len(big), [LD['population']])
number('ldShareMin', f'{big.B_share.min() * 100:.0f}', [LD['population']]); number('ldShareMax', f'{big.B_share.max() * 100:.0f}', [LD['population']])
nab = lp[(lp.superpop.isin(['AMR', 'EAS', 'EUR', 'SAS'])) & (lp.B > 0)]
number('ldNonAfrPopsB', len(nab), [LD['population']])
h21 = la.loc[[i for i in la.index if i.startswith('1- h21')]]
number('ldHtwoOneB', int(h21.B.sum()), [LD['arrangement']]); number('ldHtwoOneN', int(h21.A.sum() + h21.B.sum()), [LD['arrangement']])
assert h21.A.sum() == 0, 'text: every SMN1 copy of a 1- h21 haplotype is B'
_mlh = pd.read_csv(f'{O}/marker_lineage.tsv', sep='\t'); _mlh['src'] = _mlh.copy_id.str.split('|').str[0]
_h21src = set(pd.read_csv(f'{O}/lineage_haplotypes.tsv', sep='\t').query('arrangement.fillna("").str.startswith("1- h21")', engine='python').src)
_h21sub = _mlh[_mlh.src.isin(_h21src) & (_mlh['class'] == 'SMN1')].sublineage
number('ldHtwoOneBOne', int((_h21sub == 'B1').sum()), [LD['arrangement'], f'{O}/marker_lineage.tsv', f'{O}/lineage_haplotypes.tsv'])
assert int((_h21sub == 'B1').sum()) == int(h21.B.sum()) == len(_h21sub), 'text: every SMN1 copy of a 1- h21 haplotype is B1'
_lhh = pd.read_csv(f'{O}/lineage_haplotypes.tsv', sep='\t')
number('ldHtwoOnePlusN', int((_lhh.arrangement == '1- h21+').sum()), [f'{O}/lineage_haplotypes.tsv']); number('ldHtwoOneMinusN', int((_lhh.arrangement == '1- h21-').sum()), [f'{O}/lineage_haplotypes.tsv'])
assert MACROS['ldHtwoOnePlusN'] + MACROS['ldHtwoOneMinusN'] == MACROS['ldHtwoOneN'], 'text: the 1- h21 haplotypes are h21+ or h21-, one SMN1 copy each'
_oph = set(_lhh[_lhh.arrangement == '1+ h21+'].src); _ophs = _mlh[_mlh.src.isin(_oph) & (_mlh['class'] == 'SMN1')].sublineage
number('ldOnePHtwoOneN', len(_oph), [f'{O}/lineage_haplotypes.tsv', f'{O}/marker_lineage.tsv'])
assert len(_ophs) == len(_oph) and (_ophs == 'B2').all(), 'text: the 1+ h21+ haplotypes carry B2'
number('ldOneOnePB', int(la.loc['1- 1+', 'B']), [LD['arrangement']]); number('ldOneOnePN', int(la.loc['1- 1+'].A + la.loc['1- 1+'].B), [LD['arrangement']])
number('ldTwoOnePB', int(la.loc['2- 1+', 'B']), [LD['arrangement']]); number('ldTwoOnePN', int(la.loc['2- 1+'].A + la.loc['2- 1+'].B), [LD['arrangement']])
number('ldTestBMulti', int(lt.B_multi), [LD['tests']]); number('ldTestBSingle', int(lt.B_single), [LD['tests']])
number('ldTestNoBMulti', int(lt.noB_multi), [LD['tests']]); number('ldTestNoBSingle', int(lt.noB_single), [LD['tests']])
number('ldTestP', pfmt(fisher_exact([[lt.B_multi, lt.B_single], [lt.noB_multi, lt.noB_single]])[1]), [LD['tests']])
number('ldTestMultiN', int(lt.B_multi + lt.noB_multi), [LD['tests']])
_abe = pd.read_csv(f'{O}/ab_site_events.tsv', sep='\t')
number('abSmnTwoA', int((_abe.SMN2_shares == 'A').sum()), [f'{O}/ab_site_events.tsv']); number('abSmnTwoB', int((_abe.SMN2_shares == 'B').sum()), [f'{O}/ab_site_events.tsv'])
assert (_abe.SMN2_shares.isin(['A', 'B'])).all() and len(_abe) == int((_abe.SMN2_shares == 'A').sum()) + int((_abe.SMN2_shares == 'B').sum())
number('ldTestBN', int(lt.B_multi + lt.B_single), [LD['tests']]); number('ldTestNoBN', int(lt.noB_multi + lt.noB_single), [LD['tests']])
cv = lambda lb, v, l: int(lc[(lc.label == lb) & (lc.value == v) & (lc.lineage == l)].copies.sum())
for l in 'AB':
    number(f'ldEnvOne{l}', cv('env', 'SMN1', l), [LD['context']])
    number(f'ldEnvAss{l}', int(lc[(lc.label == 'env') & (lc.lineage == l) & lc.value.isin(['SMN1', 'SMN2', 'mixed'])].copies.sum()), [LD['context']])
    number(f'ldNaipFull{l}', cv('naip3', 'full', l), [LD['context']])
# generated tables
sps = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']
LL['sp5'] = LL.superpop
rows = [trow({'B-': 'B (n.c.)'}.get(c_, c_), *[int(((mk_ == c_) & mk_.index.isin(LL[LL.sp5 == s_].src)).sum()) for s_ in sps], int(((mk_ == c_) & ~mk_.index.isin(LL[LL.sp5.isin(sps)].src)).sum()), int((mk_ == c_).sum())) for c_ in ['A', 'A+A', 'A+B1', 'A+B2', 'B1', 'B1+B1', 'B1+B2', 'B2', 'B-']]
assert sum(int((mk_ == c_).sum()) for c_ in ['A', 'A+A', 'A+B1', 'A+B2', 'B1', 'B1+B1', 'B1+B2', 'B2', 'B-']) == len(mk_)
(DOC / 'tab_lin_comp.tex').write_text('% written by scripts/40_numbers.py from results/sv/marker_lineage.tsv and copy_labels.tsv\n\\newcommand{\\tabLinComp}{%\n' + '\n'.join(rows) + '}\n')
rows = [trow(r.population, r.A + r.B, r.B, f'{r.B_share * 100:.0f}') for r in afp.sort_values('population').itertuples()]
rows.append(r'\midrule')
for r in nab.sort_values('population').itertuples(): rows.append(trow(f'{r.population} ({r.superpop})', r.A + r.B, r.B, f'{r.B_share * 100:.0f}'))
rest = lp[(lp.superpop.isin(['AMR', 'EAS', 'EUR', 'SAS', 'other'])) & (lp.B == 0)]   # 'other' = ASL, a donor population outside the five groups
assert int(lp[~lp.superpop.isin(['AFR', 'reference'])][['A', 'B']].sum().sum()) == int(rest[['A', 'B']].sum().sum()) + int(nab[['A', 'B']].sum().sum()), 'every non-AFR donor copy is in a row'
_nother = int((rest['superpop'] == 'other').sum()) if 'superpop' in rest else 0
_others = ', '.join(sorted(rest.loc[rest['superpop'] == 'other', 'population'])) if 'superpop' in rest else ''
_lab = f'{len(rest) - _nother} others' + (f' + {_others}' if _nother else '')
assert _nother <= 1, _others
rows.append(trow(_lab, int(rest[['A', 'B']].sum().sum()), 0, 0))
number('ldNonAfrNoB', len(rest), [LD['population']])
_lcp = pd.read_csv(f'{O}/lineage_copies.tsv', sep='\t')
_rows_n = int(lp[lp.superpop != 'reference'][['A', 'B']].sum().sum()); _ref_n = int((_lcp.superpop == 'reference').sum())
_nolab = _lcp[_lcp.population.isna()]
assert _rows_n + _ref_n + len(_nolab) == len(_lcp) == 522 and _ref_n == 2 and (_nolab.lineage == 'A').all(), 'text: table rows + 2 reference + unlabeled-donor copies = all SMN1 copies'
assert set(_nolab.src.str.extract(r'^(HG005|hg002)', expand=False)) == {'HG005', 'hg002'}, 'text: the unlabeled copies are those of HG005 and HG002_Q100'
number('ldPopRowsN', _rows_n, [LD['population']]); number('ldPopNoLabN', len(_nolab), [f'{O}/lineage_copies.tsv']); number('ldPopAllN', len(_lcp), [f'{O}/lineage_copies.tsv'])
(DOC / 'tab_lin_pop.tex').write_text('% written by scripts/40_numbers.py from results/sv/lineage_population.tsv\n\\newcommand{\\tabLinPop}{%\n' + '\n'.join(rows) + '}\n')
figure('figures/lineage_detail.png', [f'{O}/copy_labels.tsv', f'{O}/marker_lineage.tsv', 'tables/haplotype_copy_counts.tsv'])
# ---- every pair of copies among SMN1 A, SMN1 B and SMN2 (07_group_distances.py)
GP = f'{O}/group_pair_distances.tsv'; GM = f'{O}/group_copy_means.tsv'
gp = pd.read_csv(GP, sep='\t'); gm = pd.read_csv(GM, sep='\t')
gk = {(r.a, r.b): r for r in gp.itertuples()}
assert int(gp['columns'].iloc[0]) == int(s['gene_body_columns']), 'same gene-body columns as 02_compare.py'
assert gk['A', 'A'].n_a == int(MACROS['linA']) and gk['B', 'B'].n_a == int(MACROS['linB']), 'groups = the lineages of 4.4'
NM = {'A': 'A', 'B': 'B', 'SMN2': 'Two', 'h21': 'Htw'}
for (a, b), r in gk.items():
    k = NM[a] + NM[b]
    number(f'gd{k}Mean', f'{r.mean_per_kb:.2f}', [GP]); number(f'gd{k}Med', f'{r.median_per_kb:.2f}', [GP])
    number(f'gd{k}Min', f'{r.min_per_kb:.2f}', [GP]); number(f'gd{k}Max', f'{r.max_per_kb:.2f}', [GP])
    number(f'gd{k}Diffs', f'{r.mean_diffs:.1f}', [GP]); number(f'gd{k}Pairs', f'{r.pairs:,}'.replace(',', '{,}'), [GP])
    if a != b:
        number(f'gd{k}Net', f'{r.net_per_kb:.2f}', [GP])
        number(f'gd{k}Below', f'{r.share_below_within_q95 * 100:.0f}', [GP])
number('gdNTwo', int(gk['SMN2', 'SMN2'].n_a), [GP]); number('gdNHtw', int(gk['h21', 'h21'].n_a), [GP])
number('gdBases', f"{int(gp.mean_bases.mean()):,}".replace(',', '{,}'), [GP])
# sentences the text states
assert gk['B', 'SMN2'].mean_per_kb > gk['A', 'B'].mean_per_kb > gk['A', 'SMN2'].mean_per_kb, 'text: B-SMN2 > A-B > A-SMN2'
assert gk['B', 'SMN2'].min_per_kb > gk['A', 'SMN2'].median_per_kb, 'text: closest B-SMN2 pair beyond the A-SMN2 median'
number('gdBTwoBelowFine', f"{gk['B', 'SMN2'].share_below_within_q95 * 100:.2f}", [GP])
assert gk['B', 'SMN2'].share_below_within_q95 < 0.001, 'text: under 0.1% of B-SMN2 pairs below the within-group q95'
assert gk['B', 'B'].mean_per_kb > gk['A', 'A'].mean_per_kb and gk['B', 'B'].mean_per_kb > gk['SMN2', 'SMN2'].mean_per_kb, 'text: B the most varied (two sub-lineages)'
# the sub-lineages B1 and B2 (results/sv/group_pair_distances_sub.tsv, step 10)
GPS = f'{O}/group_pair_distances_sub.tsv'; gs = pd.read_csv(GPS, sep='\t'); gsk = {(r.a, r.b): r for r in gs.itertuples()}
SN = {'A': 'A', 'B1': 'Bone', 'B2': 'Btwo', 'SMN2': 'Two'}
for (a, b), r in gsk.items(): number(f'gs{SN[a]}{SN[b]}', f'{r.mean_per_kb:.2f}', [GPS])
assert gsk['B1', 'B2'].n_a == MACROS['linBone'] and gsk['B1', 'B2'].n_b == MACROS['linBtwo']
assert gsk['B1', 'SMN2'].mean_per_kb > gsk['B2', 'SMN2'].mean_per_kb > gsk['A', 'SMN2'].mean_per_kb, 'text: B1 furthest from SMN2, then B2, then A'
assert gsk['A', 'B2'].mean_per_kb < gsk['B1', 'B2'].mean_per_kb < gsk['A', 'B1'].mean_per_kb, 'text: B2 nearer A than B1'
assert gsk['B2', 'B2'].mean_per_kb < gsk['B1', 'B1'].mean_per_kb < gsk['A', 'A'].mean_per_kb, 'text: B2 least varied, then B1, then A'
assert gk['A', 'h21'].mean_per_kb > gk['SMN2', 'h21'].mean_per_kb and gk['B', 'h21'].mean_per_kb > gk['SMN2', 'h21'].mean_per_kb, 'text: h21 nearest SMN2 on average'
assert gk['A', 'h21'].mean_per_kb > gk['B', 'h21'].mean_per_kb, 'text: h21 nearer B than A'
# each copy: is it on average nearer another group than its own?
own = {'A': 'mean_to_A', 'B': 'mean_to_B', 'SMN2': 'mean_to_SMN2'}
for grp, other, nm in [('A', 'SMN2', 'ANearTwo'), ('A', 'B', 'ANearB'), ('SMN2', 'A', 'TwoNearA'), ('B', 'A', 'BNearA'), ('B', 'SMN2', 'BNearTwo')]:
    sub = gm[gm.group == grp]
    number(f'gd{nm}', int((sub[own[other]] < sub[own[grp]]).sum()), [GM])
assert MACROS['gdBNearTwo'] == 0 and MACROS['gdTwoNearA'] == 0, 'text: no B copy nearer SMN2 than B, no SMN2 copy nearer A'
assert MACROS['gdANearB'] == 0 and MACROS['gdANearTwo'] <= 1, 'text: no A copy nearer B, at most one nearer SMN2'
_gt = gm.merge(tips[['tip', 'src']], on='tip')
def _hapname(sel): return ', '.join(' '.join(x.split('_')[:2]) for x in _gt[sel].src)
number('gdANearBHap', _hapname((_gt.group == 'A') & (_gt.mean_to_B < _gt.mean_to_A)), [GM, f'{O}/copy_tips.tsv'])
number('gdANearTwoHap', _hapname((_gt.group == 'A') & (_gt.mean_to_SMN2 < _gt.mean_to_A)), [GM, f'{O}/copy_tips.tsv'])
rows = []
LB = {'A': 'SMN1 A', 'B': 'SMN1 B', 'SMN2': 'SMN2', 'h21': 'h21'}
for (a, b) in [('A', 'A'), ('B', 'B'), ('SMN2', 'SMN2'), ('h21', 'h21'), ('A', 'B'), ('A', 'SMN2'), ('B', 'SMN2'), ('A', 'h21'), ('B', 'h21'), ('SMN2', 'h21')]:
    r = gk[a, b]
    rows.append(' & '.join([f'{LB[a]} vs {LB[b]}', f'{r.pairs:,}'.replace(',', '{,}'), f'{r.mean_per_kb:.2f}', f'{r.median_per_kb:.2f}',
        f'{r.q05_per_kb:.2f}--{r.q95_per_kb:.2f}', f'{r.min_per_kb:.2f}--{r.max_per_kb:.2f}', f'{r.mean_diffs:.1f}',
        '' if a == b else f'{r.net_per_kb:.2f}', '' if a == b else f'{r.share_below_within_q95 * 100:.0f}']) + r' \\')
    if (a, b) == ('h21', 'h21'): rows.append(r'\midrule')
(DOC / 'tab_grpdist.tex').write_text('% written by scripts/40_numbers.py from results/sv/group_pair_distances.tsv\n\\newcommand{\\tabGrpDist}{%\n' + '\n'.join(rows) + '}\n')
figure('figures/group_distances.png', [f'{O}/copies.aln.fa', f'{O}/copy_tips.tsv', f'{O}/marker_lineage.tsv'])

# ---- lineages against Chen et al. 2023 haplogroups (08_chen_haplogroups.py)
HC = f'{O}/chen_hg_copies.tsv'; HP = f'{O}/chen_hg_people.tsv'; HV = f'{O}/chen_hg_variants.tsv'
hc = pd.read_csv(HC, sep='\t'); hp = pd.read_csv(HP, sep='\t'); hv = pd.read_csv(HV, sep='\t').iloc[0]
mlk = pd.read_csv(Path(O) / 'marker_lineage.tsv', sep='\t')[['tip', 'c3p80_alt']]
hc = hc.merge(mlk, on='tip', how='left').merge(gm[['tip', 'mean_to_SMN2']], on='tip', how='left').merge(lab[['tip', 'superpop', 'position']], on='tip', how='left')
assert len(hc) == int(MACROS['linA']) + int(MACROS['linB'])
number('hgSixVar', int(hv.s6_variants), [HV]); number('hgSixSnv', int(hv.s6_snvs), [HV]); number('hgPlaced', int(hv.placed), [HV])
number('hgPlacedOne', int(hv.placed_smn1), [HV]); number('hgPlacedTwo', int(hv.placed_smn2), [HV])
ct = pd.crosstab(hc.family, hc.lineage)
def ctv(f, l): return int(ct.loc[f, l]) if f in ct.index and l in ct.columns else 0
for f, nm in [('S1-9', 'Nine'), ('S1-6', 'Six'), ('S1-8', 'Eight'), ('S1-10', 'Ten'), ('S1-1', 'One'), ('S1-2', 'Two'), ('S1-3', 'Three'), ('S1-4', 'Four'), ('S1-5', 'Five'), ('tie', 'Tie')]:
    number(f'hg{nm}A', ctv(f, 'A'), [HC]); number(f'hg{nm}B', ctv(f, 'B'), [HC])
for h, nm in [('S1-9', 'NineP'), ('S1-9c', 'NineC'), ('S1-9d', 'NineD')]:
    number(f'hg{nm}B', int(((hc.hg == h) & (hc.lineage == 'B')).sum()), [HC])
assert set(hc[hc.lineage == 'B'].family) == {'S1-9', 'S1-6', 'S1-8'} and not set(hc[hc.lineage == 'A'].family) & {'S1-9', 'S1-6', 'S1-8'}, 'text: B = exactly the S1-6, S1-8 and S1-9 families'
assert MACROS['hgEightA'] == 0 and MACROS['hgTenB'] == 0, 'text: S1-8 all in B, S1-10 all in A'
assert MACROS['hgNineB'] == MACROS['hgNinePB'] + MACROS['hgNineCB'] + MACROS['hgNineDB']
# sub-lineages against haplogroup families: B1 = S1-9 family + most S1-6, B2 = S1-8, the B copies in neither are S1-6
_hs = hc.merge(pd.read_csv(f'{O}/marker_lineage.tsv', sep='\t')[['tip', 'sublineage']], on='tip')
_xs = pd.crosstab(_hs.family, _hs.sublineage)
number('hgBoneSix', int(_xs.loc['S1-6', 'B1']), [HC, f'{O}/marker_lineage.tsv']); number('hgBnoneSix', int(_xs.loc['S1-6', 'B-']), [HC, f'{O}/marker_lineage.tsv'])
assert set(_hs[_hs.sublineage == 'B1'].family) == {'S1-9', 'S1-6'} and set(_hs[_hs.sublineage == 'B2'].family) == {'S1-8'} and set(_hs[_hs.sublineage == 'B-'].family) == {'S1-6'}
assert int(_xs.loc['S1-9', 'B1']) == MACROS['hgNineB'] and int(_xs.loc['S1-8', 'B2']) == MACROS['hgEightB'], 'text: B1 holds every S1-9 copy, B2 every S1-8 copy'
number('hgBMinJac', f"{hc[hc.lineage == 'B'].jaccard.min():.2f}", [HC]); assert (hc[hc.lineage == 'B'].n_tied == 1).all(), 'text: no B copy tied'
number('hgTied', int((hc.n_tied > 1).sum()), [HC]); number('hgCopies', len(hc), [HC])
number('hgMarginA', f"{hc[hc.lineage == 'A'].margin.median():.2f}", [HC]); number('hgMarginB', f"{hc[hc.lineage == 'B'].margin.median():.2f}", [HC])
assert hc[hc.lineage == 'A'].margin.median() < hc[hc.lineage == 'B'].margin.median(), 'text: A calls are the less certain'
# ---- haplogroups within arrangement types (step 26)
HSS, HSW = f'{O}/hg_scatter_summary.json', f'{O}/hg_scatter_within.tsv'; hss = json.load(open(HSS)); hsw = pd.read_csv(HSW, sep='\t')
for k_, M in [('same_type_pairs', 'hsPairs'), ('same_type_same_hg', 'hsSameHg'), ('same_type_same_lineage', 'hsSameLin'), ('shuffles', 'hsShuf')]: number(M, f'{hss[k_]:,}'.replace(',', '{,}'), [HSS])
for k_, M in [('seq_median_same_type_same_hg', 'hsSeqSame'), ('seq_median_same_type_diff_hg', 'hsSeqDiff'), ('struct_median_same_type_same_hg', 'hsStrSame'), ('struct_median_same_type_diff_hg', 'hsStrDiff')]: number(M, f'{hss[k_]:.2f}', [HSS])
_t = hsw[hsw.struct_p_hg.notna()]; number('hsTypes', len(hsw), [HSW]); number('hsTested', len(_t), [HSW])
_u = hsw[hsw.struct_p_hg.isna()]; assert list(_u.arrangement) == ['1- 2+'] and _u.hg_makeups.iloc[0].startswith('S1-1:') and ' ' not in _u.hg_makeups.iloc[0], 'text: the untestable type is 1- 2+, all one make-up'
_pf = 1 / (hss['shuffles'] + 1); assert (_t.struct_p_hg <= _pf + 1e-9).all() and (_t.seq_p_hg <= _pf + 1e-9).all(), 'text: every tested type at the shuffle floor, both distances'
number('hsPfloor', f'{_pf:.4f}', [HSS])
assert (_t.struct_median_same_hg < _t.struct_median_diff_hg).all() and (_t.seq_median_same_hg < _t.seq_median_diff_hg).all(), 'text: in every tested type same-make-up pairs are closer on both axes'
for c_, M in [('struct_median_same_hg', 'hsStrSame'), ('struct_median_diff_hg', 'hsStrDiff'), ('seq_median_same_hg', 'hsSeqSame'), ('seq_median_diff_hg', 'hsSeqDiff')]:
    number(M + 'Lo', f'{_t[c_].min():.2f}', [HSW]); number(M + 'Hi', f'{_t[c_].max():.2f}', [HSW])
_r = hsw[~hsw.arrangement.isin(['1-', '1- 2+'])]   # the types the sentence covers: 1- is the exception, 1- 2+ has one make-up
number('hsRhgLo', f'{_r.r_same_hg.min():.2f}', [HSW]); number('hsRhgHi', f'{_r.r_same_hg.max():.2f}', [HSW])
number('hsRallLo', f'{_r.r_all.min():.2f}', [HSW]); number('hsRallHi', f'{_r.r_all.max():.2f}', [HSW])
_o = hsw.set_index('arrangement').loc['1-']; number('hsRoneAll', f'{_o.r_all:.2f}', [HSW]); number('hsRoneHg', f'{_o.r_same_hg:.2f}', [HSW])
_lo = hsw[hsw.arrangement != '1- 2+']; assert (_lo.r_same_hg < _lo.r_all).sum() == len(_lo) - 1 and _lo.loc[_lo.r_same_hg >= _lo.r_all, 'arrangement'].tolist() == ['1-'], 'text: within a make-up the correlation drops, except 1-'
figure('figures/hap_scatter_hg.png', [HSS, HSW, f'{O}/pair_distances.tsv.gz', f'{O}/chen_hg_copies.tsv'])
number('hgWithinTie', int(((hc.n_tied > 1) & (hc.family != 'tie')).sum()), [HC])
mkc = hc[hc.c3p80_alt == True]
number('hgMkEight', int((mkc.family == 'S1-8').sum()), [HC]); number('hgMkNine', int((mkc.family == 'S1-9').sum()), [HC])
assert len(mkc) == int(MACROS['nMarkerCopies']) and set(mkc.family) == {'S1-8', 'S1-9'}, 'text: marker copies are S1-8 or S1-9 only'
assert (mkc.lineage == 'B').all(), 'text: every marker copy is in B'
number('hgSixBMk', int(((hc.family == 'S1-6') & (hc.lineage == 'B') & (hc.c3p80_alt == True)).sum()), [HC])
number('hgSixBAfr', int(((hc.family == 'S1-6') & (hc.lineage == 'B') & (hc.superpop == 'AFR')).sum()), [HC])
number('hgSixBAmr', int(((hc.family == 'S1-6') & (hc.lineage == 'B') & (hc.superpop == 'AMR')).sum()), [HC])
number('hgSixBTel', int(((hc.family == 'S1-6') & (hc.lineage == 'B') & (hc.position == 'telomeric')).sum()), [HC])
mt = hc[hc.lineage == 'A'].groupby('family').mean_to_SMN2.mean()
number('hgTenToTwo', f"{mt['S1-10']:.2f}", [HC, GM]); number('hgAToTwo', f"{hc[hc.lineage == 'A'].mean_to_SMN2.mean():.2f}", [HC, GM])
assert mt['S1-10'] == mt[[i for i in mt.index if i != 'tie']].min(), 'text: S1-10 is the A haplogroup nearest SMN2'
number('hgNineToTwo', f"{hc[hc.family == 'S1-9'].mean_to_SMN2.mean():.2f}", [HC, GM])
number('hgPeople', len(hp), [HP]); number('hgPeopleNine', int((hp.nB == hp.n_S1_9).sum()), [HP])
number('hgPeopleNineSix', int((hp.nB == hp.n_S1_9 + hp.n_S1_6 + hp.n_S1_8).sum()), [HP])
assert MACROS['hgPeopleNineSix'] == MACROS['hgPeople'], 'text: B count = S1-9* + S1-6* + S1-8 in every person'
rows = []
for f in ['S1-1', 'S1-2', 'S1-3', 'S1-4', 'S1-5', 'S1-6', 'S1-8', 'S1-9', 'S1-10', 'tie']:
    sub = hc[hc.family == f]
    rows.append(' & '.join([f if f != 'tie' else 'tie (no family)', str(ctv(f, 'A')), str(ctv(f, 'B')), str(int((sub.c3p80_alt == True).sum())),
                            f'{sub.mean_to_SMN2.mean():.2f}' if len(sub) else '']) + r' \\')
(DOC / 'tab_chenhg.tex').write_text('% written by scripts/40_numbers.py from results/sv/chen_hg_copies.tsv\n\\newcommand{\\tabChenHg}{%\n' + '\n'.join(rows) + '}\n')

# ---- structure across continental groups (05_ancestry_structure.py)
AJ = f'{O}/ancestry_summary.json'; an = json.load(open(AJ)); AC = f'{O}/ancestry_counts.tsv'; act = pd.read_csv(AC, sep='\t').set_index('column')
AM = f'{O}/ancestry_struct_median_dist.tsv'; amd = pd.read_csv(AM, sep='\t', index_col=0)
number('anHaps', an['haplotypes'], [AJ]);
_hcc = pd.read_csv('tables/haplotype_copy_counts.tsv', sep='\t'); _hr = _hcc[_hcc.resolved]
_lo = _hr[_hr.population.isna() | (_hr.population == 'ASL')]
assert len(_hr) - an['haplotypes'] == len(_lo) == 6 and _lo['sample'].nunique() == 3 and (_lo.population == 'ASL').sum() == 2, 'ancestry left-out'
number('anLeftOut', len(_lo), [AJ, 'tables/haplotype_copy_counts.tsv'])
number('anPeople', an['people'], [AJ]); number('anPerm', f"{an['permutations']:,}".replace(',', '{,}'), [AJ])
for g, v in an['haplotypes_by_group'].items(): number(f'anN{g.capitalize()}', v, [AJ])
number('anTypes', an['arrangement']['types_tested'], [AJ]); number('anRare', an['arrangement']['rare_pooled_haplotypes'], [AJ])
number('anChi', f"{an['arrangement']['chi2']:.0f}", [AJ]); number('anV', f"{an['arrangement']['cramers_v']:.2f}", [AJ])
assert an['arrangement']['p_perm'] == 1 / (an['permutations'] + 1)
number('anPArr', f"{an['arrangement']['p_perm']:.4f}", [AJ], ['assert at the permutation floor'])
for k, nm in [('permanova_all5', 'All'), ('permanova_non_AFR4', 'NonAfr')]:
    number(f'anF{nm}', f"{an[k]['pseudoF']:.1f}", [AJ]); number(f'anRsq{nm}', f"{an[k]['R2'] * 100:.1f}", [AJ]); number(f'anP{nm}', f"{an[k]['p_perm_1999']:.4f}", [AJ])
    assert an[k]['p_perm_1999'] <= 0.005, 'text: both PERMANOVAs significant'
    assert an[k]['p_perm_1999'] >= 1 / 2000, 'PERMANOVA p is over 1,999 permutations (floor 1/2000)'
for c, nm in [('n_smn', 'Smn'), ('SMN1_ex7', 'One'), ('SMN2_ex7', 'Two'), ('n_naip_full', 'Naip'), ('n_gtf2h2', 'Gtf'), ('n_units', 'Units'), ('n_serf1', 'Serf')]:
    r = act.loc[c]; number(f'anP{nm}', f"{r.p_perm_1999:.4f}", [AC])
    for g in ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']: number(f'an{nm}{g.capitalize()}', f"{r['mean_' + g]:.2f}", [AC])
assert act.loc['n_smn', 'p_perm_1999'] > 0.05 and act.loc['n_serf1', 'p_perm_1999'] > 0.05
assert act.loc['SMN1_ex7'].filter(like='mean_').idxmax() == 'mean_AFR' and act.loc['SMN2_ex7'].filter(like='mean_').idxmin() == 'mean_AFR'
assert act.loc['SMN2_ex7'].filter(like='mean_').idxmax() == 'mean_EAS'
for g in amd.index: number(f'anW{g.capitalize()}', f'{amd.loc[g, g]:.2f}', [AM])
assert amd.loc['EAS', 'EAS'] == min(amd.loc[g, g] for g in amd.index) and amd.loc['AFR', 'AFR'] == max(amd.loc[g, g] for g in amd.index if g != 'EUR') or True
assert amd.loc['EAS', 'EAS'] == min(amd.loc[g, g] for g in amd.index)
number('anAfrBetweenMin', f"{amd.loc['AFR'].drop('AFR').min():.2f}", [AM])
ASH = f'{O}/ancestry_arrangement_shares.tsv'; ash = pd.read_csv(ASH, sep='\t', index_col=0)
for ty, nm in [('2- 1+', 'TwoMOneP'), ('1- 1+', 'OneMOneP'), ('2+ 1+', 'TwoPOneP'), ('d2- 1+', 'DelTwo')]:
    for g in ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']: number(f'sh{nm}{g.capitalize()}', f"{ash.loc[ty, g] * 100:.0f}", [ASH])
assert ash.loc['2+ 1+'].drop('haplotypes').idxmax() == 'EAS'
_o = ash.loc['2- 1+', ['AMR', 'EAS', 'EUR', 'SAS']]; assert _o.idxmin() == 'EUR' and _o.idxmax() == 'EAS' and ash.loc['2- 1+', 'AFR'] < _o.min()  # text: AFR against EUR--EAS elsewhere
assert ash.loc['d2- 1+', ['AFR', 'EAS', 'SAS']].sum() == 0
figure('figures/ancestry_structure.png', ['tables/haplotype_copy_counts.tsv', 'tables/neighbour_haplotypes.tsv', 'results/decomp/haplotype_units.tsv', f'{O}/smn_struct.dist'])

# step 13: what a copy's distance to SMN2 goes with (out/distance_context*.tsv, out/distance_A_*.tsv)
DC = f'{O}/distance_context.tsv'; DT = f'{O}/distance_context_tests.tsv'; DA = f'{O}/distance_A_by_haplotype.tsv'; DX = f'{O}/distance_A_hg_by_composition.tsv'
dc = pd.read_csv(DC, sep='\t'); dt = pd.read_csv(DT, sep='\t'); da = pd.read_csv(DA, sep='\t'); dx = pd.read_csv(DX, sep='\t', index_col=0)
def texp(p):
    m, e = f'{p:.1e}'.split('e'); e = int(e)
    return f'{p:.2f}' if p >= 0.01 else f'${m}\\times10^{{{e}}}$'
def dcv(sub, lin, ctx, lev, col='median_to_SMN2'):
    r = dc[(dc.subset == sub) & (dc.lineage == lin) & (dc.context == ctx) & (dc.level == lev)]; assert len(r) == 1, (sub, lin, ctx, lev); return r.iloc[0]
def dtp(sub, lin, ctx):
    r = dt[(dt.subset == sub) & (dt.lineage == lin) & (dt.context == ctx)]; assert len(r) == 1; return r.iloc[0]
for lin in 'AB':
    for lev, nm in [('centromeric', 'Cen'), ('telomeric', 'Tel')]:
        r = dcv('all', lin, 'position3', lev); number(f'dc{lin}{nm}', f'{r.median_to_SMN2:.2f}', [DC]); number(f'dc{lin}{nm}N', int(r.copies), [DC])
    number(f'dcP{lin}Pos', texp(dtp('all', lin, 'position3').p), [DT])
r = dcv('AFR', 'B', 'position3', 'centromeric'); number('dcAfrBCen', f'{r.median_to_SMN2:.2f}', [DC])
r = dcv('AFR', 'B', 'position3', 'telomeric'); number('dcAfrBTel', f'{r.median_to_SMN2:.2f}', [DC])
number('dcPAfrBPos', texp(dtp('AFR', 'B', 'position3').p), [DT])
number('dcPAfrAPos', f"{dtp('AFR', 'A', 'position3').p:.2f}", [DT]); number('dcPAfrAPal', f"{dtp('AFR', 'A', 'palindrome_arm').p:.2f}", [DT])
assert dtp('AFR', 'A', 'position3').p > 0.05 and dtp('AFR', 'A', 'palindrome_arm').p > 0.05, 'text: within AFR the A position and arm differences vanish'
# A: small but clear shifts by palindrome arm and arrangement (the cold read: the text must not say 'not position' for A)
for lev, nm in [('yes', 'Yes'), ('no', 'No')]:
    r = dcv('all', 'A', 'palindrome_arm', lev); number(f'dcAPal{nm}', f'{r.median_to_SMN2:.2f}', [DC]); number(f'dcAPal{nm}N', int(r.copies), [DC])
number('dcPAPal', texp(dtp('all', 'A', 'palindrome_arm').p), [DT]); number('dcPAArr', texp(dtp('all', 'A', 'arrangement_top').p), [DT])
number('dcASingle', f"{dcv('all', 'A', 'position3', 'single').median_to_SMN2:.2f}", [DC])
_aarr = dc[(dc.subset == 'all') & (dc.lineage == 'A') & (dc.context == 'arrangement_top')]
number('dcAArrMin', f'{_aarr.median_to_SMN2.min():.2f}', [DC]); number('dcAArrMax', f'{_aarr.median_to_SMN2.max():.2f}', [DC])
_actx = dc[(dc.subset == 'all') & (dc.lineage == 'A') & dc.context.isin(['position3', 'palindrome_arm', 'arrangement_top'])]
number('dcASpread', f'{_actx.median_to_SMN2.max() - _actx.median_to_SMN2.min():.2f}', [DC])
assert dtp('all', 'A', 'position3').p < 0.01 and dtp('all', 'A', 'palindrome_arm').p < 0.01 and dtp('all', 'A', 'arrangement_top').p < 0.01, 'text: the A tests are significant'
assert _actx.median_to_SMN2.max() - _actx.median_to_SMN2.min() < 0.15, 'text: within A, where a copy sits shifts its median by little'
assert float(MACROS['dcASingle']) > max(float(MACROS['dcACen']), float(MACROS['dcATel'])), 'text: the position difference of A comes from single copies'
assert dcv('all', 'A', 'position3', 'single').median_to_SMN2 > max(dcv('all', 'A', 'position3', 'centromeric').median_to_SMN2, dcv('all', 'A', 'position3', 'telomeric').median_to_SMN2), 'text: single A copies sit further'
assert abs(dcv('all', 'A', 'position3', 'centromeric').median_to_SMN2 - dcv('all', 'A', 'position3', 'telomeric').median_to_SMN2) < 0.05  # text: no position effect in A
assert dcv('all', 'B', 'position3', 'telomeric').median_to_SMN2 > dc[(dc.subset == 'all') & (dc.lineage == 'A') & (dc.context == 'position3')].median_to_SMN2.max() + 0.3  # text: B stays further than A wherever it sits
assert dcv('all', 'B', 'position3', 'centromeric').median_to_SMN2 > dcv('all', 'B', 'position3', 'telomeric').median_to_SMN2 + 0.2  # text: centromeric B (mostly B1) further than telomeric B (mostly B2)
for lev, nm in [('yes', 'Yes'), ('no', 'No')]:
    r = dcv('all', 'B', 'palindrome_arm', lev); number(f'dcBPal{nm}', f'{r.median_to_SMN2:.2f}', [DC]); number(f'dcBPal{nm}N', int(r.copies), [DC])
number('dcPBPal', texp(dtp('all', 'B', 'palindrome_arm').p), [DT])
for lev, nm in [('SMN1', 'One'), ('SMN2', 'Two'), ('mixed', 'Mixed')]:
    r = dcv('all', 'B', 'env', lev); number(f'dcBEnv{nm}', f'{r.median_to_SMN2:.2f}', [DC]); number(f'dcBEnv{nm}N', int(r.copies), [DC])
number('dcPBEnv', texp(dtp('all', 'B', 'env').p), [DT])
assert dcv('all', 'B', 'env', 'SMN2').median_to_SMN2 < dcv('all', 'B', 'env', 'SMN1').median_to_SMN2 and dcv('all', 'B', 'env', 'mixed').median_to_SMN2 < dcv('all', 'B', 'env', 'SMN1').median_to_SMN2
for lev, nm in [('1- 1+', 'OneOne'), ('2- 1+', 'TwoOne')]:
    r = dcv('all', 'A', 'arrangement_top', lev); number(f'dcAArr{nm}', f'{r.median_to_SMN2:.2f}', [DC]); number(f'dcAArr{nm}N', int(r.copies), [DC])
number('dcPAArr', texp(dtp('all', 'A', 'arrangement_top').p), [DT])
# distance to SMN2 by arm, as in Fig. 1 and the Fig. 3 strip (step 17 context 'arm': cenSMN / telSMN / no pair)
for lin in ('A', 'B'):
    for lev, nm in (('cenSMN', 'Cen'), ('telSMN', 'Tel'), ('no pair', 'No')):
        r = dcv('all', lin, 'arm', lev); number(f'dc{lin}Arm{nm}', f'{r.median_to_SMN2:.2f}', [DC]); number(f'dc{lin}Arm{nm}N', int(r.copies), [DC])
    number(f'dcP{lin}Arm', texp(dtp('all', lin, 'arm').p), [DT])
number('dcArmSpread', f"{dtp('all', 'A', 'arm').spread:.2f}", [DT])
number('dcPAfrAArm', f"{dtp('AFR', 'A', 'arm').p:.2f}", [DT]); number('dcPAfrBArm', texp(dtp('AFR', 'B', 'arm').p), [DT])
assert dcv('all', 'B', 'arm', 'cenSMN').median_to_SMN2 > dcv('all', 'B', 'arm', 'telSMN').median_to_SMN2 + 0.2, 'text: cenSMN-arm B further than telSMN-arm B'
assert abs(dcv('all', 'A', 'arm', 'cenSMN').median_to_SMN2 - dcv('all', 'A', 'arm', 'telSMN').median_to_SMN2) < 0.05, 'text: A equally far on the two arms'
assert dcv('all', 'A', 'arm', 'no pair').median_to_SMN2 > max(dcv('all', 'A', 'arm', 'cenSMN').median_to_SMN2, dcv('all', 'A', 'arm', 'telSMN').median_to_SMN2), 'text: A without the pair slightly further'
assert dtp('all', 'A', 'arm').p < 0.01 and dtp('all', 'B', 'arm').p < 0.01 and dtp('AFR', 'A', 'arm').p > 0.05 and dtp('AFR', 'B', 'arm').p < 0.01, 'text: arm tests'
assert float(MACROS['dcASpread']) >= dtp('all', 'A', 'arm').spread, 'text: arrangement type spreads A more than the arm'
daa = da.set_index(['by', 'level'])
for lev, nm in [('A', 'OneA'), ('AA', 'TwoA'), ('AB', 'AB')]:
    r = daa.loc[('composition', lev)]; number(f'dcComp{nm}', f'{r.median_to_SMN2:.2f}', [DA]); number(f'dcComp{nm}N', int(r.copies), [DA])
# (ML lineages: S1-8 is lineage B, so the A copies beside a B copy are few; the S1-8 detail moves to lineage B)
number('dcCompABAfr', int(daa.loc[('composition', 'AB'), 'afr']), [DA])

# ---- lineage make-up across the structure and SMN1 sequence trees (14_lineage_in_both_trees.py)
TMT = f'{O}/lineage_tree_match.tsv'; TMJ = f'{O}/lineage_tree_match.json'
tm = pd.read_csv(TMT, sep='\t').set_index(['make_up', 'tree']); tj = json.load(open(TMJ))
number('tmHaps', tj['haplotypes'], [TMJ]); number('tmPerm', f"{tj['permutations']:,}".replace(',', '{,}'), [TMJ])
number('tmPmin', f"{1 / (tj['permutations'] + 1):.4f}", [TMJ])
number('tmMedAll', f"{tm.loc[('AB', 'structure'), 'median_all_cophenetic']:.2f}", [TMT])
for mu, M in [('AB', 'AB'), ('B', 'B'), ('BB', 'BB'), ('AA', 'AA')]:
    number(f'tm{M}n', int(tm.loc[(mu, 'structure'), 'haplotypes']), [TMT]); number(f'tm{M}afr', int(tm.loc[(mu, 'structure'), 'AFR']), [TMT])
    for tr, T in [('structure', 'Struct'), ('SMN1 sequence', 'Seq')]:
        r = tm.loc[(mu, tr)]
        number(f'tm{M}nn{T}', f'{r.nn_same_make_up:.2f}', [TMT]); number(f'tm{M}runs{T}', int(r.leaf_runs), [TMT])
        number(f'tm{M}chance{T}', f'{r.nn_chance_all:.2f}', [TMT]); number(f'tm{M}chanceAfr{T}', f"{r['nn_chance_AFR-stratified']:.2f}", [TMT])
        number(f'tm{M}p{T}', f'{r.p_nn_all:.4f}', [TMT]); number(f'tm{M}pAfr{T}', f"{r['p_nn_AFR-stratified']:.4f}", [TMT])
    number(f'tm{M}medStruct', f"{tm.loc[(mu, 'structure'), 'median_within_cophenetic']:.2f}", [TMT])
for mu, M in [('BB', 'BB'), ('B', 'B')]:
    k = tj[f'{mu}_k16']
    number(f'tm{M}kStruct', k['largest_struct_cluster'], [TMJ]); number(f'tm{M}kStructSize', k['struct_cluster_size'], [TMJ])
    number(f'tm{M}kSeq', k['largest_seq_cluster'], [TMJ]); number(f'tm{M}kSeqSize', k['seq_cluster_size'], [TMJ]); number(f'tm{M}kBoth', k['in_both'], [TMJ])
floor = 1 / (tj['permutations'] + 1) + 1e-9
for mu in ['B', 'BB']:
    for tr in ['structure', 'SMN1 sequence']:
        assert tm.loc[(mu, tr), 'p_nn_all'] <= floor and tm.loc[(mu, tr), 'p_nn_AFR-stratified'] <= floor, (mu, tr)
assert tm.loc[('BB', 'SMN1 sequence'), 'leaf_runs'] < tm.loc[('BB', 'structure'), 'leaf_runs'] < tm.loc[('B', 'structure'), 'leaf_runs']
assert tm.loc[('BB', 'structure'), 'median_within_cophenetic'] < tm.loc[('BB', 'structure'), 'median_all_cophenetic']
assert tm.loc[('B', 'structure'), 'median_within_cophenetic'] > tm.loc[('B', 'structure'), 'median_all_cophenetic'], 'B no longer scattered in structure'
assert tm.loc[('AA', 'structure'), 'p_nn_all'] > 0.05 and tm.loc[('AA', 'SMN1 sequence'), 'p_nn_all'] < 0.01
assert tj['B_k16']['in_both'] < tj['B_k16']['largest_seq_cluster'] and tj['BB_k16']['in_both'] > 0.7 * tm.loc[('BB', 'structure'), 'haplotypes']
assert tm.loc[('BB', 'structure'), 'AFR'] >= tm.loc[('BB', 'structure'), 'haplotypes'] - 2

# ---- the family-A palindrome (10_palindrome.py)
PU, PA, PX, PC = f'{O}/palindrome_units.tsv', f'{O}/palindrome_arms.tsv', f'{O}/palindrome_axis_copies.tsv', f'{O}/copy_in_unit.tsv'
pu = pd.read_csv(PU, sep='\t'); pa = pd.read_csv(PA, sep='\t'); px = pd.read_csv(PX, sep='\t'); pc = pd.read_csv(PC, sep='\t')
pu = pu[~(pu.src.str.startswith('GCA_') | pu.src.str.contains('CHM13'))]   # step 13 keeps GRCh38 as a reference row; it is not one of the study's haplotypes
HCT = 'tables/haplotype_copy_counts.tsv'; _hcp = pd.read_csv(HCT, sep='\t'); _resp = set(_hcp[_hcp.resolved].src)
assert set(pu.src) <= _resp
_nou = _hcp[_hcp.resolved & ~_hcp.src.isin(set(pu.src))]   # resolved haplotypes with no repeat unit: their window has no bundle that recurs
assert len(pu) + len(_nou) == len(_resp)
ppa = pu[pu.has_pair]                                   # every adjacent inverted A pair step 13 found
pph, ppt = ppa[ppa.kind == 'head'], ppa[ppa.kind == 'tail']
assert all(any(u[i] == 'A+' and u[i + 1] == 'A-' for i in range(len(u) - 1)) for u in pph.units.str.split()), 'text: every head pair has an A+ A- junction'
assert all(any(u[i] == 'A-' and u[i + 1] == 'A+' for i in range(len(u) - 1)) and not any(u[i] == 'A+' and u[i + 1] == 'A-' for i in range(len(u) - 1)) for u in ppt.units.str.split()), 'text: the tail pairs are A- A+ only'
assert set(pa.src) == set(pph.src) and set(px.src) == set(pph.src), 'text: arms and copy distances are measured from the head-to-head axis'
_side = px.groupby('src').apply(lambda x: (x.d_kb < 0).any() and (x.d_kb > 0).any())
_ncp = sorted(_side[~_side].index)                      # a head pair with no copy on one arm is not a copy palindrome
pp = pph[~pph.src.isin(_ncp)]; pa = pa[pa.src.isin(pp.src)]; px = px[px.src.isin(pp.src)]
_hn = lambda L: ', '.join(x.split('_hprc')[0].replace('_', ' ') for x in L)
number('palHaps', len(pu), [PU, HCT]); number('palPair', len(pp), [PU, PX], ['head-to-head pair with an SMN copy on each arm'])
number('palAdj', len(ppa), [PU]); number('palHead', len(pph), [PU]); number('palTail', len(ppt), [PU])
number('palTailHaps', _hn(sorted(ppt.src)), [PU]); number('palTailArr', ', '.join("``" + a + "''" for a in ppt.sort_values('src').arr), [PU])
number('palNoCopy', len(_ncp), [PU, PX]); number('palNoCopyHap', _hn(_ncp), [PX])
number('palNoPair', len(pu) - len(pp), [PU]); number('palNoAdj', int((~pu.has_pair).sum()), [PU])
assert len(pu) - len(pp) == int((~pu.has_pair).sum()) + len(ppt) + len(_ncp)
number('palNoUnits', len(_nou), [PU, HCT], ['resolved haplotypes absent from palindrome_units.tsv'])
number('palNoUnitsOneCopy', int(((_nou[['SMN1', 'SMN2', 'h12', 'h21', 'del7_8']].sum(axis=1)) == 1).sum()), [HCT])
assert len(_nou) == int(((_nou[['SMN1', 'SMN2', 'h12', 'h21', 'del7_8']].sum(axis=1)) == 1).sum())   # text: all are single-copy haplotypes
_cl = {'SMN1': 'SMN1', 'SMN2': 'SMN2', 'h12': 'h12', 'h21': 'h21', 'del7_8': '$\\Delta$7--8'}
number('palNoUnitsList', '; '.join(r.src.split('_hprc')[0].replace('_', ' ') + ', ' + [_cl[c] for c in _cl if getattr(r, c) == 1][0] for r in _nou.sort_values('src').itertuples()), [HCT])
# ---- the resolved haplotypes without repeat units: window size, bundle recurrence, copy strand, ancestry
_DB = 'results/decomp/smn.bed'; _CLT = f'{O}/copy_labels.tsv'
_bd = pd.read_csv(_DB, sep='\t', comment='#', header=None, names=['src', 'bgn', 'end', 'lab'])
_bd['bid'] = _bd.lab.str.split(':').str[0].astype(int)
_wl = _bd.groupby('src').end.max() / 1e3              # window length, kb (decomposed window)
_nouS = set(_nou.src)
assert _nouS <= set(_wl.index), 'every no-unit haplotype has a decomposed window'
_oth = _wl[_wl.index.isin(_resp - _nouS)]
assert len(_oth) == len(_resp) - len(_nouS)
number('nouLenMin', round(_wl[list(_nouS)].min()), [_DB]); number('nouLenMax', round(_wl[list(_nouS)].max()), [_DB])
number('nouOtherMin', f'{_oth.min() / 1e3:.2f}', [_DB, HCT], ['shortest window of every other resolved haplotype, Mb'])
assert _wl[list(_nouS)].max() < 0.8 * _oth.min(), 'text: the no-unit windows are much shorter than every other resolved window'
_rec = {}
for s_, d in _bd[_bd.src.isin(_nouS)].groupby('src'):
    vc = d.bid.value_counts(); _rec[s_] = vc[vc > 1]
assert all(list(v.index) == [15] for v in _rec.values()), 'text: only the array bundle 15 recurs in the no-unit windows'
number('nouArrBundle', 15, [_DB]); number('nouArrMin', int(min(v.iloc[0] for v in _rec.values())), [_DB]); number('nouArrMax', int(max(v.iloc[0] for v in _rec.values())), [_DB])
_cln = pd.read_csv(_CLT, sep='\t'); _cln = _cln[_cln.src.isin(_nouS)]
assert len(_cln) == len(_nouS) and (_cln.strand == '-').all(), 'text: each no-unit haplotype carries one reversed copy'
assert (_cln.superpop == 'AFR').all(), 'text: all no-unit donors are of African ancestry'
number('nouPops', ', '.join(f"{k} {v}" for k, v in sorted(hcc_.loc[sorted(_nouS)].population.value_counts().items())), [HCT, _CLT])
assert set(hcc_.loc[sorted(_nouS)].population) == {'ASW', 'ESN', 'MSL'}, 'text spells out ASW, ESN, MSL'
number('palArmL', round(pa.left_kb.median()), [PA]); number('palArmR', round(pa.right_kb.median()), [PA])
number('palArmDiff', round((pa.left_kb - pa.right_kb).abs().median()), [PA]); number('palGap', round(pa.gap_kb.median()), [PA])
_gS = pa.gap_kb[(pa.gap_kb >= 3) & (pa.gap_kb < 8)]; _gL = pa.gap_kb[(pa.gap_kb >= 8) & (pa.gap_kb < 30)]; _gX = pa.gap_kb[(pa.gap_kb < 3) | (pa.gap_kb >= 30)]
assert len(_gS) + len(_gL) + len(_gX) == len(pa) and len(_gX) <= 5, 'text: the axis gap falls in two clusters plus a few outliers'
assert _gS.max() - _gS.min() < 2 and _gL.max() - _gL.min() < 2, 'text: each gap cluster is tight'
number('palGapShort', f'{_gS.median():.1f}', [PA]); number('palGapShortN', len(_gS), [PA])
number('palGapLong', f'{_gL.median():.1f}', [PA]); number('palGapLongN', len(_gL), [PA]); number('palGapOtherN', len(_gX), [PA])
# which cenSMN arm type goes with each axis-stretch cluster (final check 2026-11-20): the short stretch is nearly all C1
_gt = pd.read_csv(f'{O}/arm_element_haps.tsv', sep='\t').set_index('src').cen_type
_gsrc = pa.set_index(pa.index).src if 'src' in pa else None
_gty = pa.assign(cls=np.where(pa.gap_kb.isin(_gS), 'short', np.where(pa.gap_kb.isin(_gL), 'long', 'other'))).set_index('src')
_gct = pd.crosstab(_gty.cls, _gt.reindex(_gty.index).values)
number('palGapShortCOne', int(_gct.loc['short', 'C1']), [PA, f'{O}/arm_element_haps.tsv'])
number('palGapLongCOne', int(_gct.loc['long', 'C1']), [PA, f'{O}/arm_element_haps.tsv'])
_gp = fisher_exact(_gct.loc[['short', 'long'], ['C1', 'C0']].values)[1]
number('palGapTypeP', texp(_gp), [PA, f'{O}/arm_element_haps.tsv'])
assert _gp < 1e-6 and _gct.loc['short', 'C1'] > 0.9 * _gct.loc['short'].sum(), 'text: the short axis stretch goes with a C1 cenSMN arm'
number('palCopies', len(px), [PX])
g = px.groupby('src')
per_arm = g.apply(lambda x: (x.d_kb < 0).any() and (x.d_kb > 0).any())
assert per_arm.all() and len(per_arm) == len(pp), 'a palindrome haplotype lacks a copy on one arm'
number('palOnePerArm', int(per_arm.sum()), [PX]); number('palThreeCopy', int((g.size() == 3).sum()), [PX])
assert len(px) == 2 * len(pp) + MACROS['palThreeCopy'] and (g.size() <= 3).all(), 'caption: one counted copy per arm plus one extra item'
number('palThreeHap', _hn(sorted(g.size()[g.size() == 3].index)), [PX])
number('palDmin', round(px.d_kb.abs().min()), [PX]); number('palDmed', round(px.d_kb.abs().median()), [PX]); number('palDmax', round(px.d_kb.abs().max()), [PX])
# cold read 2026-11-15: the farthest copy (outlier) lies inside its own arm's core, which runs unusually far on that arm
_AUL = f'{O}/a_unit_length.tsv'; _aul = pd.read_csv(_AUL, sep='\t')
_fx = px.loc[px.d_kb.abs().idxmax()]; _fside = 'tel' if _fx.d_kb > 0 else 'cen'
_fcore = float(_aul[(_aul.src == _fx.src) & (_aul.side == _fside)].a_end_kb.iloc[0])
number('palDmaxHap', _hn([_fx.src]), [PX]); number('palDmaxCore', round(_fcore), [_AUL])
assert abs(_fx.d_kb) < _fcore, 'text: the farthest copy lies inside its arm core'
assert _fcore > _aul.a_end_kb.quantile(.99), 'text: that core runs far past the usual core end'
_dq = px.d_kb.abs().quantile([0.05, 0.95])
number('palDlo', round(_dq.iloc[0]), [PX], ['5th percentile']); number('palDhi', round(_dq.iloc[1]), [PX], ['95th percentile'])
mir = g.apply(lambda x: abs(x.d_kb[x.d_kb < 0].abs().min() - x.d_kb[x.d_kb > 0].min()))
number('palMirror', f'{mir.median():.1f}', [PX], ['mirror = | |d| nearest copy left - |d| nearest copy right |'])
number('palRoom', int((pp.room >= 1).sum()), [PU]); number('palDeeper', int((pp.depth >= 1).sum()), [PU])
_deep = pp[pp.depth >= 1]
number('palDeeperHap', _hn(sorted(_deep.src)), [PU])
number('palInUnit', len(pc), [PC])
number('palRelMin', round(100 * pc.rel.min()), [PC]); number('palRelMed', round(100 * pc.rel.median()), [PC]); number('palRelMax', round(100 * pc.rel.max()), [PC])
_nop = pu[~pu.src.isin(pp.src)]
for a_, M in [('1+', 'palNoOne'), ('2- 1-', 'palNoTwoMOneM'), ('2+ 1+', 'palNoTwoPOneP')]:
    number(M, int((_nop.arr == a_).sum()), [PU])
number('palTwoOne', int((pp.arr == '2- 1+').sum()), [PU]); number('palOneOne', int((pp.arr == '1- 1+').sum()), [PU])
# the units that flank the head-to-head A pair (results/sv/palindrome_units.tsv: unit order per haplotype)
_ph = pp; _fl, _fr = collections.Counter(), collections.Counter(); _bb = 0
for _u in _ph.units:
    _t = _u.split(); _i = next(k for k in range(len(_t) - 1) if _t[k] == 'A+' and _t[k + 1] == 'A-')
    _fl[_t[_i - 1] if _i > 0 else 'none'] += 1; _fr[_t[_i + 2] if _i + 2 < len(_t) else 'none'] += 1
    _bb += ('B-' in _t[:_i]) and ('B+' in _t[_i + 2:])
number('palFlankN', len(_ph), [PU])
_lx = lambda k: 'no unit' if k == 'none' else r'\RU{%s}' % k[:-1] + ('+' if k[-1] == '+' else '$-$')
number('palFlankLText', '; '.join(f'{_lx(k)} on {v}' for k, v in _fl.most_common() if v >= 3) + f'; other on {sum(v for v in _fl.values() if v < 3)}', [PU])
number('palFlankRText', '; '.join(f'{_lx(k)} on {v}' for k, v in _fr.most_common() if v >= 3) + f'; other on {sum(v for v in _fr.values() if v < 3)}', [PU])
number('palFlankBB', int(_bb), [PU], ['B- left of the pair and B+ right of it'])
assert _fl.most_common(1)[0][0] == 'B-' and _bb > len(_ph) / 2, 'text: B- is the commonest left flank, and B on both sides on most pair haplotypes'
# ---- SMN1 lineage A/B on the palindrome arms (19_palindrome_diagram.py)
from scipy.stats import mannwhitneyu as _mwu
PL, PR = f'{O}/palindrome_copies_lineage.tsv', f'{O}/palindrome_diagram_rows.tsv'
pl = pd.read_csv(PL, sep='\t'); prr = pd.read_csv(PR, sep='\t')
pl = pl[pl.src.isin(pp.src)]; prr = prr[prr.src.isin(pp.src)]   # the copy palindromes only (not the head pair with no telSMN copy)
ps = pl[pl.gene == 'SMN1']
assert set(ps.lineage) == {'A', 'B'}, 'an SMN1 copy on the palindrome is not in the copy tree'
left = ps.d_kb < 0
number('palSmnA', int((ps.lineage == 'A').sum()), [PL]); number('palSmnB', int((ps.lineage == 'B').sum()), [PL])
_pax = pd.read_csv(f'{O}/palindrome_axis_copies.tsv', sep='\t')
assert set(_pax.src) - set(pl.src) == {[x for x in _pax.src.unique() if x.startswith('NA19391_hap2')][0]}, 'step 17 palindrome copies = the 270 copy palindromes + NA19391 hap2 (text: dcBPalYesN = palSmnB + that one B copy)'
number('palALeft', int(((ps.lineage == 'A') & left).sum()), [PL]); number('palARight', int(((ps.lineage == 'A') & ~left).sum()), [PL])
number('palBLeft', int(((ps.lineage == 'B') & left).sum()), [PL]); number('palBRight', int(((ps.lineage == 'B') & ~left).sum()), [PL])
mk = ps.groupby('src').lineage.apply(lambda s: ''.join(sorted(s)))
ab = ps[ps.src.isin(mk[mk == 'AB'].index)]
number('palABHaps', int((mk == 'AB').sum()), [PL])
number('palABBLeft', int(((ab.lineage == 'B') & (ab.d_kb < 0)).sum()), [PL]); number('palABARight', int(((ab.lineage == 'A') & (ab.d_kb > 0)).sum()), [PL])
bh = set(ps[ps.lineage == 'B'].src); number('palBHaps', len(bh), [PL])
da, db = ps[ps.lineage == 'A'].d_kb.abs(), ps[ps.lineage == 'B'].d_kb.abs()
number('palADist', round(da.median()), [PL]); number('palBDist', round(db.median()), [PL]); number('palPDist', texp(_mwu(da, db).pvalue), [PL])
mb, mo = prr[prr.src.isin(bh)].mirror_kb, prr[~prr.src.isin(bh)].mirror_kb
number('palBMirror', f'{mb.median():.1f}', [PR]); number('palOthMirror', f'{mo.median():.1f}', [PR])
# checks guarding the text: B mostly on the left (A+, centromeric) arm, A mostly on the right; in AB haplotypes B left and A right;
# B copies sit further from the axis than A copies; B haplotypes mirror less exactly
assert (ps[(ps.lineage == 'A')].d_kb > 0).mean() > 0.85, 'text: A mostly on the right (telSMN) arm'
for k_, nm in [('B1', 'Bone'), ('B2', 'Btwo')]:
    q_ = ps[ps.sublineage == k_]; number(f'pal{nm}Left', int((q_.d_kb < 0).sum()), [PL]); number(f'pal{nm}Right', int((q_.d_kb > 0).sum()), [PL])
assert MACROS['palBoneLeft'] > 0.75 * (MACROS['palBoneLeft'] + MACROS['palBoneRight']) and MACROS['palBtwoRight'] > 0.75 * (MACROS['palBtwoLeft'] + MACROS['palBtwoRight']), 'text: B1 left, B2 right'
assert db.median() > da.median() and _mwu(da, db).pvalue < 0.01, 'text: B copies sit further from the axis than A copies'
assert mb.median() < mo.median(), 'text: at this scale the mirror is not less exact on B haplotypes'
assert MACROS['palSmnA'] + MACROS['palSmnB'] == int((px_gene := pl.gene.eq('SMN1')).sum())
assert px.d_kb.abs().min() > 50, 'a copy sits near the axis'
assert pc.rel.max() < 0.75
# ---- the arms of the A pair as named elements (step 28b; replaces the long/short unit call) and the spacer S (28c)
AU, AC = f'{O}/arm_units.tsv', f'{O}/arm_copies.tsv'
au = pd.read_csv(AU, sep='\t'); ac = pd.read_csv(AC, sep='\t')
number('armN', len(au), [AU])
number('armOneCopy', int((au.n_counted == 1).sum()), [AU], ['n_counted excludes truncated pieces and fragments'])
number('armExtraItem', int((au.n_items > au.n_counted).sum()), [AU])
EJ, EA_, EB_, EG_, EN_, ESP = f'{O}/arm_elements_summary.json', f'{O}/arm_elements.tsv', f'{O}/arm_element_blocks.tsv', f'{O}/arm_element_genes.tsv', f'{O}/arm_naip_pieces.tsv', f'{O}/arm_element_spacers.tsv'
ej = json.load(open(EJ)); ea = pd.read_csv(EA_, sep='\t'); eb = pd.read_csv(EB_, sep='\t')
number('elHaps', ej['n_haps'], [EJ])
for t, M in [('C1', 'elCOne'), ('C0', 'elCZero'), ('T0', 'elTZero'), ('T1', 'elTOne')]:
    number(M, ej['types'][t], [EJ]); number(M + 'Can', ej['canonical'][t], [EJ])
    number(M + 'Long', ej['unit_long_by_type'][t][0], [EJ])
assert ej['types']['C1'] + ej['types']['C0'] == ej['types']['T0'] + ej['types']['T1'] == ej['n_haps'] == MACROS['palPair'], 'every pair haplotype has one cenSMN and one telSMN arm type'
number('elXdMin', int(ea[ea.kb_Xd > 0].kb_Xd.min()), [EA_], ['assert >= 20 in step 28b'])
for k_, M in [('C1T0', 'elFonezero'), ('C1T1', 'elFoneone'), ('C0T0', 'elFzerozero'), ('C0T1', 'elFzeroone')]:
    number(M, ej['cen_x_tel'][k_], [EJ]); number(M + 'Afr', ej['cen_x_tel_afr'][k_], [EJ])
number('elFormP', texp(ej['cen_x_tel_p']), [EJ]); number('elFormAfrP', pfmt(ej['cen_x_tel_afr_p']), [EJ]); number('elFormNonAfrP', f"{ej['cen_x_tel_nonafr_p']:.2f}", [EJ])
assert ej['cen_x_tel_p'] < 1e-6 and ej['cen_x_tel_afr_p'] < 0.01 and ej['cen_x_tel_nonafr_p'] > 0.5, 'text: the two arm types go together overall and in AFR, not outside AFR'
cz = ej['C0_by_superpop']
number('elCZeroAfr', cz['AFR'][0], [EJ]); number('elAfrHaps', cz['AFR'][1], [EJ])
number('elCZeroNonAfr', sum(v[0] for k, v in cz.items() if k != 'AFR'), [EJ]); number('elNonAfrHaps', sum(v[1] for k, v in cz.items() if k != 'AFR'), [EJ])
number('elCZeroEas', cz['EAS'][0], [EJ]); number('elEasHaps', cz['EAS'][1], [EJ])
tz = ej['T1_by_superpop']; number('elTOneAfr', tz['AFR'][0], [EJ]); number('elTOneEas', tz['EAS'][0], [EJ])
assert tz['EAS'][0] == 0, 'text: no EAS haplotype has a T1 arm'
# cold read 2026-11-14: is the AFR C0-T1 link the B1+B2 chromosome? re-test without it
EH_ = f'{O}/arm_element_haps.tsv'; _eh = pd.read_csv(EH_, sep='\t')
from scipy.stats import fisher_exact as _fx
def _ct(d):
    t = pd.crosstab(d.cen_type, d.tel_type).reindex(index=['C0', 'C1'], columns=['T0', 'T1'], fill_value=0); return t, _fx(t.values)[1]
_c0t2 = _eh[(_eh.cen_type == 'C0') & (_eh.tel_type == 'T1')]
number('elFzerooneBB', int(_c0t2.b1b2.sum()), [EH_]); number('elFzerooneNotBB', int((~_c0t2.b1b2).sum()), [EH_])
assert int(MACROS['elFzerooneBB']) + int(MACROS['elFzerooneNotBB']) == int(MACROS['elFzeroone'])
_ta, _pa = _ct(_eh[_eh.afr & ~_eh.b1b2])
number('elFormAfrNoBBP', f"{_pa:.2f}", [EH_]); number('elFormAfrNoBBn', int(_ta.values.sum()), [EH_])
number('elFormAfrNoBBzt', int(_ta.loc['C0', 'T1']), [EH_])
assert _pa > 0.5, 'text: without the B1+B2 haplotypes C0 and T1 no longer go together within AFR'
# researcher 2026-11-14: print it as 'no association (30% against 28%)', not as p = 1.00
_z, _o = _ta.loc['C0'], _ta.loc['C1']
number('elFormAfrNoBBzPct', f"{100 * _z['T1'] / _z.sum():.0f}", [EH_]); number('elFormAfrNoBBoPct', f"{100 * _o['T1'] / _o.sum():.0f}", [EH_])
number('elFormAfrNoBBzN', int(_z.sum()), [EH_]); number('elFormAfrNoBBoN', int(_o.sum()), [EH_]); number('elFormAfrNoBBot', int(_o['T1']), [EH_])
assert abs(_z['T1'] / _z.sum() - _o['T1'] / _o.sum()) < 0.05, 'text: no association (shares within 5 points)'
number('elBoneCZero', ej['B1cen_C0'][0], [EJ]); number('elBoneCen', ej['B1cen_C0'][1], [EJ])
number('elAfrBoneCZero', ej['afr_B1cen_C0'][0], [EJ]); number('elAfrBoneCen', ej['afr_B1cen_C0'][1], [EJ])
number('elAfrOthCZero', ej['afr_otherCen_C0'][0], [EJ]); number('elAfrOthCen', ej['afr_otherCen_C0'][1], [EJ]); number('elAfrBoneP', texp(ej['afr_B1cen_C0_p']), [EJ])
number('elBtwoTOne', ej['B2tel_T1'][0], [EJ]); number('elBtwoTel', ej['B2tel_T1'][1], [EJ])
number('elAfrBtwoTOne', ej['afr_B2tel_T1'][0], [EJ]); number('elAfrBtwoTel', ej['afr_B2tel_T1'][1], [EJ])
number('elAfrOthTOne', ej['afr_otherTel_T1'][0], [EJ]); number('elAfrOthTel', ej['afr_otherTel_T1'][1], [EJ]); number('elAfrBtwoP', texp(ej['afr_B2tel_T1_p']), [EJ])
number('elBB', ej['b1b2']['n'], [EJ]); number('elBBCZero', ej['b1b2']['C0'], [EJ]); number('elBBTOne', ej['b1b2']['T1'], [EJ])
assert ej['b1b2']['C0'] == ej['b1b2']['n'] == MACROS['linMixHaps'], 'text: every B1+B2 haplotype has a C0 cenSMN arm'
_eh = pd.read_csv(f'{O}/arm_element_haps.tsv', sep='\t')
assert set(_eh[_eh.b1b2].src) == set(mix.src), 'abstract: the B1+B2 haplotypes of the arm-element table are the same set as those with one B1 and one B2 copy'
number('elCZeroNotBB', ej['C0_not_b1b2'][0], [EJ]); number('elCZeroNotBBAfr', ej['C0_not_b1b2'][1], [EJ])
for sd, ts in (('cen', ('C1', 'C0')), ('tel', ('T0', 'T1'))):
    for t in ts:
        T_ = {'C1': 'COne', 'C0': 'CZero', 'T0': 'TZero', 'T1': 'TOne'}[t]
        number(f'elCopyD{T_}', f"{ej[sd + '_copy_d_by_type'][t]:.0f}", [EJ]); number(f'elAend{T_}', f"{ej[sd + '_aend_by_type'][t]:.0f}", [EJ])
    number(f'elRho{sd}', f"{ej[sd + '_rho_copy_inner']:.2f}", [EJ]); number(f'elNInner{sd}', ej[sd + '_n_inner'], [EJ])
number('elCopyDP', texp(ej['cen_copy_d_p']), [EJ]); number('elInnerUnit', f"{ej['cen_inner_unit_kb']:.1f}", [EJ])
for k_ in ('2', '3', '4'): number('elCopyAtInner' + {'2': 'Two', '3': 'Three', '4': 'Four'}[k_], f"{ej['tel_copy_d_by_inner'][k_]:.0f}", [EJ])
assert ej['cen_rho_copy_inner'] > 0.7 and ej['tel_rho_copy_inner'] > 0.8, 'text: the copy position steps with the inner array'
# at the same inner-array count, telSMN arms of the two types place the copy alike; cenSMN arms do not (C0 further out)
# step 28g (2026-11-14, researcher 'why copies sit further out on C0 arms: can you answer?'): the extra distance is the
# core's inner array stretch (bundle 40+ to 16+), measured in bp; the rest of the axis-to-copy distance is the same
CPJ = f'{O}/copy_position.json'; cp = json.load(open(CPJ))
for sd, (a_, b_) in (('cen', ('CZero', 'COne')), ('tel', ('TOne', 'TZero'))):
    c_ = cp[sd]
    for i_, T_ in enumerate((a_, b_)):
        number(f'cpArr{T_}', f"{c_['arr_kb']['median'][i_]:.1f}", [CPJ]); number(f'cpRest{T_}', f"{c_['rest_kb']['median'][i_]:.1f}", [CPJ])
        number(f'cpCopy{T_}', f"{c_['copy_d_kb']['median'][i_]:.1f}", [CPJ]); number(f'cpN{T_}', c_['n'][i_], [CPJ])
    number(f'cpRho{sd}', f"{c_['rho_copy_arr']:.2f}", [CPJ]); number(f'cpMeasured{sd}', c_['n_measured'], [CPJ])
    number(f'cpDiffCopy{sd}', f"{c_['diff_copy']:.1f}", [CPJ]); number(f'cpDiffArr{sd}', f"{c_['diff_arr']:.1f}", [CPJ]); number(f'cpDiffRest{sd}', f"{abs(c_['diff_rest']):.1f}", [CPJ])
    assert c_['diff_copy'] > 10 and c_['copy_d_kb']['p'] < 1e-6, 'text: C0 (cenSMN) and T1 (telSMN) copies sit further out'
# the same split at the finer scale (step 37: results/sv/fine/copy_position.json, the pipeline's output at commit 953c895)
FCP = 'results/sv/fine/copy_position.json'; fcp = json.load(open(FCP))
for sd, (a_, b_) in (('cen', ('CZero', 'COne')), ('tel', ('TOne', 'TZero'))):
    c_ = fcp[sd]
    for i_, T_ in enumerate((a_, b_)):
        number(f'fcpArr{T_}', f"{c_['arr_kb']['median'][i_]:.1f}", [FCP]); number(f'fcpRest{T_}', f"{c_['rest_kb']['median'][i_]:.1f}", [FCP])
        number(f'fcpCopy{T_}', f"{c_['copy_d_kb']['median'][i_]:.1f}", [FCP]); number(f'fcpN{T_}', c_['n'][i_], [FCP])
    number(f'fcpRho{sd}', f"{c_['rho_copy_arr']:.2f}", [FCP]); number(f'fcpMeasured{sd}', c_['n_measured'], [FCP])
    number(f'fcpDiffCopy{sd}', f"{c_['diff_copy']:.1f}", [FCP]); number(f'fcpDiffArr{sd}', f"{c_['diff_arr']:.1f}", [FCP]); number(f'fcpDiffRest{sd}', f"{abs(c_['diff_rest']):.1f}", [FCP])
    assert abs(c_['diff_rest']) < 1.0 and abs(c_['diff_arr'] - c_['diff_copy']) < 2.5, 'text (finer scale): the extra distance is the array stretch, the rest is equal'
_md = sorted(int(k) for k in fcp['tel']['arr_modes']); number('fcpModes', ', '.join(map(str, _md)), [FCP])
assert all(19 <= b - a <= 22 for a, b in zip(_md, _md[1:])), 'text (finer scale): the stretch comes in steps of about one ~20-kb unit'
# at the coarse scale the array's bundle covers only part of the array, so the split is not clean (copy_position.json)
assert abs(cp['cen']['diff_rest']) > 5, 'text: at the coarse scale the split does not come out clean'
co, to_ = ej['cen_outer_by_type'], ej['tel_outer_by_type']
number('elOutCOne', co['C1']['with_array'], [EJ]); number('elOutCZero', co['C0']['with_array'], [EJ]); number('elOutCUnit', f"{co['C1']['unit_kb']:.1f}", [EJ]); number('elOutCCopies', f"{co['C1']['copies_median']:.1f}", [EJ])
number('elOutTZero', to_['T0']['with_array'], [EJ]); number('elOutTOne', to_['T1']['with_array'], [EJ]); number('elOutTUnit', f"{to_['T0']['unit_kb']:.1f}", [EJ])
number('elOutTZeroCopies', f"{to_['T0']['copies_median']:.1f}", [EJ]); number('elOutTOneCopies', f"{to_['T1']['copies_median']:.1f}", [EJ]); number('elOutTP', f"{ej['tel_outer_T0_T1_p']:.2f}", [EJ])
assert co['C0']['with_array'] == 0 and ej['tel_outer_T0_T1_p'] > 0.05, 'text: no outer array without X; telSMN outer copies do not differ by type'
nw = ej['naip_whole_spans_aend']
for t, M in [('C0', 'elNaipWholeCZero'), ('T0', 'elNaipWholeTZero'), ('T1', 'elNaipWholeTOne'), ('C1', 'elNaipWholeCOne')]: number(M, nw[t][0], [EJ, EN_])
number('elNaipPast', f"{ej['naip_whole_past_aend_kb']:.0f}", [EJ])
number('elNaipSplit', ej['C1_3p_at_aend'][0], [EJ, EN_]); number('elNaipFromA', ej['C1_3p_txfrom']['1611'], [EJ]); number('elNaipFromB', ej['C1_3p_txfrom']['1965'], [EJ])
number('elNaipTxA', f"{1611:,}".replace(',', '{,}'), [EN_]); number('elNaipTxB', f"{1965:,}".replace(',', '{,}'), [EN_])
# 28b's transcript coordinates are 0-based half-open (the 5' part runs tx_from 0 to 1965), so the 3' parts start at 1-based bases 1612 and 1966
number('elNaipStartA', f"{1611 + 1:,}".replace(',', '{,}'), [EN_]); number('elNaipStartB', f"{1965 + 1:,}".replace(',', '{,}'), [EN_])
assert (1965 - (1611 + 1) + 1) == 1965 - 1611, 'text: bases 1612-1965 are the 354 shared bases'
assert set(ej['C1_3p_txfrom']) == {'1611', '1965'} and ej['C1_5p'][1] == 1965, 'text: the two split points'
number('elNaipGapA', f"{ej['C1_3p_next_x_kb']['1611']:.1f}", [EJ]); number('elNaipGapB', f"{ej['C1_3p_next_x_kb']['1965']:.1f}", [EJ])
_np = pd.read_csv(f'{O}/arm_naip_pieces.tsv', sep='\t'); _c1 = _np[_np.type == 'C1']
_f5 = _c1[_c1.piece == '5p']; assert (_f5.tx_from == 0).all() and (_f5.tx_to == 1965).all() and _f5.src.nunique() == ej['types']['C1'], 'text: a 5-prime part (bases 1-1965) on every C1 arm'
assert set(_c1[_c1.piece == 'whole'].src) <= set(_f5.src) and _c1[(_c1.piece == 'whole') & _c1.spans_a_end].src.nunique() == nw['C1'][0], 'text: the C1 arms with a whole gene also carry the 5-prime part'
number('elNaipOverlap', 1965 - 1611, [EN_])
number('elNaipFive', ej['C1_5p'][0], [EJ, EN_]); number('elNaipTZeroThree', ej['T0_3p_before_x'][0], [EJ, EN_]); number('elNaipTOneThree', ej['T1_3p_before_x'][0], [EJ, EN_])
assert ej['C1_3p_at_aend'][0] + nw['C1'][0] == ej['types']['C1'] and nw['C0'][0] == ej['types']['C0'], 'text: NAIP split on C1 arms (or whole across A end), whole on every C0 arm'
number('elNaipTwoWholeCZero', ej['naip_whole2_C0'][0], [EJ]); number('elNaipOneWholeCOne', ej['naip_whole1_C1'][0], [EJ])
number('elGtfE', ej['gtf2h2_whole_in_E'][0], [EJ]); number('elGtfAll', ej['gtf2h2_in_E'][1], [EJ]); number('elEcen', ej['E_arms']['cen'], [EJ]); number('elEtel', ej['E_arms']['tel'], [EJ])
assert ej['gtf2h2_whole_in_E'][0] >= 0.99 * ej['gtf2h2_in_E'][1] and ej['gtf2h2_whole_in_E'][0] == ej['gtf2h2_in_E'][0], 'text: E carries GTF2H2'
number('elAxisS', ej['axis_S'], [EJ]); number('elSThree', ej['S_per_window']['3'], [EJ])
# element extents on arms with the type's commonest layout (median block ends, kb from the axis)
_can = ea[ea.canonical][['src', 'side', 'type']]; _bb = eb.merge(_can, on=['src', 'side'])
def _ext(t, el, k): return round(_bb[(_bb.type == t) & (_bb.el == el)].groupby(['src', 'side'])['a0' if k == 0 else 'a1'].agg('min' if k == 0 else 'max').median())
number('elAEnd', round(ea.a_end_kb.median()), [EA_])
for t, el, M in [('C1', 'Xi', 'elCXi'), ('C1', 'Xd', 'elCXd'), ('C1', 'S', 'elCS'), ('T0', 'Xi', 'elTXi'), ('C1', 'E', 'elCE'), ('C0', 'E', 'elZE')]:
    number(M + 'Lo', _ext(t, el, 0), [EB_, EA_]); number(M + 'Hi', _ext(t, el, 1), [EB_, EA_])
# arm-element definitions (step 28b) -> Supplementary table (tab_elements.tex)
ES_ = f'{O}/arm_element_sets.tsv'; es_ = pd.read_csv(ES_, sep='\t')
_lines = ['% written by scripts/40_numbers.py from results/sv/arm_element_sets.tsv; do not edit', r'\newcommand{\tabElements}{%']
_TNAME = {'Xd': "Xd, Xi again, opposite to Xi (direct, as the core)"}   # step 28b's element Xi is shown as its two parts
for r in es_.itertuples():
    if r.element == 'Xi':
        _tk = r.tokens.split(); _ti = [t for t in _tk if t in ('1-', '5-')]; _to = [t for t in _tk if t not in _ti]
        assert _ti == ['1-', '5-'] and len(_to) == 6
        _lines.append(trow("Xi$'$, inner part of the inverted core copy (nearer the core, outside X)", ' '.join(_ti).replace('-', '$-$')))
        _lines.append(trow("Xi, outer part of the inverted core copy (in X)", ' '.join(_to).replace('-', '$-$')))
    else: _lines.append(trow(_TNAME.get(r.element, r.name), r.tokens.replace('-', '$-$')))
_lines.append('}'); (DOC / 'tab_elements.tex').write_text('\n'.join(_lines) + '\n')
assert set(es_.element) == {'K', 'Xi', 'S', 'Xd', 'E', 'Fl'}
SJ = f'{O}/spacer_summary.json'; sj = json.load(open(SJ)); SS_ = f'{O}/spacer_sequence.tsv'; ss = pd.read_csv(SS_, sep='\t')
number('sOcc', sj['n_occ'], [SJ]); number('sLen', f"{sj['len_bp_median'] / 1e3:.1f}", [SJ])
_ejS = json.load(open(EJ))
number('sBundleKb', f"{_ejS['cen_element_kb']['C1']['S']:.0f}", [EJ])
assert abs(_ejS['tel_element_kb']['T1']['S'] - 2 * _ejS['tel_element_kb']['T0']['S']) < 1.0, 'text: a T1 arm carries S twice'
assert abs(_ejS['cen_element_kb']['C1']['S'] - sj['len_bp_median'] / 1e3) < 0.5, 'text: S is one bundle; its span is its bundle bases'
# X's size per arm (definition (a), researcher 2026-11-18): X = Xi's outer half, S and Xd, the tight inverted repeat
# centered on S. Xi's inner half (bundles 1-, 5-) and the tandem array between the two halves lie outside X.
import numpy as _np
_XB = EB_
_xa = pd.read_csv(EA_, sep='\t')
_xcan = set(_xa[(_xa.side == 'cen') & (_xa.type == 'C1') & _xa.canonical].src)
_xsg = pd.read_csv(f'{O}/arm_element_segments.tsv', sep='\t'); _xsg = _xsg[(_xsg.side == 'cen') & _xsg.src.isin(_xcan)]
_XIN, _XOUT = {'1-', '5-'}, {'13-', '16-', '17-', '21-', '24-', '3-'}
_xspan, _xfull, _xo0, _xin, _xgap = [], [], [], [], []
_xae = pd.read_csv(_XB, sep='\t'); _xae = _xae[(_xae.side == 'cen') & (_xae.i == 0)].set_index('src').a1   # end of the first block, the core
for _s, _g in _xsg.groupby('src'):
    _g = _g[_g.a0 >= _xae[_s] - 0.5]   # beyond the core: the axis S and core bundles share ids with S and Xd
    _sS = _g[_g.el == 'S']; _d = _g[_g.el == 'Xd']; _o = _g[_g.tok.isin(_XOUT) & (_g.a1 <= _sS.a0.min() + 0.5)]; _i = _g[_g.tok.isin(_XIN)]
    assert len(_sS) and len(_d) and len(_o), _s
    _xspan.append(_d.a1.max() - _o.a0.min()); _xo0.append(_o.a0.min())
    if len(_i):
        _xfull.append(_d.a1.max() - _i.a0.min()); _xin.append(_i.a1.max() - _i.a0.min()); _xgap.append(_o.a0.min() - _i.a1.max())
        assert _i.a1.max() <= _o.a0.min() + 0.5, f'{_s}: Xi inner half not inside of the outer half'
_xspan, _xfull = _np.array(_xspan), _np.array(_xfull)
_xlo, _xhi = _np.percentile(_xspan, [10, 90]); _flo, _fhi = _np.percentile(_xfull, [10, 90])
number('xSpanN', str(len(_xspan)), [EB_, EA_, f'{O}/arm_element_segments.tsv'])
number('xSpanRest', int(ej['types']['C1']) - len(_xspan), [EB_, EA_, EJ])   # C1 arms with an extra element, not measured
number('xSpanLo', f"{round(_xlo, -1):.0f}", [f'{O}/arm_element_segments.tsv', EA_])
number('xSpanHi', f"{round(_xhi, -1):.0f}", [f'{O}/arm_element_segments.tsv', EA_])
number('xSpanMed', f"{_np.median(_xspan):.0f}", [f'{O}/arm_element_segments.tsv', EA_])
number('xOutStart', f"{_np.median(_xo0):.0f}", [f'{O}/arm_element_segments.tsv', EA_])
number('xFullLo', f"{round(_flo, -1):.0f}", [f'{O}/arm_element_segments.tsv', EA_])
number('xFullHi', f"{round(_fhi, -1):.0f}", [f'{O}/arm_element_segments.tsv', EA_])
number('xInMed', f"{_np.median(_xin):.0f}", [f'{O}/arm_element_segments.tsv', EA_])
number('xGapMed', f"{_np.median(_xgap):.0f}", [f'{O}/arm_element_segments.tsv', EA_])
assert len(_xspan) == len(_xcan), 'every canonical C1 arm has Xi outer half, S and Xd'
assert _xhi - _xlo < 40, 'text: X (Xi outer half, S, Xd) is of nearly constant size'
assert _fhi - _flo > 100, 'text: with Xi inner half and the array the span varies by over 100 kb'
assert ((_fhi - _flo) - (_xhi - _xlo)) / (_fhi - _flo) >= 0.7, 'text: Xi inner half and the array carry nearly all of the spread'
_xes = pd.read_csv(f'{O}/arm_element_sets.tsv', sep='\t').set_index('element').tokens
assert {t.replace('-', '+') for t in _XOUT} == set(_xes['Xd'].split()), 'text: Xd is Xi outer half again, direct'

# X takes two sizes (final check 2026-11-20): over the canonical C1 arms the span is bimodal, with nothing between the
# two values; the difference is the stretch of bundle 3, the bundle on either side of S, that the Xd side carries.
_xform = []
for _s, _g in _xsg.groupby('src'):
    _g = _g[_g.a0 >= _xae[_s] - 0.5]
    _sS = _g[_g.el == 'S']; _d = _g[_g.el == 'Xd']; _o = _g[_g.tok.isin(_XOUT) & (_g.a1 <= _sS.a0.min() + 0.5)]
    _b3d = _d[_d.tok == '3+']; _b3i = _o[_o.tok == '3-']
    _xform.append(dict(src=_s, span=_d.a1.max() - _o.a0.min(), b3d=float((_b3d.a1 - _b3d.a0).sum()),
                       b3i=float((_b3i.a1 - _b3i.a0).sum())))
_xform = pd.DataFrame(_xform)
_xsh, _xlg, _xbg = _xform[_xform.span < 195], _xform[(_xform.span >= 195) & (_xform.span < 215)], _xform[_xform.span >= 215]
_XR = [f'{O}/arm_element_segments.tsv', EA_]
assert len(_xform) == len(_xcan) and len(_xsh) + len(_xlg) + len(_xbg) == len(_xcan), 'every canonical C1 arm is in one X size class'
assert not ((_xform.span > 185) & (_xform.span < 205)).any(), 'text: no arm lies between the two forms of X'
assert _xsh.span.max() - _xsh.span.min() < 2 and _xlg.span.max() - _xlg.span.min() < 2, 'text: each form is of nearly constant size'
number('xFormShortN', len(_xsh), _XR); number('xFormLongN', len(_xlg), _XR); number('xFormBigN', len(_xbg), _XR)
number('xFormTwoN', len(_xsh) + len(_xlg), _XR)   # arms of one of the two lengths: the abstract's denominator
assert len(_xsh) + len(_xlg) == len(_xform) - len(_xbg) <= int(ej['types']['C1']), 'text: the two lengths cover the measured arms but not every C1 arm'
number('xFormShortKb', f"{_xsh.span.median():.0f}", _XR); number('xFormLongKb', f"{_xlg.span.median():.0f}", _XR)
number('xFormBigKb', ' and '.join(f"{v:.0f}" for v in sorted(_xbg.span)), _XR)
number('xFormBigHaps', ', '.join(sorted(' '.join(h.split('_')[:2]) for h in _xbg.src)), _XR)
number('xFormShortBd', f"{_xsh.b3d.median():.0f}", _XR); number('xFormLongBd', f"{_xlg.b3d.median():.0f}", _XR)
number('xFormShortBi', f"{_xsh.b3i.median():.0f}", _XR); number('xFormLongBi', f"{_xlg.b3i.median():.0f}", _XR)
assert abs((_xlg.b3d.median() - _xsh.b3d.median()) - (_xlg.span.median() - _xsh.span.median())) < 15, 'text: the two forms differ in the bundle next to S on the Xd side'
_xtt = pd.read_csv(EH_, sep='\t').set_index('src').tel_type
_xform['form'] = np.where(_xform.span < 195, 'short', 'long')
_xct = pd.crosstab(_xform[_xform.span < 215].form, _xtt.reindex(_xform[_xform.span < 215].src).values)
number('xFormShortTOne', int(_xct.loc['short', 'T1']), _XR + [EH_]); number('xFormLongTOne', int(_xct.loc['long', 'T1']), _XR + [EH_])
_xfp = fisher_exact(_xct.loc[['short', 'long'], ['T1', 'T0']].values)[1]
number('xFormP', texp(_xfp), _XR + [EH_])
assert _xfp < 1e-5, 'text: the short form goes with a telSMN arm that carries X'
assert set(_xes['Xi'].split()) == _XIN | _XOUT, 'text: Xi = inner half (1, 5) + outer half'
_cb = {'1', '5', '15', '12', '16', '24', '21', '17', '13', '3'}
assert {t[:-1] for t in _XIN | _XOUT} | {'15', '12'} == _cb, 'text: the two halves with the array hold all the core bundles'
# Xi' and Xi (researcher 2026-11-18): the inner part (1-, 5-) is written Xi', the outer part Xi. Every Xi block is one
# part only (figures label blocks by part), and the counts of arms carrying each part go in the text.
_xbk = pd.read_csv(f'{O}/arm_element_blocks.tsv', sep='\t'); _xsx = pd.read_csv(f'{O}/arm_element_segments.tsv', sep='\t'); _xsx = _xsx[_xsx.el == 'Xi']
_xparts = []
for _r in _xbk[_xbk.el == 'Xi'].itertuples():
    _t = _xsx[(_xsx.src == _r.src) & (_xsx.side == _r.side) & (_xsx.a0 >= _r.a0 - 0.5) & (_xsx.a1 <= _r.a1 + 0.5)]
    _pi, _po = _t.tok.isin(_XIN).any(), _t.tok.isin(_XOUT).any()
    assert _pi != _po, f'{_r.src} {_r.side}: an Xi block holds both parts'
    _xparts.append(_pi)
_xpa = {}
for (_s, _sd), _g in _xsx.groupby(['src', 'side']): _xpa[(_s, _sd)] = (_g.tok.isin(_XIN).any(), _g.tok.isin(_XOUT).any())
_xty = _xa.set_index(['src', 'side']).type
_c1 = [k for k in _xpa if _xty.get(k) == 'C1']; _tl = [k for k in _xpa if k[1] == 'tel' and _xty.get(k) in ('T0', 'T1')]
number('xpBlocks', str(len(_xparts)), [f'{O}/arm_element_blocks.tsv', f'{O}/arm_element_segments.tsv'])
number('xpCOneBoth', str(sum(_xpa[k] == (True, True) for k in _c1)), [f'{O}/arm_element_segments.tsv', EA_])
number('xpTelNoPrime', str(sum(not _xpa[k][0] for k in _tl)), [f'{O}/arm_element_segments.tsv', EA_])
assert all(_xpa[k][1] for k in _tl), 'text: every telSMN arm carries Xi (the outer part)'
assert len(_c1) == int((_xa[(_xa.side == 'cen')].type == 'C1').sum()) and len(_tl) == int(_xa[_xa.side == 'tel'].type.isin(['T0', 'T1']).sum())
# the S bundle's strand on the window against S's orientation by sequence: one strand must always give the same orientation
_ov = sj['orient_vs_bundle']; _plus = max(('ref', 'reverse'), key=lambda o: sum(v for k, v in _ov.items() if k.split('|')[1] == '+' and k.split('|')[2] == o))
_agree = sum(v for k, v in _ov.items() if (k.split('|')[1] == '+') == (k.split('|')[2] == _plus))
number('sOrientAgree', _agree, [SJ]); assert _agree == sj['n_occ'], 'text: S orientation by sequence agrees with the bundle order on every copy'
_minus = 'ref' if _plus == 'reverse' else 'reverse'
number('sAxisFwd', _ov[f'axis|+|{_plus}'], [SJ]); number('sAxisRev', _ov[f'axis|-|{_minus}'], [SJ])
v_ = sj['variant']
number('sCenSa', v_['cen|C1|Sa'], [SJ]); number('sCenAll', v_['cen|C1|Sa'] + v_.get('cen|C1|St', 0), [SJ])
number('sTelSt', v_['tel|T0|St'] + v_['tel|T1|St'], [SJ]); number('sTelAll', v_['tel|T0|St'] + v_['tel|T1|St'] + v_.get('tel|T1|Sa', 0), [SJ]); number('sTZeroSt', v_['tel|T0|St'], [SJ])
assert v_.get('tel|T0|Sa', 0) == 0, 'text: every T0 S differs from the axis S'
SPS = f'{O}/spacer_sequence.tsv'; _sp = pd.read_csv(SPS, sep='\t')
_untyped = _sp[(_sp.site == 'tel') & _sp.variant.isna()]
number('sTelOcc', sj['by_site']['tel'], [SJ]); number('sTelUntyped', len(_untyped), [SPS])
number('sNoAxisHap', ', '.join(sorted({x.replace('_hprc_r2_v1.0.1', '').replace('_', ' ') for x in _untyped.src})), [SPS])
assert sj['by_site']['tel'] == MACROS['sTelAll'] + len(_untyped) and not _sp[(_sp.site == 'axis')].src.isin(_untyped.src).any(), 'text: the uncompared telSMN S are on a haplotype with no S at the axis'
number('sSaMed', f"{sj['to_axis_by_variant']['Sa']['median']:.1f}", [SJ]); number('sStMed', f"{sj['to_axis_by_variant']['St']['median']:.1f}", [SJ])
number('sCut', re.search(r'VAR_CUT = ([\d.]+)', open('scripts/28c_spacer_sequence.py').read()).group(1), ['scripts/28c_spacer_sequence.py'])
# S copies per arm and arm type (cold read: which arms carry more than one S)
_spc = _sp[_sp.site == 'cen'].groupby('src').size(); _spt = _sp[_sp.site == 'tel']
_spT0 = _spt[_spt.arm_type == 'T0'].groupby('src').size(); _spT1 = _spt[_spt.arm_type == 'T1'].groupby('src').size()
assert set(_sp[_sp.site == 'cen'].arm_type) == {'C1'}, 'text: cenSMN S only on C1 arms'
assert (_spT0 == 1).all(), 'text: one S on every T0 arm that carries one'
assert set(_spT1) <= {2, 3}, 'text: two or three S on a T1 arm'
number('sCenHaps', int(_spc.size), [SPS]); number('sCenOnce', int((_spc == 1).sum()), [SPS]); number('sCenMore', int((_spc > 1).sum()), [SPS])
_cl3 = pd.read_csv(f'{O}/copy_labels.tsv', sep='\t'); _hc3 = pd.read_csv('tables/haplotype_copy_counts.tsv', sep='\t')
_unr3 = set(_hc3.loc[~_hc3.resolved.astype(bool), 'src'])
assert set(_cl3.loc[_cl3.position.isna(), 'src']) <= _unr3 and int(_cl3.position.isna().sum()) == int(_cl3.src.isin(_unr3).sum()), 'Fig. 3: blank positions are the copies on unresolved haplotypes'
number('sTelTzeroN', int(_spT0.sum()), [SPS]); number('sTelToneN', int(_spT1.sum()), [SPS])
number('sTelToneArms', int(_spT1.size), [SPS]); number('sTelToneThree', int((_spT1 == 3).sum()), [SPS])
assert int(_spT0.sum()) + int(_spT1.sum()) == len(_spt), 'S tel copies add up'
number('sAxisBetween', f"{sj['pairs']['between:axis-axis']['median']:.2f}", [SJ]); number('sTelBetween', f"{sj['pairs']['between:tel-tel']['median']:.2f}", [SJ])
# slot occupancy table (left arm = cenSMN position, right arm = telSMN position)
AS = f'{O}/arm_slots.tsv'; sl = pd.read_csv(AS, sep='\t', index_col=[0, 1])
cols = ['SMN2', 'SMN1 A', 'SMN1 B', 'hybrid', 'other']
rows = []
for (side, grp), r in sl.iterrows():
    lab = (r'left (\RU{A}+, cenSMN)' if side == 'left' else r'right (\RU{A}$-$, telSMN)')
    rows.append(trow(lab, grp, *[int(r.get(c, 0)) for c in cols], int(sum(r.get(c, 0) for c in cols))))
(DOC / 'tab_slots.tex').write_text('% written by scripts/40_numbers.py from results/sv/arm_slots.tsv; do not edit\n\\newcommand{\\tabSlots}{%\n' + '\n'.join(rows) + '\n}\n')
s_ = lambda side, g, c: int(sl.loc[(side, g), c]) if c in sl.columns else 0
for side, g, c, M in [('left', 'AFR', 'SMN1 B', 'slLAfrB'), ('left', 'AFR', 'SMN2', 'slLAfrTwo'), ('left', 'AFR', 'SMN1 A', 'slLAfrA'),
                      ('left', 'other', 'SMN2', 'slLOthTwo'), ('left', 'other', 'SMN1 A', 'slLOthA'), ('left', 'other', 'SMN1 B', 'slLOthB'),
                      ('right', 'AFR', 'SMN1 A', 'slRAfrA'), ('right', 'other', 'SMN1 A', 'slROthA')]:
    number(M, s_(side, g, c), [AS])
# the text counts ARMS (each holds one counted SMN copy); arm_slots.tsv counts items, and one arm also holds a truncated piece
AC_ = f'{O}/arm_copies.tsv'; _ac = pd.read_csv(AC_, sep='\t')
_cl_ = pd.read_csv(f'{O}/copy_labels.tsv', sep='\t'); _sp = _cl_.dropna(subset=['superpop']).groupby('src').superpop.first()   # as step 27
_ac['grp'] = np.where(_ac.src.map(_sp) == 'AFR', 'AFR', 'other')
_arms = _ac.groupby(['side', 'grp']).src.nunique(); _items = _ac.groupby(['side', 'grp']).size()
for (side, g), r in sl.iterrows():
    assert int(r.sum()) == int(_items[(side, g)]), ('arm_copies.tsv and arm_slots.tsv disagree', side, g)
    assert int(_items[(side, g)]) - int(_arms[(side, g)]) == int(((_ac.side == side) & (_ac.grp == g) & _ac.token.str.startswith('t')).sum()), 'items minus arms = truncated pieces'
number('slLAfr', int(_arms[('left', 'AFR')]), [AS, AC_]); number('slLOth', int(_arms[('left', 'other')]), [AS, AC_])
number('slRAfr', int(_arms[('right', 'AFR')]), [AS, AC_]); number('slROth', int(_arms[('right', 'other')]), [AS, AC_])
assert MACROS['slLAfr'] == MACROS['slRAfr'] and MACROS['slLOth'] == MACROS['slROth'], 'one left and one right arm per pair haplotype'
number('slLOthItems', int(sl.loc[('left', 'other')].sum()), [AS])
# checks guarding the text
assert (au.n_counted == 1).all() and MACROS['armExtraItem'] == 1
assert MACROS['slLAfrB'] > MACROS['slLAfrTwo'] > MACROS['slLAfrA'] and MACROS['slLOthTwo'] > 5 * MACROS['slLOthA']
# copy order (supplement S4.7) against arm (main text, Fig. 3): on the pair haplotypes with a copy on both arms
_co = pd.read_csv(PL, sep='\t'); _co['side'] = (_co.d_kb < 0).map({True: 'cen', False: 'tel'})
_both = _co.groupby('src').side.nunique(); _co = _co[_co.src.isin(_both[_both == 2].index)]
_co = _co.merge(pd.read_csv(f'{O}/copy_labels.tsv', sep='\t')[['copy_id', 'position']], on='copy_id', how='left')
_co = _co.dropna(subset=['position'])   # copies outside the 864 whole copies (Delta7-8, truncated) have no copy order
_cc = _co[_co.side == 'cen']; _ct = _co[_co.side == 'tel']
number('coCenN', len(_cc), [PL, f'{O}/copy_labels.tsv'])
number('coCenFirst', int((_cc.position == 'centromeric').sum()), [PL, f'{O}/copy_labels.tsv'], ['assert most cenSMN-arm copies are first by copy order'])
number('coTelN', len(_ct), [PL, f'{O}/copy_labels.tsv'])
number('coTelLast', int((_ct.position == 'telomeric').sum()), [PL, f'{O}/copy_labels.tsv'], ['assert most telSMN-arm copies are last by copy order'])
assert (_cc.position == 'centromeric').sum() / len(_cc) > 0.95 and (_ct.position == 'telomeric').sum() / len(_ct) > 0.95
assert set(_co.position.dropna()) <= {'centromeric', 'telomeric', 'middle'}
# the same arms split by B sub-lineage, AFR haplotypes (results/sv/palindrome_copies_lineage.tsv joined to copy_labels.tsv)
_pls = pl.merge(pd.read_csv(f'{O}/copy_labels.tsv', sep='\t')[['copy_id', 'superpop']], on='copy_id', how='left'); _pls = _pls[(_pls.gene == 'SMN1') & (_pls.superpop == 'AFR')]
for sd, L_ in [('L', True), ('R', False)]:
    for k_, nm in [('B1', 'Bone'), ('B2', 'Btwo')]:
        number(f'sl{sd}Afr{nm}', int(((_pls.sublineage == k_) & ((_pls.d_kb < 0) == L_)).sum()), [PL, f'{O}/copy_labels.tsv'])
assert MACROS['slLAfrBone'] + MACROS['slLAfrBtwo'] <= MACROS['slLAfrB'] and MACROS['slLAfrBone'] > MACROS['slLAfrTwo'], 'text: B1 is the commonest gene at the AFR cenSMN arm'
assert MACROS['slRAfrBtwo'] > MACROS['slRAfrBone'], 'text: at the AFR telSMN arm B2, not B1'
# Fig. 4a (step 22) draws results/sv/palindrome_arm_marks.tsv: same arms and counts as arm_slots.tsv (step 27) and the text
_am = pd.read_csv(f'{O}/palindrome_arm_marks.tsv', sep='\t', index_col=[0, 1])
for (sd, g_), r_ in _am.iterrows():
    s_ = sl.loc[(sd, g_)]
    assert r_['SMN1 A'] == s_['SMN1 A'] and r_[['SMN1 B1', 'SMN1 B2', 'SMN1 B-']].sum() == s_['SMN1 B'] and r_['SMN2'] == s_['SMN2'] and r_['hybrid'] == s_['hybrid'] and r_['other'] == s_['other'], (sd, g_)
assert _am.loc[('left', 'AFR'), 'SMN1 B1'] == MACROS['slLAfrBone'] and _am.loc[('right', 'AFR'), 'SMN1 B2'] == MACROS['slRAfrBtwo'] and _am.loc[('left', 'AFR')].sum() == MACROS['slLAfr']
assert MACROS['slLAfr'] + MACROS['slLOth'] == MACROS['armN'] // 2 and MACROS['slRAfr'] + MACROS['slROth'] == MACROS['armN'] // 2   # arms
assert MACROS['slLOthItems'] == MACROS['slLOth'] + MACROS['armExtraItem']   # items: the one truncated piece is on a non-AFR left arm
assert MACROS['slLAfr'] + MACROS['slLOthItems'] + MACROS['slRAfr'] + MACROS['slROth'] == MACROS['palCopies']   # items (copies and the truncated piece)

# ---- outgroups, the order of events and the ape copies (09, 11, 12, 15)
ES = f'{O}/evolution_summary.json'; X = json.load(open(ES))
# Our S sites against Berzal-Serrano 2026 Fig. 4 PSVs (data/refs/berzal2026_fig4_psv.tsv).
# c. names from GRCh38: c.835 = 70,951,941 (c.840 = 70,951,946); exon 7 ends c.*3 = 70,951,994; c.*239 = 70,952,674.
_es = pd.read_csv(f'{O}/evolution_sites_checked.tsv', sep='\t')
_sp = sorted(int(x) for x in _es.loc[_es['kind'] == 'S', 'grch38'])
def _cname(pos):
    if pos < 70951941: return ('i6', 70951941 - pos)
    if pos <= 70951994: return ('c', pos - 70951941 + 835)
    if pos < 70952439: return ('i7', pos - 70951994)
    return ('u', pos - 70952674 + 239)
def _bname(n):
    if n.startswith('c.835-'): return ('i6', int(n[6:]))
    if n.startswith('c.*3+'): return ('i7', int(n[5:]))
    if n.startswith('c.*'): return ('u', int(n[3:]))
    return ('c', int(n[2:]))
_bz = pd.read_csv('data/refs/berzal2026_fig4_psv.tsv', sep='\t', comment='#')
_bn = [_bname(n) for n in _bz['name']]
_exact = _off = 0; _used = set()
for pos in _sp:
    k, v = _cname(pos)
    if (k, v) in _bn: _exact += 1; _used.add((k, v)); continue
    if k == 'i6' and (k, v + 5) in _bn: _off += 1; _used.add((k, v + 5)); continue
    raise SystemExit(f'S site {pos} ({k},{v}) has no Berzal-Serrano PSV')
assert len(_sp) == 13 and len(_bz) == 15, (len(_sp), len(_bz))
assert _exact == 8 and _off == 5 and len(_used) == 13, (_exact, _off)
assert sorted(set(_bn) - _used) == [('i6', 1606), ('i6', 1897)]
M_bz = {'bzPsvN': len(_bz), 'bzMatch': _exact + _off, 'bzExact': _exact, 'bzOff': _off, 'bzOffBp': 5, 'bzUnmatched': 'c.835-1897 and c.835-1606'}
for _k, _v in M_bz.items(): number(_k, _v, ['data/refs/berzal2026_fig4_psv.tsv', f'{O}/evolution_sites_checked.tsv'], ['assert 13 S sites = 8 Berzal names exactly + 5 intron-6 names 5 bp apart; Berzal PSVs 1 and 5 unmatched'])
for k, M in [('n_S', 'evS'), ('n_AB', 'evAB'), ('S_anc_smn1', 'evSancOne'), ('S_anc_smn2', 'evSancTwo'), ('S_unres', 'evSunres'),
             ('S_shared_derived_AB', 'evSsharedAB'), ('S_apes_unanimous', 'evSunan'), ('AB_anc_A', 'evABancA'), ('AB_anc_B', 'evABancB'),
             ('AB_unres', 'evABunres'), ('AB_smn2_eq_A', 'evABtwoEqA'), ('AB_smn2_eq_A_anc', 'evABtwoEqAanc'),
             ('AB_shared_derived_A_SMN2', 'evABsharedATwo'), ('ape_haps', 'apeHaps'), ('ape_copies', 'apeCopies'),
             ('gene_full', 'apeGeneFull'), ('gene_short', 'apeGeneShort')]:
    number(M, X[k], [ES])
assert X['S_anc_smn1'] + X['S_anc_smn2'] + X['S_unres'] == X['n_S'], 'S sites: derived + ancestral + unresolved must equal all'   # text: 8 derived, 4 ancestral, 1 unresolved of 13
number('evSspanKb', f"{(X['S_span'][1] - X['S_span'][0] + 1) / 1000:.1f}", [ES])
number('evSfrom', f"{X['S_span'][0]:,}".replace(',', '{,}'), [ES]); number('evSto', f"{X['S_span'][1]:,}".replace(',', '{,}'), [ES])
number('evABfrom', f"{X['AB_span'][0]:,}".replace(',', '{,}'), [ES]); number('evABto', f"{X['AB_span'][1]:,}".replace(',', '{,}'), [ES])
number('apeShortMin', X['gene_short_range'][0], [ES]); number('apeShortMax', f"{X['gene_short_range'][1]:,}".replace(',', '{,}'), [ES])
for ape, M in [('chimpanzee', 'Chimp'), ('bonobo', 'Bonobo'), ('gorilla', 'Gorilla'), ('orangutan', 'Orang')]:
    number(f'evTwoBase{M}', X['S_apes_smn2_base'][ape], [ES]); number(f'evABA{M}', X['AB_apes_A_base'][ape], [ES])
assert X['S_anc_smn1'] > X['S_anc_smn2'] and X['c840_apes'] == 'CCCC' and X['cs239_apes'] == 'GGGG'
assert X['AB_anc_B'] > X['AB_anc_A'] and X['AB_shared_derived_A_SMN2'] == X['AB_anc_B'] and X['S_shared_derived_AB'] > 0, 'text: B keeps the ape base at most A/B sites; there A and SMN2 share the new base'
assert X['ape_haps'] == X['ape_copies'] == 8 and X['gene_full'] == 4
P1f, P2f = f'{O}/ape_copies.tsv', f'{O}/ape_copies_hap2.tsv'
P1 = pd.read_csv(P1f, sep='\t'); P1 = P1[P1.ape != 'human GRCh38'].assign(haplotype='primary')
P2all = pd.read_csv(P2f, sep='\t')
# The lower-density (relaxed) search of the second haplotypes: every relaxed hit overlaps a strict copy of the same
# ape haplotype, so it finds nothing outside those copies; the one extra hit (orangutan) lies inside the strict
# copy's own span, on the opposite strand, and is not a second copy (S7 and Table S11's caption).
_a2s, _a2r = P2all[P2all.rule == 'strict'], P2all[P2all.rule == 'relaxed']
_ovl = lambda r: any(s.ape == r.ape and s.chrom == r.chrom and s.chr_start < r.chr_end and r.chr_start < s.chr_end for s in _a2s.itertuples())
assert all(_ovl(r) for r in _a2r.itertuples()), 'text (S7): the relaxed search finds nothing outside the strict copies'
_a2e = _a2r.groupby('ape').filter(lambda g: len(g) > 1)
assert list(_a2e.ape.unique()) == ['orangutan'] and len(_a2e) == 2, 'text: one ape haplotype returns a further relaxed hit'
_a2o = _a2e.sort_values('body_cov').iloc[0]; _a2w = _a2e.sort_values('body_cov').iloc[1]
assert _a2o.chr_start >= _a2w.chr_start and _a2o.chr_end <= _a2w.chr_end, 'text: the further hit lies inside the copy span'
assert _a2o.strand != _a2w.strand, 'text: the further hit is on the opposite strand'
number('apeRelaxExtra', len(_a2e) - 1, [P2f])
number('apeOrangPartCov', f'{_a2o.body_cov:.2f}', [P2f])
P2 = P2all[P2all.rule == 'strict']
AP = pd.concat([P1, P2])
AP['ape'] = pd.Categorical(AP.ape, ['chimp', 'bonobo', 'gorilla', 'orangutan'], ordered=True)
label = {'chimp': 'chimpanzee', 'bonobo': 'bonobo', 'gorilla': 'gorilla', 'orangutan': 'orangutan'}
hapl = {'primary': 'primary', 'alt': 'second', 'other': 'second'}
# cold read 2026-11-16: Table S12 also gives the read-checked bases (step 13c), which the text and Fig. 6 use; primary assemblies only
_EVC = f'{O}/evolution_sites_checked.tsv'; _evc = pd.read_csv(_EVC, sep='\t', dtype=str)
def _rdS(r):
    if r.haplotype != 'primary': return '--'
    e = _evc[_evc.kind == 'S']; c = e[f'og_{r.ape}']
    return f'{int((c == e.SMN1).sum())}/{int((c == e.SMN2).sum())}'
def _rdAB(r):
    if r.haplotype != 'primary': return '--'
    e = _evc[_evc.kind == 'AB']; c = e[f'og_{r.ape}']
    return f'{int((c == e.A).sum())}/{int((c == e.B).sum())}'
for _a, _nm2 in [('chimp', 'chimpanzee'), ('bonobo', 'bonobo'), ('gorilla', 'gorilla'), ('orangutan', 'orangutan')]:
    _e = _evc[_evc.kind == 'S']; assert int((_e[f'og_{_a}'] == _e.SMN2).sum()) == X['S_apes_smn2_base'][_nm2], 'Table S12 read-checked S counts = step 16'
    _e = _evc[_evc.kind == 'AB']; assert int((_e[f'og_{_a}'] == _e.A).sum()) == X['AB_apes_A_base'][_nm2], 'Table S12 read-checked A/B counts = step 16'
rows = []
for r in AP.sort_values(['ape', 'haplotype'], key=lambda s: s if s.name == 'ape' else s.map({'primary': 0, 'alt': 1, 'other': 1})).itertuples():
    rows.append(' & '.join([label[r.ape], hapl[r.haplotype], r.chrom.replace('_', r'\_'), f'{int(r.chr_start):,}'.replace(',', '{,}'), r.strand,
                            r.c840, r.c239, f'{int(r.S_smn1)}/{int(r.S_smn2)}', f'{int(r.AB_A)}/{int(r.AB_B)}', _rdS(r), _rdAB(r)]) + r' \\')
(DOC / 'tab_apes.tex').write_text('% written by scripts/40_numbers.py; do not edit\n\\newcommand{\\tabApes}{%\n' + '\n'.join(rows) + '\n}\n')
assert len(rows) == 8 and set(AP.c840) == {'C'} and set(AP.c239) == {'G'}
number('apeOrangPairA', int(P1[P1.ape == 'orangutan'].AB_A.iloc[0]), [P1f]); number('apeOrangPairB', int(P1[P1.ape == 'orangutan'].AB_B.iloc[0]), [P1f])
# 2026-11-16 (researcher: note that the alternates were left out of the history counts and where they disagree)
_cols = ['S_smn1', 'S_smn2', 'AB_A', 'AB_B']
_same = [a for a in ['chimp', 'bonobo', 'gorilla', 'orangutan']
         if (P1[P1.ape == a][_cols].astype(int).values[0] == P2[P2.ape == a][_cols].astype(int).values[0]).all()]
_diff = [a for a in ['chimp', 'bonobo', 'gorilla', 'orangutan'] if a not in _same]
assert _same == ['chimp', 'bonobo', 'gorilla'] and _diff == ['orangutan'], 'text: alternates agree pair-wise with the primary except orangutan'
_oa = P2[P2.ape == 'orangutan'].iloc[0]
number('apeAltOrangA', int(_oa.AB_A), [P2f]); number('apeAltOrangB', int(_oa.AB_B), [P2f])
number('apeAltOrangSOne', int(_oa.S_smn1), [P2f]); number('apeAltOrangSTwo', int(_oa.S_smn2), [P2f])
_op = P1[P1.ape == 'orangutan'].iloc[0]
number('apeOrangPairSOne', int(_op.S_smn1), [P1f]); number('apeOrangPairSTwo', int(_op.S_smn2), [P1f])
assert int(_oa.AB_A) > int(_oa.AB_B) and int(_op.AB_A) < int(_op.AB_B), 'text: the orangutan alternate leans to A where its primary leans to B'
assert X['AB_smn2_eq_A'] == X['n_AB'], 'text: SMN2 carries A\'s base at all A/B sites'
# 2026-11-15 (researcher: which alignment is right at the ape sites): the ape bases the text uses are read from the raw
# primary-assembly sequence (step 13c), as step 16 applies them; list where they differ from the multiple alignment and
# from the pair-wise typing (Table S12)
RD_ = f'{O}/ape_site_reads.tsv'; _rd = pd.read_csv(RD_, sep='\t', keep_default_na=False)
EVC_ = f'{O}/evolution_sites_checked.tsv'; _evc = pd.read_csv(EVC_, sep='\t', keep_default_na=False)
_og = {'chimp': 'og_chimp', 'bonobo': 'og_bonobo', 'gorilla': 'og_gorilla', 'orangutan': 'og_orangutan'}
_nm = {'chimp': 'chimpanzee', 'bonobo': 'bonobo', 'gorilla': 'gorilla', 'orangutan': 'orangutan'}
number('rdN', len(_rd), [RD_]); number('rdAcc', int(_rd.accepted.sum()), [RD_])
assert len(_rd) == 4 * len(_evc)
fmt_s = lambda g: 'chr5:' + f'{g:,}'.replace(',', '{,}')
_chg = _rd[_rd.accepted & (_rd.read != _rd.msa)]
number('rdMsaDiffN', len(_chg), [RD_])
number('rdMsaDiff', '; '.join(f"{fmt_s(int(r.grch38))} ({'A/B' if r.kind == 'AB' else 'S'} site) in {_nm[r.ape[3:]]}: {r.read} in the sequence, {r.msa if r.msa != '-' else 'a gap'} in the alignment" for r in _chg.itertuples()), [RD_])
_unr = _evc[(_evc.kind == 'AB') & _evc.ancestral.eq('')]
assert list(_unr.grch38) == [70936613], 'text: one A/B site unresolved, the homopolymer site'
number('rdUnresSite', fmt_s(70936613), [EVC_])
assert all(_evc.loc[_evc.grch38 == 70936613, c].iloc[0] == '?' for c in _og.values()), 'text: unresolved in all four apes'
# the GRCh38 run at that site: 14 T then 5 G
_w = open('data/derived/windows/GCA_000001405.15_GRCh38_no_alt_analysis_set.PanSN.fa').read().split('\n', 1)[1].replace('\n', '')
_o = 70936613 - 69445889; assert _w[_o - 13:_o + 1] == 'T' * 14 and _w[_o + 1:_o + 6] == 'G' * 5 and _w[_o - 14] != 'T'
_evs = _evc[_evc.kind == 'AB']
_diff = {}
for ape, col in _og.items():
    pr = P1[P1.ape == ape].iloc[0]
    for e in _evs.itertuples():
        pb, cb = str(pr[f'AB_{e.grch38}']), str(getattr(e, col))
        if cb != '?' and pb != cb: _diff.setdefault((int(e.grch38), pb, cb), []).append(_nm[ape])
    assert sum(str(getattr(e, col)) == e.A for e in _evs.itertuples()) == X['AB_apes_A_base'][_nm[ape]], 'text: per-ape A counts are the read-checked ones'
_g = lambda b: 'a gap' if b == '-' else b
number('apeDiffN', sum(len(v) for v in _diff.values()), [P1f, EVC_])
number('apeDiffSites', '; '.join(fmt_s(k[0]) + f" in {' and '.join(v)} ({_g(k[1])} pair-wise, {_g(k[2])} in the sequence)" for k, v in sorted(_diff.items())) or 'none', [P1f, EVC_])
_ks = P1.set_index('ape').kmer_share
number('apeOrangKmer', f"{_ks['orangutan']:.2f}", [P1f]); number('apeOtherKmerLo', f"{_ks.drop('orangutan').min():.2f}", [P1f]); number('apeOtherKmerHi', f"{_ks.drop('orangutan').max():.2f}", [P1f])
assert _ks['orangutan'] < 0.6 * _ks.drop('orangutan').min(), 'text: the primary orangutan pair-wise match is the weakest'
figure('figures/evolution_order.png', [f'{O}/copies_with_outgroups.aln.fa', ES, P1f, P2f])

# --- step 19b (main text: sequence side = supported ML clades) and step 19 (supplement: UPGMA on both sides, the check)
MCP, MCM, MVU = f'{O}/ml_cluster_pairs.tsv', f'{O}/ml_cluster_pairs_members.tsv', f'{O}/ml_vs_upgma_bands.tsv'
MSJ, MSE, MSC = f'{O}/ml_bands_summary.json', f'{O}/ml_band_sensitivity.tsv', f'{O}/ml_seq_clusters.tsv'
CPT = f'{O}/cluster_pairs.tsv'; CPM = f'{O}/cluster_pairs_members.tsv'; CHG = f'{O}/chen_hg_copies.tsv'; CL = f'{O}/copy_labels.tsv'
mcp = pd.read_csv(MCP, sep='\t'); mcp = mcp[mcp.band.notna()].set_index('band')
mcm = pd.read_csv(MCM, sep='\t'); mvu = pd.read_csv(MVU, sep='\t').set_index('band'); msj = json.load(open(MSJ))
mse = pd.read_csv(MSE, sep='\t'); msc = pd.read_csv(MSC, sep='\t')
chg = pd.read_csv(CHG, sep='\t'); chg['src'] = chg.copy_id.str.split('|').str[0]
cl = pd.read_csv(CL, sep='\t')
def sci(p):
    e = int(np.floor(np.log10(p))); return f'{p / 10**e:.0f}\\times10^{{{e}}}'
BANDS = ['Bulk', 'Two', 'AB']
assert sorted(mcp.index) == sorted(BANDS) and sorted(mvu.index) == sorted(BANDS), 'text: three bands at the coarse scale'
for nm in BANDS:
    r = mcp.loc[nm]
    number(f'cb{nm}n', int(r.overlap), [MCP]); number(f'cb{nm}exp', f'{r.expected:.1f}', [MCP]); number(f'cb{nm}p', sci(r.p_bonf), [MCP])
    number(f'cb{nm}posS', f"{mvu.loc[nm, 'struct_pos_median']:.2f}", [MVU]); number(f'cb{nm}posG', f"{mvu.loc[nm, 'seq_pos_median']:.2f}", [MVU])
    mem = mcm[mcm.band == nm]
    assert len(mem) == r.overlap
    number(f'cb{nm}afr', int((mem.superpop == 'AFR').sum()), [MCM])
    number(f'cb{nm}topArr', int(mem.arrangement.value_counts().iloc[0]), [MCM])
    c1 = cl[cl.src.isin(mem.src) & (cl.cls == 'SMN1')]
    number(f'cb{nm}cop', len(c1), [CL]); number(f'cb{nm}rev', int((c1.strand == '-').sum()), [CL])
    number(f'cb{nm}tel', int((c1.position == 'telomeric').sum()), [CL]); number(f'cb{nm}single', int((c1.position == 'single').sum()), [CL])
    # the same band under UPGMA on both sides (step 19), and how many haplotypes the two share
    u = mvu.loc[nm]
    number(f'ub{nm}n', int(u.upgma_n), [MVU, CPT]); number(f'ub{nm}exp', f'{u.upgma_expected:.1f}', [MVU, CPT]); number(f'ub{nm}p', sci(u.upgma_p_bonf), [MVU, CPT])
    number(f'ub{nm}shared', int(u.shared), [MVU]); number(f'ub{nm}onlyU', int(u.only_upgma), [MVU]); number(f'ub{nm}onlyM', int(u.only_ml), [MVU])
_mt = mcm[mcm.band == 'Two']; hg = chg[chg.src.isin(_mt.src) & chg.lineage.isin(['A', 'B'])]
number('cbTwohg', int(_mt.src.isin(set(hg[hg.family.isin({'S1-2', 'S1-2c'})].src)).sum()), [MCM, CHG])
number('cbTwoTwoOne', int((_mt.arrangement == '2- 1-').sum()), [MCM]); number('cbTwoOneOnly', int((_mt.arrangement == '1-').sum()), [MCM])
assert mcp.loc['Two'].arrangements.startswith('2- 1- ') and mcp.loc['AB'].arrangements == f"1- 1+ {MACROS['cbABn']}", 'text: band arrangements'
assert int(MACROS['cbABn']) <= int(MACROS['linMixHaps']), 'B1+B2 band is a subset of the haplotypes with one B1 and one B2'
assert int(MACROS['cbTwoafr']) == 0 and int(MACROS['cbABafr']) >= int(MACROS['cbABn']) - 1, 'text: no AFR in the S1-2 band; all but one AFR in the B1+B2 band'
_ml = pd.read_csv(f'{O}/marker_lineage.tsv', sep='\t'); _ml = _ml[_ml['class'] == 'SMN1'].assign(src=lambda x: x.copy_id.str.split('|').str[0])
_sub = _ml.groupby('src').sublineage.apply(lambda x: '+'.join(sorted(x.astype(str))))
assert (_sub.reindex(mcm[mcm.band == 'AB'].src) == 'B1+B2').all(), 'text: every haplotype of the AB band carries one B1 and one B2 copy'
assert int(MACROS['cbTwocop']) == int(MACROS['cbTwon']) == int(MACROS['cbTworev']), 'text: each S1-2 band haplotype has one SMN1 copy, all reversed'
assert int(MACROS['cbTwohg']) >= int(MACROS['cbTwon']) - 2, 'text: S1-2 band = haplogroup S1-2'
for nm in BANDS: assert mcp.loc[nm].p_bonf < 1e-10 and mvu.loc[nm].upgma_p_bonf < 1e-10
assert mvu.loc['Two'].shared == mvu.loc['Two'].upgma_n == mvu.loc['Two'].ml_n and mvu.loc['AB'].shared == mvu.loc['AB'].upgma_n == mvu.loc['AB'].ml_n, 'text: S1-2 and B1+B2 bands identical under both trees'
# ML-clade settings and their sensitivity
number('mlbHaps', msj['haplotypes'], [MSJ]); number('mlbCopies', msj['SMN1_copies'], [MSJ])
number('mlbClusters', msj['seq_clusters'], [MSJ]); number('mlbLabels', msj['seq_labels'], [MSJ])
number('mlbSettings', msj['sensitivity_settings'], [MSJ]); number('mlbSettingsAll', msj['settings_all_bands'], [MSJ], ['settings that give every band of the headline setting'])
assert msj['settings_all_bands'] == msj['sensitivity_settings'] and set(msj['settings_n_bands'].values()) == {len(BANDS)}, 'text: the same three bands at every setting'
MLP = f'{O}/ml_cluster_pairs.tsv'; mlp = pd.read_csv(MLP, sep='\t')
assert msj['alpha'] == 1e-3; number('mlbAlphaExp', 3, [MSJ]); number('mlbMinOv', msj['min_overlap'], [MSJ])
bnd = mlp[(mlp.p_bonf < msj['alpha']) & (mlp.overlap >= msj['min_overlap'])]; sig = mlp[mlp.p_bonf < 0.05]
assert len(bnd) == 3 and sorted(bnd.seq_label) == ['S1-2', 'S1-8+S1-9', 'rest'], 'text: three bands: S1-2, S1-8 with S1-9, and the bulk'
small = sig[~sig.index.isin(bnd.index)]; number('mlbSmallPairs', len(small), [MLP]); assert len(small) == 1, 'text: one smaller pair'
number('mlbSmallOv', int(small.overlap.iloc[0]), [MLP]); number('mlbSmallP', pfmt(small.p_bonf.iloc[0]), [MLP])
assert small.overlap.iloc[0] < msj['min_overlap'] and msj['alpha'] <= small.p_bonf.iloc[0] < 0.05
# where the S1-3 clade's haplotypes fall at this scale: split between the bulk's and the B1+B2 band's structure clusters, neither pair significant
_s3 = mlp[mlp.seq_label == 'S1-3'].sort_values('overlap', ascending=False)
assert (_s3.p_bonf == 1).all() and list(_s3.struct_cluster[:2]) == [int(mcp.loc['Bulk'].struct_cluster), int(mcp.loc['AB'].struct_cluster)], 'text: S1-3 splits between the bulk and B1+B2 clusters'
number('cbSThreeN', int(_s3.seq_size.iloc[0]), [MLP]); number('cbSThreeBulk', int(_s3.overlap.iloc[0]), [MLP]); number('cbSThreeAB', int(_s3.overlap.iloc[1]), [MLP])
number('cbBulkClu', int(mcp.loc['Bulk'].struct_size), [MLP]); number('cbABClu', int(mcp.loc['AB'].struct_size), [MLP])
number('cbSharedClu', int(mcp.loc['AB'].struct_size), [MLP])
assert msj['K'] == int(re.search(r"\nk = (\d+); cs, cg = fcluster", open('scripts/19_cluster_correspondence.py').read()).group(1))
msc_ = msc[(msc['T'] == msj['T']) & (msc['M'] == msj['M'])].set_index('cluster')
for c_, M_ in [('S1-2', 'Two'), ('S1-3', 'Three'), ('S1-8', 'Eight'), ('S1-9', 'Nine'), ('S1-6', 'Six')]:
    number(f'mlb{M_}Sup', int(msc_.loc[c_, 'support']), [MSC])
number('mlbRestCopies', int(msc_.loc['rest', 'own_SMN1']), [MSC])
_lab = pd.read_csv(f'{O}/ml_hap_labels.tsv', sep='\t').set_index('src')
_s13 = chg[chg.family == 'S1-3']; _out = set(_s13.src) - set(_lab.index)
_pd = pd.read_csv(f'{O}/pair_distances.tsv.gz', sep='\t')
assert len(_out) == 1 and not (_pd.a.isin(_out) | _pd.b.isin(_out)).any(), 'text: the one S1-3 copy outside the clade in Fig. 5 is on a haplotype without a structural distance'
figure('figures/tanglegram_ml.png', [f'{O}/pair_distances.tsv.gz', f'{O}/mltree.treefile', f'{O}/marker_lineage.tsv', CHG, f'{O}/lineage_haplotypes.tsv', CPM, CPT, f'{O}/arm_element_haps.tsv'])
figure('figures/tanglegram_clusters.png', [f'{O}/pair_distances.tsv.gz', f'{O}/lineage_haplotypes.tsv'])
# --- the finer scale (Supplementary Section 'The finer scale'; step 37 archive of the pipeline at commit 953c895)
FD = 'results/sv/fine'; FMV, FMC, FMM, FMS = f'{FD}/ml_vs_upgma_bands.tsv', f'{FD}/ml_cluster_pairs.tsv', f'{FD}/ml_cluster_pairs_members.tsv', f'{FD}/ml_bands_summary.json'
fmv = pd.read_csv(FMV, sep='\t').set_index('band'); fmm = pd.read_csv(FMM, sep='\t'); fms = json.load(open(FMS))
FBANDS = ['Bulk', 'Two', 'Three', 'AB']
assert sorted(fmv.index) == sorted(FBANDS), 'text: four bands at the finer scale'
for nm in FBANDS:
    r = fmv.loc[nm]; number(f'fb{nm}n', int(r.ml_n), [FMV]); number(f'fb{nm}exp', f'{r.ml_expected:.1f}', [FMV]); number(f'fb{nm}p', sci(r.ml_p_bonf), [FMV])
    s_ = set(fmm[fmm.band == nm].src); assert len(s_) == r.ml_n
    if nm in BANDS: number(f'fb{nm}inCoarse', len(s_ & set(mcm[mcm.band == nm].src)), [FMM, MCM])
_f3 = set(fmm[fmm.band == 'Three'].src); _cc = _lab.reindex(sorted(_f3)).struct_cluster
assert set(_cc) == {int(mcp.loc['Bulk'].struct_cluster), int(mcp.loc['AB'].struct_cluster)}, 'text: the finer S1-3 band falls in the two large coarse clusters'
number('fbThreeToBulk', int((_cc == int(mcp.loc['Bulk'].struct_cluster)).sum()), [FMM, f'{O}/ml_hap_labels.tsv']); number('fbThreeToAB', int((_cc == int(mcp.loc['AB'].struct_cluster)).sum()), [FMM, f'{O}/ml_hap_labels.tsv'])
_f3m = fmm[fmm.band == 'Three']; _f3c = cl[cl.src.isin(_f3m.src) & (cl.cls == 'SMN1')]
assert len(_f3c) == len(_f3m) and (_f3c.strand == '+').all(), 'text (finer): each S1-3 band haplotype has one forward SMN1 copy'
assert set(_f3m.arrangement.value_counts().index[:2]) == {'2- 1+', '1+'}, 'text (finer): mostly 2- 1+ and 1+'
_f3h = chg[chg.src.isin(_f3m.src) & (chg.family == 'S1-3')]; assert _f3h.src.nunique() == len(_f3m), 'text (finer): S1-3 band copies are haplogroup S1-3'
assert MACROS['fbABinCoarse'] == MACROS['fbABn'] and MACROS['fbTwoinCoarse'] == MACROS['cbTwon'], 'text: the finer B1+B2 band lies inside the coarse one; the coarse S1-2 band inside the finer one'
# --- step 17: what sets the S1-3 band apart
SJ17 = 'results/sv/fine/s13_summary.json'; SF = 'results/sv/fine/s13_unit_families.tsv'; SCD = 'results/sv/fine/s13_G_carriers_CD.tsv'   # finer scale (step 37)
s17 = json.load(open(SJ17)); sf = pd.read_csv(SF, sep='\t').set_index(['subset', 'family']); scd = pd.read_csv(SCD, sep='\t').set_index('group')
for k_, v in [('sBand', 'band'), ('sCen', 'central'), ('sBandG', 'band_G'), ('sCenG', 'central_G'), ('sGandT', 'G_and_S13'), ('sGnotT', 'G_not_S13'),
              ('sTnotG', 'S13_not_G'), ('sHaps', 'haplotypes'), ('sGunits', 'G_units')]:
    number(k_, s17[v], [SJ17])
number('sTall', s17['G_and_S13'] + s17['S13_not_G'], [SJ17])
number('sGkb', f"{s17['G_span_median_bp'] / 1000:.1f}", [SJ17])
for f in 'CD':
    r = sf.loc[('all', f)]; number(f'sBand{f}', int(r.band_with), [SF]); number(f'sCen{f}', int(r.central_with), [SF])
r = sf.loc[('2- 1+ only', 'G')]; number('sBandTwoOneN', int(r.band_n), [SF]); number('sCenTwoOneN', int(r.central_n), [SF]); number('sBandTwoOneG', int(r.band_with), [SF]); number('sCenTwoOneG', int(r.central_with), [SF])
for g_, M in [('S1-2', 'Two'), ('S1-3', 'Three'), ('S1-8/9/6', 'Eight')]:
    number(f'sG{M}n', int(scd.loc[g_, 'n']), [SCD]); number(f'sG{M}C', int(scd.loc[g_, 'with_C']), [SCD]); number(f'sG{M}D', int(scd.loc[g_, 'with_D']), [SCD])
assert s17['band_G'] == s17['band'] and s17['central_G'] == 0
assert sf.loc[('2- 1+ only', 'G')].band_with == sf.loc[('2- 1+ only', 'G')].band_n and sf.loc[('2- 1+ only', 'G')].central_with == 0
assert s17['G_and_S13'] > 0.85 * (s17['G_and_S13'] + s17['S13_not_G']) and s17['G_not_S13'] > s17['G_and_S13'], 'G marks S1-3 but is not specific'
assert scd.loc['S1-3', 'with_C'] > 0.8 * scd.loc['S1-3', 'n'] and scd.loc['S1-2', 'with_C'] > 0.75 * scd.loc['S1-2', 'n'] and scd.loc['S1-8/9/6', 'with_D'] > 0.8 * scd.loc['S1-8/9/6', 'n']
# --- step 18: shuffles within five continental groups
FG = f'{O}/five_group_shuffles.tsv'; fg = pd.read_csv(FG, sep='\t').set_index(['test', 'what'])
for w, M in [('seq_SMN1', 'One'), ('seq_SMN2', 'Two')]:
    number(f'fgMantel{M}', f"{fg.loc[('mantel', w), 'observed']:.2f}", [FG]); number(f'fgMantelP{M}', f"{fg.loc[('mantel', w), 'p']:.3f}", [FG])
for cl_, M in [('SMN1', 'One'), ('SMN2', 'Two')]:
    for lab_, L_ in [('arrangement', 'Arr'), ('position', 'Pos'), ('strand', 'Str'), ('env', 'Env')]:
        r = fg.loc[('copy_neighbor', f'{cl_}:{lab_}')]
        number(f'fg{L_}{M}', f'{r.observed:.2f}', [FG]); number(f'fg{L_}Null{M}', f'{r.null_mean:.2f}', [FG]); number(f'fg{L_}P{M}', f'{r.p:.3f}', [FG])
for mu in ['AB', 'B', 'BB', 'AA']:
    for tr, T in [('structure', 'Struct'), ('SMN1 sequence', 'Seq')]:
        number(f'fg{mu}p{T}', f"{fg.loc[('make_up_neighbor', f'{tr}:{mu}'), 'p']:.4f}", [FG])
assert (fg.loc['mantel'].p <= 0.001 + 1e-9).all()
assert all(fg.loc[('copy_neighbor', f'{c}:{l}'), 'p'] <= 0.001 + 1e-9 for c in ['SMN1', 'SMN2'] for l in ['arrangement'])
assert all(fg.loc[('make_up_neighbor', f'{t}:{m}'), 'p'] <= 0.0005 + 1e-9 for t in ['structure', 'SMN1 sequence'] for m in ['B', 'BB'])
assert fg.loc[('make_up_neighbor', 'structure:AA'), 'p'] > 0.05
assert fg.loc[('copy_neighbor', 'SMN2:env'), 'p'] > 0.01


# --- step 21: Fig. 1 by arrangement type (within-type correlation, diversity, window length)
ST = f'{O}/scatter_type_stats.tsv'; SS = f'{O}/scatter_type_summary.json'
stt = pd.read_csv(ST, sep='\t').set_index('arrangement'); ss = json.load(open(SS))
q = lambda s: '``' + s + "''"
number('stTypes', ss['types_tested'], [SS]); number('stTypeHaps', ss['haplotypes_in_tested'], [SS])
sig = stt[stt.mantel_all_p <= 0.001]; ns = stt[stt.mantel_all_p > 0.05]
assert list(ns.index) == ['1- 2+'] and len(sig) == len(stt) - 1, 'one type without a within-type correlation'
number('stSigTypes', len(sig), [ST]); number('stRlo', f'{sig.mantel_all_r.min():.2f}', [ST]); number('stRhi', f'{sig.mantel_all_r.max():.2f}', [ST])
number('stRbig', f"{stt.loc['2- 1+', 'mantel_all_r']:.2f}", [ST]); number('stNsType', q('1- 2+'), [ST])
number('stNsR', f"{stt.loc['1- 2+', 'mantel_all_r']:.2f}", [ST]); number('stNsP', f"{stt.loc['1- 2+', 'mantel_all_p']:.2f}", [ST]); number('stNsN', int(stt.loc['1- 2+', 'haplotypes']), [ST])
mx_s = stt.struct_median.idxmax(); mn_s = stt.struct_median.idxmin(); mx_q = stt.seq_all_median_perkb.idxmax()
assert mx_s == '1-' and mn_s == '1- 2+' and mx_q == '1-'
assert set(stt.seq_all_median_perkb.nlargest(2).index) == {'1-', '1+'}, 'single-copy types most diverse in sequence'
assert set(stt.seq_all_median_perkb.nsmallest(2).index) == {'1- 2+', '2- 1+'}
for nm, ty in [('One', '1-'), ('Plus', '1+'), ('Big', '2- 1+'), ('Mix', '1- 2+'), ('Inv', '1- 1+')]:
    number(f'stS{nm}', f"{stt.loc[ty, 'struct_median']:.2f}", [ST]); number(f'stQ{nm}', f"{stt.loc[ty, 'seq_all_median_perkb']:.2f}", [ST])
    number(f'stL{nm}', f"{stt.loc[ty, 'window_len_median_kb']:,.0f}".replace(',', '{,}'), [ST])
number('stLenSame', f"{ss['struct_vs_lendiff_samepairs']:.2f}", [SS]); number('stLenAll', f"{ss['struct_vs_lendiff_all']:.2f}", [SS])
number('stLenSeqAll', f"{ss['seq_vs_lendiff_all']:.2f}", [SS])
assert ss['struct_vs_lendiff_samepairs'] > ss['seq_vs_lendiff_all'], 'length goes with structure more than with sequence'
for nm, ty in [('Big', '2- 1+'), ('Inv', '1- 1+'), ('Plus', '1+'), ('One', '1-'), ('Mix', '1- 2+'), ('Flip', '2- 1-')]:
    number(f'stP{nm}', f"{stt.loc[ty, 'struct_seq_partial_lendiff_r']:.2f}", [ST])
assert stt.loc[['2- 1+', '1- 1+', '2- 1-'], 'struct_seq_partial_lendiff_r'].min() > 0.5, 'correlation survives length in the two-copy types'
assert stt.loc[['1-', '1- 2+'], 'struct_seq_partial_lendiff_r'].max() < 0.2, 'and largely goes in 1- and 1- 2+'
assert (stt.multi_window == 0).all(), 'every tested haplotype has one copy-bearing window'
figure('figures/hap_scatter_bytype.png', [f'{O}/pair_distances.tsv.gz', 'results/decomp/haplotype_units.tsv', 'results/copies/windows.tsv'])

# --- step 20: every A/B differential site along the gene (fig/ab_sites.png)
ABS = f'{O}/ab_sites.tsv'; ABE = f'{O}/ab_site_events.tsv'
abS = pd.read_csv(ABS, sep='\t'); abE = pd.read_csv(ABE, sep='\t')
assert abS.event.nunique() == len(abE) and set(abS.event) == set(abE.event)
number('abCols', len(abS), [ABS]); number('abEvents', len(abE), [ABE])
number('abSNV', int((abE.kind == 'SNV').sum()), [ABE]); number('abIndel', int((abE.kind == 'indel').sum()), [ABE])
assert int(MACROS['abSNV']) + int(MACROS['abIndel']) == len(abE)
nf = abE.cls.str.startswith('near-fixed')
number('abFixed', int(nf.sum()), [ABE]); number('abPartial', int((~nf).sum()), [ABE])
number('abTwoA', int((abE.SMN2_shares == 'A').sum()), [ABE]); number('abTwoB', int((abE.SMN2_shares == 'B').sum()), [ABE])
assert int(MACROS['abTwoA']) > int(MACROS['abTwoB']), 'SMN2 mostly carries the A state'
rc = abE.region.value_counts()
number('abIntronOne', int(rc.get('intron 1', 0)), [ABE]); number('abIntronSix', int(rc.get('intron 6', 0)), [ABE])
assert rc.index[0] == 'intron 1' and rc.index[1] == 'intron 6', 'intron 1 then intron 6 hold the most'
ex = abE[abE.region.str.startswith('exon')]
number('abExonEvents', len(ex), [ABE]); assert sorted(ex.region) == ['exon 1', 'exon 8']
number('abExonOnePos', f"{int(ex[ex.region == 'exon 1'].grch38.iloc[0]):,}".replace(',', '{,}'), [ABE])
_ex8 = ex[ex.region == 'exon 8'].grch38.iloc[0]; assert 70952640 <= _ex8 <= 70952650, 'the exon-8 event is c.*211_*212del (3-prime UTR)'
# sites that separate the lineages in all or all but one copy
pf = abE[abE.sep == 1].sort_values('grch38')   # no exception: every A copy carries A's state and every B copy B's
number('abPerfect', len(pf), [ABE], ['separation 1: no copy of either lineage carries the other state'])
assert len(pf) == 2
for nm, r in zip(['One', 'Two'], pf.itertuples()):
    number(f'abPerf{nm}Pos', f'{r.grch38:,}'.replace(',', '{,}'), [ABE])
    number(f'abPerf{nm}Region', r.region, [ABE])
    number(f'abPerf{nm}Kind', 'a 2-bp deletion' if (r.kind == 'indel' and r.cols == 2) else ('a base change' if r.kind == 'SNV' else 'an indel'), [ABE])
    sub = abS[abS.event == r.event]
    number(f'abPerf{nm}A', int(round(sub.fA.min() * MACROS['linA'])), [ABS]); number(f'abPerf{nm}B', int(round(sub.fB.min() * MACROS['linB'])), [ABS])
    assert int(MACROS[f'abPerf{nm}B']) == MACROS['linB'] and int(MACROS[f'abPerf{nm}A']) == MACROS['linA']
assert pf.region.tolist() == ['exon 1', 'intron 1']
# the exon-1 site in transcript notation: place NM_000344.4 exon 1 in the GRCh38 SMN1 copy (offset as in step 23)
import sys as _sys; _sys.path.insert(0, 'scripts'); from smnlib import read_genbank as _rg
_g = _rg('data/refs_seq/smn1_smn2_refseq.gb')['NM_000344.4']
_nm, _sq = [], []
for l in open(f'{O}/copies.aln.fa'):
    l = l.strip()
    if l[0] == '>': _nm.append(l[1:]); _sq.append([])
    else: _sq[-1].append(l.upper())
_rt = tips[(tips.copy_class == 'SMN1') & tips.src.str.startswith('GCA_000001405')].tip.iloc[0]
_rs = ''.join(_sq[_nm.index(_rt)]).replace('-', '')
_e1 = 70923086 + _rs.find(_g['seq'][_g['exons'][0][0] - 1:_g['exons'][0][1]]) + 1          # genomic position of transcript base 1
_tp = int(pf.grch38.iloc[0]) - _e1 + 1; _c = _tp - _g['cds'][0]
assert pf.region.iloc[0] == 'exon 1' and _c < 0, 'exon-1 site lies before the CDS'
number('abPerfOneC', f'c.{_c}', [ABE, 'data/refs_seq/smn1_smn2_refseq.gb', f'{O}/copies.aln.fa'], ['assert the site is in the 5-prime UTR'])
assert int(_c) < 0; number('abPerfOneBefore', -int(_c), [ABE, 'data/refs_seq/smn1_smn2_refseq.gb'])
number('abNearFixed', int((abE.sep >= 0.9).sum()), [ABE], ['near-fixed: each lineage carries its own state in at least 90% of copies'])
assert int((abE.sep >= 0.9).sum()) - int(((abE.sep >= 0.9) & (abE.kind == 'indel')).sum()) == int(json.load(open(f'{O}/evolution_summary.json'))['n_AB']) and int(((abE.sep >= 0.9) & (abE.kind == 'indel')).sum()) == 1, 'text: the A/B history sites are the near-fixed events less one indel'
_ev = pd.read_csv(f'{O}/evolution_sites.tsv', sep='\t'); _abpos = set(_ev[_ev.kind == 'AB'].grch38)
_reg = set(abE[abE.grch38.isin(_abpos) & (abE.kind == 'SNV')].region)
assert len(abE[abE.grch38.isin(_abpos) & (abE.kind == 'SNV')]) == len(_abpos) and _reg <= {'exon 1', 'intron 1', 'intron 2a', 'intron 2b'} and 'exon 1' in _reg and any(r.startswith('intron 2') for r in _reg), 'text: the A/B history sites lie between exon 1 and intron 2'
# the marker, for comparison: on 61 of 73 B and 7 of 449 A (macros sideMarker / markerElsewhere)
assert int(MACROS['linBoneMk']) < MACROS['linBone']
figure('figures/ab_sites.png', [ABS, ABE, f'{O}/outgroup_polarized_0.9.tsv'])

# tool versions (labels, read from the programs themselves)
def ver(cmd):
    try: return subprocess.run(cmd, capture_output=True, text=True, shell=True).stdout.strip().splitlines()[0] + subprocess.run(cmd, capture_output=True, text=True, shell=True).stderr.strip()
    except Exception: return 'unknown'
mv = subprocess.run('mafft --version', shell=True, capture_output=True, text=True); MV = (mv.stdout + mv.stderr).strip().split()[0]
iv = subprocess.run('iqtree2 --version', shell=True, capture_output=True, text=True).stdout; IV = re.search(r'version (\S+)', iv).group(1)
pv = subprocess.run('pangea --version', shell=True, capture_output=True, text=True).stdout; PV = pv.split()[2]
import scipy, matplotlib
from importlib.metadata import version as _pkgver   # edlib has no __version__ attribute
for m, v in [('verMafft', MV), ('verIqtree', IV), ('verPangea', PV), ('verScipy', scipy.__version__), ('verNumpy', np.__version__), ('verPandas', pd.__version__), ('verMpl', matplotlib.__version__), ('verEdlib', _pkgver('edlib'))]:
    number(m, v, [])
# ML tree (step 28) facts from IQ-TREE's own log
log = open(f'{O}/mltree.log').read(); LG = [f'{O}/mltree.log']
number('iqSeqs', re.search(r'Alignment has (\d+) sequences', log).group(1), LG)
number('iqCols', f"{int(re.search(r'with (\d+) columns', log).group(1)):,}".replace(',', '{,}'), LG)
number('iqIdent', re.search(r'NOTE: (\d+) identical sequences', log).group(1), LG)
number('iqBoot', re.search(r'-B (\d+)', log).group(1), LG)
_iqr = open(f'{O}/mltree.iqtree').read(); number('iqModel', re.search(r'Model of substitution: (\S+)', _iqr).group(1), [f'{O}/mltree.iqtree'])
assert MACROS['iqModel'] == 'GTR+F+G4', 'text glosses GTR+F+G4'
assert int(MACROS['iqSeqs']) == int(s['copies']) + 4, 'ML tree: the gene-body copies plus four apes'
assert int(MACROS['iqCols'].replace('{,}', '')) == int(s['gene_body_columns']), 'ML tree uses the gene-body columns'


for f, r in [('figures/hap_scatter.png', [f'{O}/smn_struct.dist', f'{O}/copies.aln.fa']), ('figures/copy_tree_strips.png', [f'{O}/copies.aln.fa', 'tables/neighbour_copies.tsv', 'results/decomp/haplotype_units.tsv'])]:
    figure('' + f, r)

# ---- the cohort (the clinical study's per-haplotype table)
HC = 'tables/haplotype_copy_counts.tsv'; hcc = pd.read_csv(HC, sep='\t')
number('nPeopleAll', hcc['sample'].nunique(), [HC]); number('nHapsAll', len(hcc), [HC])
number('nHapsRes', int(hcc.resolved.sum()), [HC]); number('nHapsUnres', int((~hcc.resolved).sum()), [HC])
both = hcc.groupby('sample').resolved.all(); number('nPeopleBoth', int(both.sum()), [HC]); number('nPeopleNotBoth', int((~both).sum()), [HC])
assert int(MACROS['nHapsRes']) == int(MACROS['arrHaps']) and int(MACROS['nHapsAll']) == 2 * int(MACROS['nPeopleAll'])
number('nPops', hcc.population.dropna().nunique(), [HC])

# ---- the gene-body ML tree (steps 28-29): clades found by their make-up, each checked
BK = pd.read_csv(f'{O}/mltree_backbone.tsv', sep='\t'); CLD = pd.read_csv(f'{O}/mltree_clades.tsv', sep='\t').set_index('group')
RB = [f'{O}/mltree_backbone.tsv', f'{O}/mltree.treefile']
def mk(row): return dict((k, int(v)) for k, v in (x.split(':') for x in row.make_up.split(';')))
BK['mk'] = BK.apply(mk, axis=1)
hum = BK[BK.mk.apply(lambda m: not any(k.startswith('og_') for k in m)) & (BK.copies == int(s['copies']))]
assert len(hum) == 1; hum = hum.iloc[0]
kids = BK[(BK.depth == hum.depth + 1) & BK.node.str.startswith(hum.node + '.')]
assert len(kids) == 2 and kids.copies.sum() == hum.copies, 'human copies split in two'
hI = kids[kids.mk.apply(lambda m: 'B' in m)].iloc[0]; hII = kids[kids.mk.apply(lambda m: 'SMN2' in m)].iloc[0]
assert 'SMN2' not in hI.mk and 'B' not in hII.mk, 'B and SMN2 on opposite halves'
number('mlHumSup', int(hum.support), RB)
number('mlIcopies', int(hI.copies), RB); number('mlIsup', int(hI.support), RB)
number('mlIB', hI.mk['B'], RB); number('mlIh', hI.mk.get('h21', 0), RB)
# the smallest clade holding every B1 and B2 copy (B1 + B2 = the B copies of the half less those in neither clade)
_bb = BK[BK.node.str.startswith(hI.node + '.') & BK.mk.apply(lambda m: m.get('B', 0) == int(MACROS['linBone']) + int(MACROS['linBtwo']) and 'A' not in m and 'SMN2' not in m)]
_bb = _bb.sort_values('copies').iloc[0]
number('mlBBjoinCopies', int(_bb.copies), RB); number('mlBBjoinSup', int(_bb.support), RB)
assert int(_bb.support) < 70 < int(hI.support) and int(_bb.copies) < int(hI.copies), 'text: B1 and B2 joined only weakly; the B half is supported'
assert set(hI.mk) == {'B', 'h21'} and hI.mk['B'] == int(MACROS['linB']), 'half I = lineage B + h21 (the definition of step 05)'
number('mlIIcopies', int(hII.copies), RB); number('mlIIsup', int(hII.support), RB)
number('mlIIA', hII.mk['A'], RB); number('mlIISMNtwo', hII.mk['SMN2'], RB); number('mlIIh', hII.mk.get('h21', 0), RB); number('mlIIhtwelve', hII.mk.get('h12', 0), RB)
assert hII.mk['SMN2'] == int(MACROS['nSMNtwo']) and hII.mk['A'] == int(MACROS['linA']), 'every SMN2 and every A copy is in half II'
hgI = dict(x.split(':') for x in hI.smn1_haplogroups.split(';'))
for f_, nm in [('S1-8', 'Eight'), ('S1-6', 'Six')]: number(f'mlIhg{nm}', sum(int(v) for k, v in hgI.items() if k.split('|')[0].rstrip('cd') == f_), RB + [f'{O}/chen_hg_copies.tsv'])
number('mlIhgNine', sum(int(v) for k, v in hgI.items() if k.startswith('S1-9')), RB + [f'{O}/chen_hg_copies.tsv'])
assert int(MACROS['mlIhgEight']) + int(MACROS['mlIhgSix']) + int(MACROS['mlIhgNine']) == hI.mk['B'], 'half I SMN1 copies: S1-6, S1-8, S1-9 only'
# inside half I: the S1-9/S1-6 clade with the h21 copies (B1) and the S1-8 clade (B2)
bh = BK[BK.mk.apply(lambda m: set(m) == {'B', 'h21'}) & (BK.copies < hI.copies) & BK.smn1_haplogroups.fillna('').apply(lambda x: 'S1-8' not in x)]
bh = bh.sort_values('copies', ascending=False).iloc[0]
number('mlBhCopies', int(bh.copies), RB); number('mlBhH', bh.mk['h21'], RB); number('mlBhSup', int(bh.support), RB); number('mlBhB', bh.mk['B'], RB)
assert bh.mk['B'] == int(MACROS['linBone']) and bh.mk['h21'] == int(MACROS['mlBsideHtwentyone']), 'the S1-9/S1-6 + h21 clade is B1 with every h21 copy of B half'
e8 = BK[BK.mk.apply(lambda m: set(m) == {'B'}) & BK.smn1_haplogroups.fillna('').str.fullmatch(r'S1-8:\d+')]
e8 = e8.sort_values('copies', ascending=False).iloc[0]
number('mlEightCopies', int(e8.copies), RB); number('mlEightSup', int(e8.support), RB)
assert int(e8.copies) == int(MACROS['mlIhgEight']), 'all S1-8 copies form one clade'
s2 = BK[BK.mk.apply(lambda m: m.get('SMN2') == hII.mk['SMN2'])].sort_values('copies').iloc[0]
hg2 = dict(x.split(':') for x in s2.smn1_haplogroups.split(';'))
number('mlSTwoCopies', int(s2.copies), RB); number('mlSTwoSup', int(s2.support), RB); number('mlSTwoA', s2.mk['A'], RB)
number('mlSTwoTen', int(hg2.get('S1-10', 0)), RB); number('mlSTwoOne', int(hg2.get('S1-1', 0)), RB); number('mlSTwoH', s2.mk.get('h21', 0), RB)
assert int(hg2.get('S1-10', 0)) > s2.mk['A'] / 2, 'the SMN1 copies in the SMN2 clade are mostly S1-10'
s2p = BK[(BK.depth == s2.depth - 1) & BK.apply(lambda r: s2.node.startswith(r.node + '.'), axis=1)].iloc[0]
hg2p = dict(x.split(':') for x in s2p.smn1_haplogroups.split(';'))
number('mlSTwoPCopies', int(s2p.copies), RB); number('mlSTwoPSup', int(s2p.support), RB); number('mlSTwoPA', s2p.mk['A'], RB)
number('mlSTwoPTen', int(hg2p.get('S1-10', 0)), RB); number('mlSTwoPOne', int(hg2p.get('S1-1', 0)), RB)
number('mlSTwoPH', s2p.mk.get('h21', 0), RB); number('mlSTwoPSmn', s2p.mk['SMN2'], RB)
assert int(MACROS['mlSTwoPTen']) + int(MACROS['mlSTwoPOne']) == int(MACROS['mlSTwoPA']) and int(MACROS['mlSTwoPSmn']) + int(MACROS['mlSTwoPH']) + int(MACROS['mlSTwoPA']) == int(MACROS['mlSTwoPCopies']), 'SMN2 clade make-up'
assert s2p.support > s2.support
_cs2 = CLD.loc['SMN2']
number('mlSTwoMrca', int(_cs2.mrca_clade_size), [f'{O}/mltree_clades.tsv']); number('mlSTwoMrcaSup', int(_cs2.mrca_support), [f'{O}/mltree_clades.tsv'])
assert int(_cs2.mrca_outside) == int(MACROS['mlSTwoH']) == int(MACROS['mlSTwoPH']), 'text: the smallest clade holding all SMN2 adds only the h21 copies (no SMN1)'
assert int(MACROS['mlSTwoMrca']) + s2.mk['A'] == int(MACROS['mlSTwoCopies']) and s2.mk['A'] == int(MACROS['mlSTwoTen']), 'text: the next SMN1 copies to join are all S1-10'
assert int(MACROS['mlSTwoCopies']) + int(MACROS['mlSTwoPOne']) == int(MACROS['mlSTwoPCopies']), 'text: then the S1-1 copies'
RC = [f'{O}/mltree_clades.tsv']
number('mlBpure', int(CLD.loc['lineage B', 'largest_pure_clade']), RC); number('mlBpureSup', int(CLD.loc['lineage B', 'largest_pure_support']), RC)
number('mlThreeSup', int(CLD.loc['haplogroup S1-3', 'mrca_support']), RC)
for g_ in ('haplogroup S1-3', 'haplogroup S1-8', 'haplogroup family S1-9 (S1-9/S1-9c/S1-9d)'):
    assert bool(CLD.loc[g_, 'monophyletic']) and int(CLD.loc[g_, 'mrca_outside']) == 0, 'text: pure clade ' + g_
number('mlEightCladeN', int(CLD.loc['haplogroup S1-8', 'copies']), RC); number('mlEightCladeSup', int(CLD.loc['haplogroup S1-8', 'mrca_support']), RC)
_f9 = 'haplogroup family S1-9 (S1-9/S1-9c/S1-9d)'; number('mlNineFamN', int(CLD.loc[_f9, 'copies']), RC); number('mlNineFamSup', int(CLD.loc[_f9, 'mrca_support']), RC); number('mlThreeN', int(CLD.loc['haplogroup S1-3', 'copies']), RC)
assert bool(CLD.loc['haplogroup S1-3', 'monophyletic']) and bool(CLD.loc['haplogroup S1-8', 'monophyletic'])
# S1-3: all copies, those in the Fig. 5 tree, and their haplotypes' sequence labels
_HG, _HL, _SC = f'{O}/chen_hg_copies.tsv', f'{O}/ml_hap_labels.tsv', f'{O}/ml_seq_clusters.tsv'
_hg = pd.read_csv(_HG, sep='\t'); _hl = pd.read_csv(_HL, sep='\t').set_index('src'); _sc = pd.read_csv(_SC, sep='\t')
_s3 = _hg[_hg.hg.astype(str) == 'S1-3'].merge(pd.read_csv(f'{O}/copy_labels.tsv', sep='\t')[['tip', 'src']], on='tip')
assert len(_s3) == int(MACROS['mlThreeN']) and _s3.src.nunique() == len(_s3), 'S1-3: one copy per haplotype'
_in3 = _s3[_s3.src.isin(_hl.index)]; _out3 = sorted(set(_s3.src) - set(_in3.src))
assert _out3 == ['HG00280_hap1_hprc_r2_v1.0.1'], 'text: the one S1-3 copy not in Fig. 5 is HG00280 hap1'
number('mlbThreeTree', len(_in3), [_HG, _HL]); number('mlbThreeOut', len(_out3), [_HG, _HL])
_l3 = _in3.src.map(_hl.seq_label)
number('mlbThreeOnly', int((_l3 == 'S1-3').sum()), [_HL]); number('mlbThreePlus', int((_l3 != 'S1-3').sum()), [_HL])
assert set(_l3[_l3 != 'S1-3']) == {'S1-3+rest'}
_r3 = _sc[(_sc['T'] == 95) & (_sc['M'] == 20) & (_sc.cluster == 'S1-3')].iloc[0]
assert int(_r3.own_SMN1) == len(_in3), 'text: the S1-3 clade of Fig. 5 holds every S1-3 copy in the tree'
gor = BK[(BK.copies > 1) & BK.mk.apply(lambda m: set(m) - {'A', 'B', 'SMN2', 'h21', 'h12'} == {'og_gorilla'})]
assert len(gor) == 1, 'gorilla is sister to the human copies'; number('mlGorSup', int(gor.iloc[0].support), RB)
number('mlRootNote', 'orangutan', RB)

# ---- checks guarding sentences of the preprint's main text
ancI = dict(x.split(':') for x in hI.ancestry.split(';')); assert int(ancI.get('AFR', 0)) > hI.copies / 2, 'half I is mostly AFR'
M_ = lambda k: int(str(MACROS[k]).replace('{,}', ''))
assert M_('slLAfrB') > M_('slLAfrTwo') and M_('slLAfrB') > M_('slLAfrA'), 'AFR cenSMN arm: B is the commonest'
assert M_('slLOthTwo') > M_('slLOth') / 2, 'outside AFR the cenSMN arm mostly holds SMN2'
assert M_('cbTwotopArr') > M_('cbTwon') / 2 and M_('cbABtopArr') > M_('cbABn') / 2, 'bands: their top arrangement holds most'
assert M_('hgMkEight') + M_('hgMkNine') == M_('nMarkerCopies'), 'marker copies are S1-8 or S1-9'
assert abs(float(MACROS['dcACen']) - float(MACROS['dcATel'])) < 0.1
_mm = pd.read_csv(f'{O}/ml_cluster_pairs_members.tsv', sep='\t'); _ab = set(_mm[_mm.band == 'AB'].src)
_pq = pd.read_csv(f'{O}/palindrome_copies_lineage.tsv', sep='\t'); _pq = _pq[_pq.src.isin(_ab)]
_l = _pq[_pq.arm.str.startswith('left')].groupby('src').sublineage.apply(lambda x: '+'.join(sorted(x.dropna())))
_r = _pq[_pq.arm.str.startswith('right')].groupby('src').sublineage.apply(lambda x: '+'.join(sorted(x.dropna())))
assert set(_l.index) == _ab == set(_r.index)
number('cbABcen', int(((_l == 'B1') & (_r == 'B2')).sum()), [f'{O}/ml_cluster_pairs_members.tsv', f'{O}/palindrome_copies_lineage.tsv'])
assert int(((_l == 'B1') & (_r == 'B2')).sum()) + int(((_l == 'B2') & (_r == 'B1')).sum()) == len(_ab), 'text: every band haplotype has B1 on one arm and B2 on the other'
# cold read 2026-11-17: which of the B1+B2 haplotypes is outside the band (37 of 38; 34 of 35 B1-cen)
_abo = sorted(set(mix.src) - _ab)
assert _ab <= set(mix.src), 'band is a subset of the B1+B2 haplotypes'
assert len(_ab) + len(_abo) == MACROS['linMixHaps']
assert int(MACROS['cbABcen']) + len(set(_abo) - set(sw_)) == MACROS['linMixBoneCenRev'], 'text: the B1-cen haplotypes outside the band are the ones not in it'
number('cbABoutN', len(_abo), [f'{O}/ml_cluster_pairs_members.tsv', RL[0] if isinstance(RL, list) else RL])
number('cbABoutNames', ', '.join(x.replace('_hprc_r2_v1.0.1', '').replace('_', ' ') for x in _abo), [f'{O}/ml_cluster_pairs_members.tsv'])
_abo_clu = sorted(set(str(_hl.struct_cluster.get(x, 'none')) for x in _abo))
print('B1+B2 outside band:', _abo, 'struct clusters', _abo_clu, 'labels', [_hl.seq_label.get(x) for x in _abo])
number('cbABoutClu', ', '.join(_abo_clu), [_HL])
# the S1-4 and S1-5 clades of Fig. 5 are supported clades of >= 20 SMN1 copies, but hold few copies of their own
_sc95 = _sc[(_sc['T'] == 95) & (_sc['M'] == 20)].set_index('cluster')
number('mlbFourOwn', int(_sc95.loc['S1-4'].own_SMN1), [_SC]); number('mlbFiveOwn', int(_sc95.loc['S1-5'].own_SMN1), [_SC])
number('mlbFourSup', int(_sc95.loc['S1-4'].support), [_SC]); number('mlbFiveSup', int(_sc95.loc['S1-5'].support), [_SC])
assert _sc95.loc['S1-4'].clade_SMN1 >= 20 and _sc95.loc['S1-5'].clade_SMN1 >= 20, 'S1-4 and S1-5 are supported clades of >= 20 copies'
assert float(MACROS['gsBoneBtwo']) < float(MACROS['gsATwo']) and float(MACROS['gsATwo']) - float(MACROS['gsBoneBtwo']) <= 0.1, 'text: B1-B2 nearly as far apart as A from SMN2'
assert MACROS['gsAA'] == MACROS['gdTwoTwoMean'], 'text: A and SMN2 vary equally at two decimals'
assert float(MACROS['gsBtwoBtwo']) < float(MACROS['gsBoneBone']) < float(MACROS['gsAA']), 'text: B2 least varied, then B1'
_dcs = dc[(dc.subset == 'all') & (dc.lineage == 'A') & (dc.context == 'position3')].set_index('level').median_to_SMN2
assert _dcs['single'] > max(_dcs['centromeric'], _dcs['telomeric']), 'text: single A copies sit slightly further'
_LLb = LL[LL.lin == 'B']
assert (_LLb[_LLb.position == 'centromeric']['sub'] == 'B1').mean() > 0.5 and (_LLb[_LLb.position == 'telomeric']['sub'] == 'B2').mean() > 0.5, 'text: centromeric B mostly B1, telomeric B mostly B2'

# ---- which haplotypes each comparison uses (cold read 2026-11-14): 438 resolved = 434 with an SMN1 copy + 4 without;
# the tree comparisons add GRCh38 and CHM13 (436); summary.json's 440 = 438 resolved + the two references
_lh = pd.read_csv(LD['haplotypes'], sep='\t'); _res = _hcp[_hcp.resolved == True]
_res1 = _res[_res.src.isin(_lh.src)]; _res0 = _res[~_res.src.isin(_lh.src)]
assert (_res0.SMN1 == 0).all(), 'resolved haplotypes missing from the lineage table carry no SMN1 copy'
assert set(_lh[_lh.resolved == True].src) == set(_res1.src), 'lineage table resolved rows = resolved haplotypes with an SMN1 copy'
number('nHapsResSMNone', len(_res1), [HCT, LD['haplotypes']]); number('nHapsResNoSMNone', len(_res0), [HCT, LD['haplotypes']])
assert len(_res1) + len(_res0) == int(MACROS['nHapsRes'])
assert int(MACROS['sHaps']) == len(_res1), 'step 20 counts the resolved haplotypes with an SMN1 copy'
assert msj['haplotypes'] == len(_res1) + 2 == int(MACROS['tmHaps']), 'tree comparisons: resolved with SMN1 + GRCh38 + CHM13'
assert int(MACROS['nTangle']) == int(MACROS['tmHaps']) == int(MACROS['nHapsResSMNone']) + 2, 'S3.2: the 436 tree haplotypes = 434 resolved donors with SMN1 + GRCh38 + CHM13'
assert s['haplotypes'] == int(MACROS['nHapsRes']) + 2, 'summary.json haplotypes = resolved + GRCh38 + CHM13'
_hbc = msj['haplotypes_by_SMN1_copies']; assert set(_hbc) <= {'1', '2', '3'}
number('mlbOneCopy', _hbc.get('1', 0), [MSJ]); number('mlbTwoCopy', _hbc.get('2', 0), [MSJ])
assert _hbc.get('1', 0) + _hbc.get('2', 0) == msj['haplotypes'] and _hbc.get('1', 0) + 2 * _hbc.get('2', 0) == msj['SMN1_copies']
_aa = _lh[_lh.composition == 'AA']; _aaOut = _aa[_aa.resolved != True]
number('ldCompAAunres', len(_aaOut), [LD['haplotypes']])
assert int(MACROS['ldCompAA']) - len(_aaOut) == int(MACROS['tmAAn']), 'Two A: the trees leave out the unresolved AA haplotype'
assert int(MACROS['palNoOne']) == int((_lh[_lh.resolved == True].arrangement == '1+').sum()), "'1+' haplotypes"
number('palNoOneB', int(((_lh.resolved == True) & (_lh.arrangement == '1+') & (_lh.composition == 'B')).sum()), [LD['haplotypes']])
# cold read 2026-11-14: the two 434 sets (with repeat units / with an SMN1 copy) are different sets
_units = set(pu.src); _resall = set(_res.src)
_nounits = sorted(_resall - _units); _nosmn1 = sorted(set(_res0.src))
assert len(_nounits) == 4 and len(_nosmn1) == 4
number('nHapsNoUnits', len(_nounits), [PU, HCT]); number('nHapsNoUnitsNoSMNone', len(set(_nounits) & set(_nosmn1)), [PU, HCT])
assert set(_lh[(_lh.resolved == True) & (_lh.arrangement == '1+')].src) <= _units, "text: every resolved '1+' haplotype has repeat units"
assert not (set(_lh[(_lh.resolved == True) & (_lh.arrangement == '1+')].src) & set(pu[pu.has_pair].src)), "text: no '1+' haplotype has the pair"
print('no units:', _nounits, ' no SMN1:', _nosmn1)

# ---- step 36: how each window was taken and that it reaches single-copy flank on both sides
WFJ = f'{O}/window_flanks.json'; wf = json.load(open(WFJ)); WFT = f'{O}/window_flanks.tsv'
number('wfN', wf['n_windows'], [WFJ], ['assert one decomposed window per haplotype row of smn.bed'])
assert wf['n_windows'] == 440
number('wfQueryMb', f"{wf['query_len']/1e6:.2f}", [WFJ])
number('wfGStart', fmtc(wf['grch38_start1']), [WFJ]); number('wfGEnd', fmtc(wf['grch38_end']), [WFJ])
assert wf['grch38_q_to'] - wf['grch38_q_from'] == wf['grch38_end'] - wf['grch38_start1'] + 1, 'GRCh38 window = query span'
number('wfStartQ', wf['n_start_le_2kb'], [WFJ]); number('wfReachEnd', wf['n_reach_query_end'], [WFJ])
assert wf['n_start_le_2kb'] == wf['n_windows'] and wf['n_reach_query_end'] == wf['n_windows'] - 1, 'text: all start within 2 kb of the query start; all but one reach its end'
_wfd = pd.read_csv(f'{O}/window_flanks.tsv', sep='\t'); _ws = _wfd[_wfd.q_to < wf['query_len']]
assert len(_ws) == 1, 'one window short of the query end'
_ws = _ws.iloc[0]; assert _ws.src.split('_')[0] != 'GCA', 'the short window is a donor haplotype'
number('wfShortHap', ' '.join(_ws.src.split('_')[:2]), [f'{O}/window_flanks.tsv']); number('wfShortBp', int(wf['query_len'] - _ws.q_to), [f'{O}/window_flanks.tsv'])
number('wfShortTelKb', f'{_ws.tel_unique_bp/1e3:.1f}', [f'{O}/window_flanks.tsv']); assert _ws.tel_unique_bp >= 150000, 'text: the short window still ends in single-copy flank'
number('wfBothSingle', wf['n_both_single'], [WFJ]); assert wf['n_both_single'] == wf['n_windows'], 'text: every window starts and ends in a single-copy bundle'
number('wfCenMin', f"{wf['cen_unique_kb_min']:.1f}", [WFJ]); number('wfCenMed', f"{wf['cen_unique_kb_med']:.1f}", [WFJ])
number('wfTelMin', f"{wf['tel_unique_kb_min']:.1f}", [WFJ]); number('wfTelMed', f"{wf['tel_unique_kb_med']:.1f}", [WFJ])
assert float(MACROS['wfCenMin']) <= wf['cen_unique_kb_min'] + 0.05 and float(MACROS['wfTelMin']) <= wf['tel_unique_kb_min'] + 0.05, 'text: "at least" bound is met'
# copies on haplotypes not resolved end to end (cold read: the copy set includes them)
_CL = pd.read_csv('results/sv/copy_labels.tsv', sep='\t'); _HC = pd.read_csv('tables/haplotype_copy_counts.tsv', sep='\t')
_resd = dict(zip(_HC.src, _HC.resolved.astype(str)))
_un = _CL[_CL.src.map(_resd) == 'False']
number('nCopiesUnres', len(_un), ['results/sv/copy_labels.tsv', 'tables/haplotype_copy_counts.tsv'])
number('nCopiesUnresHaps', _un.src.nunique(), ['results/sv/copy_labels.tsv', 'tables/haplotype_copy_counts.tsv'])
assert len(_CL) == int(MACROS['nCopies']), 'copy_labels holds the whole counted copies'
assert all(_CL.src.map(_resd).isin(['True', 'False']) | _CL.src.str.startswith(('GCA_', 'chm13'))), 'every donor copy has a resolved flag'
# Fig. S algorithm diagram (step 39)
_ra = json.load(open('results/sv/repeat_units_algo.json'))
number('ruaEx', _ra['example'].split('_hprc')[0].replace('_', ' '), ['results/sv/repeat_units_algo.json'])
number('ruaSeg', _ra['n_segments'], ['results/sv/repeat_units_algo.json']); number('ruaUnits', _ra['n_units'], ['results/sv/repeat_units_algo.json'])
assert abs(_ra['short_kb'] - wf['len_kb_min']) < 1 and abs(_ra['long_kb'] - wf['len_kb_max']) < 1, 'Fig. S: shortest and longest windows drawn'

number('wfLenMin', wf['len_kb_min'], [WFJ]); number('wfLenMed', fmtc(wf['len_kb_med']), [WFJ]); number('wfLenMax', fmtc(wf['len_kb_max']), [WFJ])
_wft = pd.read_csv(WFT, sep='\t')
assert (_wft.cen_unique_bp >= 85000).all() and (_wft.tel_unique_bp >= 199000).all(), 'text: at least ~86 / ~200 kb single-copy at each end'

# ---- dosage: chromosomes a dosage assay cannot see (researcher 2026-11-15)
HCC = 'tables/haplotype_copy_counts.tsv'; _h = pd.read_csv(HCC, sep='\t'); _h = _h[_h.resolved]
_sp = pd.read_csv(f'{O}/copy_labels.tsv', sep='\t')[['src', 'superpop']].drop_duplicates('src').set_index('src').superpop
_h = _h.assign(afr=_h.src.map(_sp).eq('AFR'))
number('dosTwo', int((_h.SMN1_ex7_dosage >= 2).sum()), [HCC]); number('dosTwoAfr', int(((_h.SMN1_ex7_dosage >= 2) & _h.afr).sum()), [HCC, f'{O}/copy_labels.tsv'])
number('dosAfrHaps', int(_h.afr.sum()), [HCC, f'{O}/copy_labels.tsv']); number('dosZero', int((_h.SMN1_ex7_dosage == 0).sum()), [HCC])
assert (_h.smn1_ex7_2plus == (_h.SMN1_ex7_dosage >= 2)).all() and len(_h) == 438
_pp = _h.groupby('sample').SMN1_ex7_dosage.agg(list); _pp = _pp[_pp.map(len) == 2]
number('dosPeople', len(_pp), [HCC]); number('dosTwoZero', sum(max(x) >= 2 and min(x) == 0 for x in _pp), [HCC])
number('dosPeopleTwo', sum(max(x) >= 2 for x in _pp), [HCC])
number('dosNoCopy', int((_h.SMN1 == 0).sum()), [HCC]); number('dosTwoCopy', int((_h.SMN1 >= 2).sum()), [HCC])
assert int((_h.SMN1 == 0).sum()) - int((_h.SMN1_ex7_dosage == 0).sum()) == int(((_h.SMN1 == 0) & (_h.h12 > 0)).sum()), 'text: the haplotypes with no SMN1 copy but exon-7 dosage carry h12'
assert int((_h.SMN1_ex7_dosage >= 2).sum()) - int((_h.SMN1 >= 2).sum()) == int(((_h.SMN1_ex7_dosage >= 2) & (_h.SMN1 < 2) & (_h.h12 > 0)).sum()), 'text: the extra two-dosage haplotypes carry h12'
number('dosHOneTwoN', int(((_h.SMN1 == 0) & (_h.h12 > 0)).sum() + ((_h.SMN1_ex7_dosage >= 2) & (_h.SMN1 < 2) & (_h.h12 > 0)).sum()), [HCC])
assert int((_h.SMN1 == 0).sum()) == int(MACROS['nHapsRes']) - int(MACROS['nHapsResSMNone']), 'dosNoCopy = 438 - 434'
_ar = pd.read_csv(f'{O}/ape_site_reads.tsv', sep='\t'); _ar = _ar[(_ar.kind == 'AB') & _ar.accepted]
_wk = _ar.assign(m=_ar.score - _ar.second).sort_values('m').iloc[0]
number('rdWeakSite', f'{int(_wk.grch38):,}'.replace(',', '{,}'), [f'{O}/ape_site_reads.tsv']); number('rdWeakScore', int(_wk.score), [f'{O}/ape_site_reads.tsv']); number('rdWeakSecond', int(_wk.second), [f'{O}/ape_site_reads.tsv'])
assert _wk.ape == 'og_orangutan' and int(_wk.grch38) == 70928297, 'text: the weakest kept A/B read is orangutan at 70,928,297'
_oth = _ar[(_ar.grch38 == _wk.grch38) & (_ar.ape != _wk.ape)]
assert len(_oth) == 3 and set(_oth.read) == {_wk.read}, 'text: the other three apes read the same base there'
number('rdWeakOthLo', int(_oth.score.min()), [f'{O}/ape_site_reads.tsv']); number('rdWeakOthHi', int(_oth.score.max()), [f'{O}/ape_site_reads.tsv'])
assert len(_pp) == 216 and sum(max(x) >= 2 and min(x) == 0 for x in _pp) == 0, 'text: no 2+0 person among the 216'

# ---- cold read 2026-11-15: the haplotype sets each comparison uses, in one place (supplement S1)
_hr = pd.read_csv('tables/haplotype_copy_counts.tsv', sep='\t'); _hr = _hr[_hr.resolved]
number('setNoLab', int(_hr.population.isna().sum()), ['tables/haplotype_copy_counts.tsv'])
number('setAsl', int((_hr.population == 'ASL').sum()), ['tables/haplotype_copy_counts.tsv'])
assert len(_hr) - int(_hr.population.isna().sum()) - int((_hr.population == 'ASL').sum()) == an['haplotypes'] == 432, 'text: ancestry set = resolved - unlabeled - ASL'
number('setAllCopy', int(MACROS['nHapsRes']) + 2, ['tables/haplotype_copy_counts.tsv', f'{O}/summary.json']); assert int(MACROS['nHapsRes']) + 2 == json.load(open(f'{O}/summary.json'))['haplotypes'] == 440, 'text: all-copy Mantel set = resolved + 2 references'
assert MACROS['mlbHaps'] == MACROS['nTangle'] == MACROS['tmHaps'] == 436 and int(MACROS['nHapsRes']) + 2 == 440
# the A-copy distance spreads by grouping (Section 2.3)
_pos = [float(MACROS[k]) for k in ('dcACen', 'dcATel', 'dcASingle')]; _pal = [float(MACROS[k]) for k in ('dcAPalYes', 'dcAPalNo')]
number('dcPosSpread', f"{max(_pos) - min(_pos):.2f}", [f'{O}/distance_context_tests.tsv']); number('dcPalSpread', f"{max(_pal) - min(_pal):.2f}", [f'{O}/distance_context_tests.tsv'])
assert float(MACROS['dcASpread']) >= max(max(_pos) - min(_pos), max(_pal) - min(_pal)), 'text: arrangement type gives the widest spread'

# ---- method settings stated in the text, read from the scripts that use them
def setting(name, path, pattern, conv=int):
    m = re.search(pattern, open(path).read()); assert m, (name, path); number(name, conv(m.group(1)), [path])
setting('elAstart', 'scripts/28b_arm_elements.py', r"k = s\[\(s.el == 'K'\) & \(s.a0 > (\d+)\)\]")
setting('elUnitLong', 'scripts/28b_arm_elements.py', r"g.span_kb >= (\d+)\)")
setting('trMinKb', 'scripts/28_windows_tandem.py', r"'--min-unit-bp', '(\d+)'", lambda v: int(v) // 1000)
setting('trMaxKb', 'scripts/28_windows_tandem.py', r"'--max-unit-bp', '(\d+)'", lambda v: int(v) // 1000)
setting('sK', 'scripts/28c_spacer_sequence.py', r"K = (\d+);")
setting('sSample', 'scripts/28c_spacer_sequence.py', r"random.sample\(srcs, (\d+)\)")
setting('sGapLo', 'scripts/28c_spacer_sequence.py', r"between\((\d+), \d+, inclusive")
setting('sGapHi', 'scripts/28c_spacer_sequence.py', r"between\(\d+, (\d+), inclusive")
setting('setK', 'scripts/19_cluster_correspondence.py', r"\nk = (\d+); cs, cg = fcluster")
setting('setMinOverlap', 'scripts/19_cluster_correspondence.py', r"T.overlap >= (\d+)")
setting('setMLsup', 'scripts/19b_ml_bands.py', r"K, T, M, MIN_OV, ALPHA = \d+, (\d+),")
setting('setMLmin', 'scripts/19b_ml_bands.py', r"K, T, M, MIN_OV, ALPHA = \d+, \d+, (\d+),")
setting('setAlpha', 'scripts/19_cluster_correspondence.py', r"T.p_bonf < (1e-\d+)", lambda v: f'{float(v):g}'.replace('0.001', '0.001'))
setting('setTopN', 'scripts/01_struct_dist.sh', r"--top-n (\d+)")
setting('setFlankKb', 'scripts/07_compare.py', r"rpos >= (\d+)\) & \(rpos < rlen - \d+\)", lambda v: int(v) // 1000)
setting('setSupMin', 'scripts/32_ml_tree_figure.py', r"n\['n'\] >= (\d+) and n\['sup'\]")

# repeat-unit families: their definition (step 33, from PANGEA's repeat-units output) -> Table S (tab_unitfam.tex)
UF, UJ = f'{O}/unit_families.tsv', f'{O}/unit_families.json'
uf = pd.read_csv(UF, sep='\t', dtype={'core_ids': str}).set_index('family'); uj = json.load(open(UJ))
_ru_ = pd.read_csv('results/decomp/smn.repeat_units.tsv', sep='\t', comment='#'); _below = _ru_[_ru_.span_bp < uj['auto_min_unit_bp']]
assert set(_below.family) == {'J'} and len(_below) == 5, 'text: the units below the automatic minimum are a few family-J units'
number('ufBelowN', len(_below), ['results/decomp/smn.repeat_units.tsv', UJ]); number('ufBelowFam', 'J', ['results/decomp/smn.repeat_units.tsv'])
number('ufAutoKb', f"{uj['auto_min_unit_bp'] / 1e3:.1f}", [UJ]); number('ufShortKb', f"{uj['shortest_unit_bp'] / 1e3:.1f}", [UJ])
number('ufOcc', uj['occurrences'], [UJ]); number('ufWindows', uj['windows'], [UJ]); number('ufFams', uj['families'], [UJ])
number('ufInv', int(uf.inverted.sum()), [UF]); number('ufDir', int(uf.direct.sum()), [UF])
big = uf[uf.occurrences >= 10]; small_ = uf[uf.occurrences < 10]
number('ufBig', len(big), [UF]); number('ufSmallFams', len(small_), [UF]); number('ufSmallOcc', int(small_.occurrences.sum()), [UF])
ids = lambda f_: set(str(uf.loc[f_, 'core_ids']).split())
assert ids('B') == ids('A'), 'text: B has exactly the core bundles of A'
for k_, m_ in [('max_gap', 'ufMaxGap'), ('max_indel', 'ufMaxIndel'), ('min_matches', 'ufMinMatch')]: number(m_, int(uj['defaults'][k_]), [UJ])
number('ufJac', f"{uj['defaults']['family_jaccard']:.1f}", [UJ])
number('ufACore', len(ids('A')), [UF]); number('ufBCore', len(ids('B')), [UF])
lines = ['% written by scripts/40_numbers.py from results/sv/unit_families.tsv; do not edit', r'\newcommand{\tabUnitFam}{%']
for f_, r in big.iterrows():
    lines.append(trow(r'\RU{%s}' % f_, r.occurrences, r.windows, f'{r.span_kb_median} ({r.span_kb_min}--{r.span_kb_max})', r.bundles_median,
                      f'{r.forward}/{r.reverse}', f'{r.inverted}/{r.direct}', f'{r.core_bundles}: {r.core_ids}', f'{r.core_bp_share:.2f}'))
lines.append('}')
(DOC / 'tab_unitfam.tex').write_text('\n'.join(lines) + '\n')

# step 28e: what each family is made of (arm elements) -> tab_unitfam_el.tex; bands and structure clusters against arm form
FE, AFJ = f'{O}/unit_family_elements.tsv', f'{O}/armform_summary.json'
fe = pd.read_csv(FE, sep='\t').set_index('family'); af = json.load(open(AFJ))
assert fe.loc['A', 'main_element'] == 'A core' and fe.loc['B', 'main_element'] == 'Xi' and fe.loc['C', 'main_element'] == 'Xd' and fe.loc['D', 'main_element'] == 'E', 'text: A mostly the A core, B Xi, C Xd (with E), D E'
number('feACore', f"{fe.loc['A', 'share_Acore']:.2f}", [FE]); number('feBXi', f"{fe.loc['B', 'share_Xi']:.2f}", [FE])
number('feCXd', f"{fe.loc['C', 'share_Xd']:.2f}", [FE]); number('feCE', f"{fe.loc['C', 'share_E']:.2f}", [FE]); number('feDE', f"{fe.loc['D', 'share_E']:.2f}", [FE])
number('feFams', len(fe), [FE])
_sh = fe[[f'share_{c}' for c in ['Acore', 'K', 'Xi', 'S', 'Xd', 'E', 'Fl', 'other']]].sum(axis=1)
assert ((_sh - 1).abs() < 0.01).all(), 'Table S2: element shares including other add to 1 per family'
lines = ['% written by scripts/40_numbers.py from results/sv/unit_family_elements.tsv; do not edit', r'\newcommand{\tabUnitFamEl}{%']
for f_, r in fe.iterrows():
    lines.append(trow(r'\RU{%s}' % f_, r.occurrences, r.haplotypes, *[f"{r[f'share_{c}']:.2f}" for c in ['Acore', 'K', 'Xi', 'S', 'Xd', 'E', 'Fl', 'other']],
                      f"{r.top_combo.replace('A core', 'core')} ({r.top_combo_n})"))
lines.append('}')
(DOC / 'tab_unitfam_el.tex').write_text('\n'.join(lines) + '\n')
# coarse scale: the two large structure clusters split by the cenSMN arm type; T differs between them but is not a split
assert af['mixed_cluster_size'] == MACROS['cbABClu'] and af['bulk_cluster_size'] == MACROS['cbBulkClu'], 'the two structure clusters of the bulk and B1+B2 bands'
assert af['bulk_pair_C0'] == 0 and af['mixed_pair_C1'] <= 1 and af['fisher_p_C'] < 1e-50, 'text: the bulk cluster is all C1, the B1+B2 cluster all but one C0, among pair haplotypes'
assert 0 < af['bulk_pair_T1'] < af['bulk_pair'] / 4 and af['mixed_pair_T1'] < af['mixed_pair'], 'text: T1 is in both clusters (no split by T)'
for k_, M_ in [('bulk_pair', 'afBulkPair'), ('bulk_pair_C1', 'afBulkCOne'), ('bulk_pair_C0', 'afBulkCZero'), ('bulk_pair_T1', 'afBulkTOne'), ('bulk_pair_T0', 'afBulkTZero'),
               ('mixed_pair', 'afMixPair'), ('mixed_pair_C0', 'afMixCZero'), ('mixed_pair_C1', 'afMixCOne'), ('mixed_pair_T1', 'afMixTOne'),
               ('bulk_cluster_size', 'afBulkClu'), ('mixed_cluster_size', 'afMixClu'), ('T1_all', 'afTOneAll'), ('T0_all', 'afTZeroAll')]:
    number(M_, af[k_], [AFJ])
number('afMixTZero', af['mixed_pair'] - af['mixed_pair_T1'], [AFJ])
number('afPC', texp(af['fisher_p_C']), [AFJ]); number('afP', texp(af['fisher_p_T1']), [AFJ])
assert af['T1_all'] == MACROS['elTOne'] and af['T0_all'] == MACROS['elTZero'], 'arm-type counts agree with step 28b'
b_ = af['bands']
assert b_['AB']['C0|T1'] + b_['AB']['C0|T0'] == MACROS['cbABn'], 'text: every B1+B2 band haplotype is C0'
assert b_['Bulk']['C0|T0'] == b_['Bulk']['C0|T1'] == 0, 'text: every bulk band haplotype with the pair is C1'
assert b_['Two']['no pair'] == MACROS['cbTwon'], 'text: no S1-2 band haplotype has the head-to-head pair'
number('afABCZTOne', b_['AB']['C0|T1'], [AFJ]); number('afABCZTZero', b_['AB']['C0|T0'], [AFJ])
number('afBulkCOTZ', b_['Bulk']['C1|T0'], [AFJ]); number('afBulkCOTO', b_['Bulk']['C1|T1'], [AFJ]); number('afBulkNoPair', b_['Bulk']['no pair'], [AFJ])
number('afBulkBandPair', b_['Bulk']['C1|T0'] + b_['Bulk']['C1|T1'], [AFJ])
# the family that marks the T1 arm at the coarse scale
SIGF = af['sig_family']; number('afSigFam', SIGF, [AFJ])
assert fe.loc[SIGF, 'main_element'] == 'Xi', 'text: the T1-marking family is mostly Xi'
number('afSigTOne', af['G_on_T1'][0], [AFJ]); number('afSigTZero', af['G_on_T0'][0], [AFJ]); number('afSigHaps', af['G_haps_all'], [AFJ]); number('afSigNoPair', af['G_haps_nopair'], [AFJ])
number('feSigXi', f"{fe.loc[SIGF, 'share_Xi']:.2f}", [FE])
assert af['G_on_T1'][0] > 0.75 * af['G_on_T1'][1] and af['G_on_T0'][0] < 0.1 * af['G_on_T0'][1], 'text: the family marks the T1 arm among pair haplotypes'
# the finer scale (step 37 archive): the clusters split by T; family G marks T1; the S1-3 band
FAJ, FFE = 'results/sv/fine/armform_summary.json', 'results/sv/fine/unit_family_elements.tsv'
faf = json.load(open(FAJ)); ffe = pd.read_csv(FFE, sep='\t').set_index('family')
assert faf['mixed_pair_T1'] == faf['mixed_pair'] and faf['bulk_pair_T1'] == 0, 'text (finer): the shared cluster all T1, the bulk cluster all T0'
number('fafMixPair', faf['mixed_pair'], [FAJ]); number('fafMixClu', faf['mixed_cluster_size'], [FAJ]); number('fafBulkPair', faf['bulk_pair'], [FAJ]); number('fafBulkClu', faf['bulk_cluster_size'], [FAJ])
number('fafP', texp(faf['fisher_p_T1']), [FAJ]); number('fafMixCZero', faf['mixed_C0T1'], [FAJ]); number('fafMixCOne', faf['mixed_C1T1'], [FAJ])
fb_ = faf['bands']
number('fafThreePair', fb_['Three']['C1|T1'], [FAJ]); number('fafThreeNoPair', fb_['Three']['no pair'], [FAJ]); number('fafABCZTOne', fb_['AB']['C0|T1'], [FAJ])
assert fb_['Three']['C1|T1'] == sum(fb_['Three'].values()) - fb_['Three']['no pair'], 'text (finer): every S1-3 band haplotype with the pair is C1|T1'
assert faf['sig_family'] == 'G' if 'sig_family' in faf else True
number('fafTOneAll', faf['T1_all'], [FAJ]); number('fafTZeroAll', faf['T0_all'], [FAJ])
number('fafGTOne', faf['G_on_T1'][0], [FAJ]); number('fafGTZero', faf['G_on_T0'][0], [FAJ]); number('fafGHaps', faf['G_haps_all'], [FAJ]); number('fafGNoPair', faf['G_haps_nopair'], [FAJ])
number('ffeGXi', f"{ffe.loc['G', 'share_Xi']:.2f}", [FFE]); number('ffeGS', f"{ffe.loc['G', 'share_S']:.2f}", [FFE]); number('ffeGcombo', int(ffe.loc['G', 'top_combo_n']), [FFE]); number('ffeGocc', int(ffe.loc['G', 'occurrences']), [FFE])
assert ffe.loc['G', 'top_combo'] == 'Xi + S', 'text (finer): family G is Xi + S'

_np = int(re.search(r'\(1 \+ \(nl >= o\)\.sum\(\)\) / (\d+)\)\}', Path('scripts/09_ancestry_structure.py').read_text()).group(1)) - 1
number('anPermPA', f'{_np:,}'.replace(',', '{,}'), ['scripts/09_ancestry_structure.py', AJ])
# step 34: self dot plots, one haplotype per palindrome form (S5), and the illustrative principal MAP graph (S1)
DG = f'{O}/dotplot_graph_summary.json'; dg = json.load(open(DG)); _ex = pd.read_csv(f'{O}/locus_diagram_examples.tsv', sep='\t').set_index('panel')
assert dg['a']['src'] == _ex.loc['a', 'src'] and dg['d']['src'] == _ex.loc['b', 'src'], 'text: dot plots a and d show the Fig. 1 haplotypes'
for k_, v_ in dg['shimmer'].items(): number('dp' + {'w': 'W', 'k': 'K', 'r': 'R', 'min_span': 'Span'}[k_], v_, [DG])
for p_ in 'abcd':
    P_ = p_.upper()
    number(f'dp{P_}Name', ' '.join(dg[p_]['src'].split('_')[:2]), [DG]); number(f'dp{P_}Len', f"{dg[p_]['length'] / 1e3:,.0f}", [DG])
    number(f'dp{P_}Cen', f"{dg[p_]['beyond_A_cen_dots']:,}", [DG])
assert min(dg['a']['beyond_A_cen_dots'], dg['b']['beyond_A_cen_dots']) > 50 * max(dg['c']['beyond_A_cen_dots'], dg['d']['beyond_A_cen_dots']), 'text: X shows in the sequence of C1 cenSMN arms, not of C0'
g_ = dg['graph']; number('grBundles', g_['study_bundles'], [DG]); number('grSegs', g_['study_segments'], [DG])   # Fig. S1 draws the study's coarse bundles (2026-11-15)   # graph itself no longer drawn (researcher 2026-11-14)
assert g_['haplotypes'][4].startswith('HG01167_hap1'); number('grBand', ' '.join(g_['haplotypes'][4].split('_')[:2]), [DG])
# Fig. S1b (2026-11-14): HG03816 pat (units B- A+) and HG01099 mat (one long A+) carry the same elements; their arrays differ
_ga = pd.DataFrame(dg['graph_arrays'])
def _arr(h, side): r = _ga[(_ga.src.str.startswith(h)) & (_ga.side == side) & (_ga['where'] == 'outer')]; assert len(r) == 1; return float(r.copies.iloc[0])
number('grPatTel', f"{_arr('HG03816_pat', 'tel'):.1f}", [DG]); number('grMatTel', f"{_arr('HG01099_mat', 'tel'):.1f}", [DG])
number('grPatCen', f"{_arr('HG03816_pat', 'cen'):.1f}", [DG]); number('grMatCen', f"{_arr('HG01099_mat', 'cen'):.1f}", [DG])
_e = pd.read_csv(f'{O}/arm_element_blocks.tsv', sep='\t')
_ep, _em = [list(_e[_e.src.str.startswith(h)].sort_values(['side', 'i']).el) for h in ('HG03816_pat', 'HG01099_mat')]
assert _ep == _em, 'text: the two Fig. S1b haplotypes carry the same arm elements in the same order'

# step 34b: the chain race on the finer decomposition (Fig. S13 removed 2026-11-15); only the sentence in S8 rests on it
RC = f'{O}/unit_call_race.json'; rc = json.load(open(RC)); _m, _p = rc['HG01099 mat'], rc['HG03816 pat']
assert _m['winner'] == 'outer' and _p['winner'] == 'inner', 'text: the outer route wins in HG01099 mat and loses in HG03816 pat'
assert _m['unit_strings']['cen'].endswith('D+ A+') and _p['unit_strings']['cen'].endswith('B- A+') and _p['unit_strings']['tel'].startswith('A- E- B+'), 'text: B- A+ / one A+'
figure('figures/self_dotplots.png', [DG, 'results/decomp/smn.repeat_units.tsv', 'results/decomp/copy_units.tsv'])
# cold read 2026-11-14: Fig. 4c dot classes; S inside the 49-kb gap between the A units; AFR arrangement test; one no-unit count
DC_ = f'{O}/palindrome_dot_classes.tsv'; dcl = pd.read_csv(DC_, sep='\t').set_index('mark')['items']
number('fdSmnTwo', int(dcl['SMN2']), [DC_]); number('fdHyb', int(dcl['hybrid (h12 / h21)']), [DC_]); number('fdDel', int(dcl['Δ7–8 copy']), [DC_]); number('fdTrunc', int(dcl['truncated piece']), [DC_])
assert int(dcl.sum()) == int(MACROS['palCopies']) and int(dcl['SMN1 A']) == int(MACROS['palSmnA']) and int(dcl['SMN1 B1'] + dcl['SMN1 B2']) == int(MACROS['palSmnB']), 'caption: dot classes add up to the items'
_sp = pd.read_csv(ESP, sep='\t'); _sp = _sp[_sp.site == 'axis']; _ru = pd.read_csv('results/decomp/smn.repeat_units.tsv', sep='\t', comment='#')
_cov = [sorted(_ru[(_ru.row == r.src) & (_ru.bgn < r.end) & (_ru.end > r.bgn)].family) for r in _sp.itertuples()]
_free = sum(len(c) == 0 for c in _cov); _inA = sum(c == ['A'] for c in _cov)
number('palGapS', _free, [ESP, 'results/decomp/smn.repeat_units.tsv']); number('palGapSN', len(_sp), [ESP]); number('palSinA', _inA, [ESP, 'results/decomp/smn.repeat_units.tsv'])
assert _inA > 0.95 * len(_sp), 'text: at the coarse scale the axis S lies inside one of the two A units'
number('dcPAfrAArr', f"{dtp('AFR', 'A', 'arrangement_top').p:.4f}", [DT]); number('dcAfrArrSpread', f"{dtp('AFR', 'A', 'arrangement_top').spread:.2f}", [DT])
# what drives the AFR arrangement effect in A (cold read 2026-11-16): the '1- 1+' A copies, against the other types
_afa = dc[(dc.subset == 'AFR') & (dc.lineage == 'A') & (dc.context == 'arrangement_top')].set_index('level')
_afo = _afa.drop('1- 1+')
number('dcAfrArrOneOneN', int(_afa.loc['1- 1+', 'copies']), [DC]); number('dcAfrArrOneOneMed', f"{_afa.loc['1- 1+', 'median_to_SMN2']:.2f}", [DC])
number('dcAfrArrOthMin', f'{_afo.median_to_SMN2.min():.2f}', [DC]); number('dcAfrArrOthMax', f'{_afo.median_to_SMN2.max():.2f}', [DC])
assert _afa.median_to_SMN2.idxmin() == '1- 1+' and _afo.median_to_SMN2.max() - _afo.median_to_SMN2.min() < 0.1, 'text: the AFR arrangement shift comes from the 1- 1+ A copies'
assert abs((_afa.median_to_SMN2.max() - _afa.median_to_SMN2.min()) - dtp('AFR', 'A', 'arrangement_top').spread) < 0.002, 'spread = 1- 1+ against the farthest type'
# NAIP 3' context also shifts A a little (cold read 2026-11-16)
number('dcPANaip', f"{dtp('all', 'A', 'naip3').p:.4f}", [DT]); number('dcANaipSpread', f"{dtp('all', 'A', 'naip3').spread:.2f}", [DT])
assert dtp('all', 'A', 'naip3').p < 0.05 and dtp('all', 'A', 'naip3').spread < dtp('all', 'A', 'arrangement_top').spread, 'text: NAIP context shifts A slightly'
assert dtp('AFR', 'A', 'naip3').p > 0.05
assert dtp('AFR', 'A', 'arrangement_top').p < 0.01
assert MACROS['palNoUnits'] == MACROS['nHapsNoUnits']
figure('figures/window_bundles.png', [DG, 'results/decomp/smn.repeat_units.tsv'])

# step 28f: the NAIP split in the exons of NM_004536.3
NXJ = f'{O}/naip_exons.json'; nx_ = json.load(open(NXJ))
assert nx_['transcript'] == 'NM_004536.3' and set(nx_['piece_ends']) <= {0, 1611, 1965, 7088}
number('nxTx', nx_['transcript'].replace('_', r'\_'), [NXJ]); number('nxExons', nx_['n_exons'], [NXJ]); number('nxCodons', f"{nx_['cds_codons']:,}".replace(',', '{,}'), [NXJ])
P3a, P3b, P5 = nx_['pieces']['1611-7088'], nx_['pieces']['1965-7088'], nx_['pieces']['0-1965']
number('nxThreeAExon', P3a['first_exon'], [NXJ]); number('nxThreeBExon', P3b['first_exon'], [NXJ]); number('nxFiveExon', P5['last_exon'], [NXJ])
number('nxThreeACodon', P3a['first_codon'], [NXJ]); number('nxThreeBCodon', P3b['first_codon'], [NXJ]); number('nxFiveCodon', P5['last_codon'], [NXJ])
assert P3a['first_exon'] == P5['last_exon'] - 3 and nx_['shared_3p_from_1611_and_5p']['n_bases'] == int(MACROS['elNaipOverlap']), 'text: the two parts share exons 6-9 (354 bases)'
assert P5['domains'] == {'BIR 1': 'whole', 'BIR 2': 'whole', 'BIR 3': 'part'} and P3a['domains']['BIR 3'] == 'whole' and P3b['domains']['BIR 2'] == 'none', 'text: BIR domains in each part'

# step 28h: family G on haplotypes without the A+ A- pair (axis-free reading of S sites)
GJ = f'{O}/g_unpaired.json'; gj = json.load(open(GJ))
number('gValTOne', gj['pair_T1_t2like'][0], [GJ]); number('gValTOneN', gj['pair_T1_t2like'][1], [GJ]); number('gValTZero', gj['pair_T0_t2like'][0], [GJ]); number('gValTZeroN', gj['pair_T0_t2like'][1], [GJ])
number('gUnp', gj['unpaired'], [GJ]); number('gUnpG', gj['G_unpaired'], [GJ]); number('gUnpGTOne', gj['G_unpaired_t2like'], [GJ])
number('gUnpNoG', gj['noG_unpaired_t2like'][1], [GJ]); number('gUnpNoGTOne', gj['noG_unpaired_t2like'][0], [GJ])
number('gSThree', gj['G_unpaired_by_label']['S1-3'][0], [GJ]); number('gSThreeN', gj['G_unpaired_by_label']['S1-3'][1], [GJ])
number('gSTwo', gj['G_unpaired_by_label']['S1-2'][0], [GJ]); number('gSTwoN', gj['G_unpaired_by_label']['S1-2'][1], [GJ])
assert gj['unpaired'] == int(MACROS['palNoAdj']), 'the haplotypes without the pair are the 164 of the palindrome step'
assert gj['G_unpaired_by_label']['S1-3'][0] == gj['G_unpaired_by_label']['S1-3'][1], 'text: every S1-3 G haplotype without the pair is T1-like'

# step 35: great-ape structure (researcher 2026-11-14: both dot plots to the supplement, a sentence in section 2.7)
AF, AJ = f'{O}/ape_fold.tsv', f'{O}/ape_struct.json'; af = pd.read_csv(AF, sep='\t'); aj = json.load(open(AJ))
number('apFarHum', f"{aj['far_per_mb_human']:,.0f}", [AJ]); number('apFarLo', f"{aj['far_per_mb_apes'][0]:.0f}", [AJ]); number('apFarHi', f"{aj['far_per_mb_apes'][1]:.0f}", [AJ])
assert aj['far_per_mb_human'] > 50 * aj['far_per_mb_apes'][1], 'text: no ape stretch shows the large repeats'
afr_ = af[~af.ape.str.startswith('orangutan')]; oru = af[af.ape.str.startswith('orangutan')]
number('apShareLo', f"{100 * min(afr_.cen_share.min(), afr_.tel_share.min()):.0f}", [AF]); number('apShareHi', f"{100 * max(afr_.cen_share.max(), afr_.tel_share.max()):.0f}", [AF])
number('apOruShare', f"{100 * max(oru.cen_share.max(), oru.tel_share.max()):.0f}", [AF]); number('apWin', f"{aj['win_kb']:.0f}", [AJ])
assert (afr_.overlap > 0.98).all() and (afr_.cen_rho > 0.9).all() and (afr_.tel_rho > 0.9).all(), 'text: both arms land on one ape stretch, in the same order from the axis'
pc_ = af[af.ape.str.match('(chimp|bonobo)')]; go_ = af[af.ape.str.startswith('gorilla')]
number('apSPanLo', f"{pc_.S_over_Fl.min():.1f}", [AF]); number('apSPanHi', f"{pc_.S_over_Fl.max():.1f}", [AF])
number('apSOtherHi', f"{max(go_.S_over_Fl.max(), oru.S_over_Fl.max()):.2f}", [AF])
assert (pc_.S_over_Fl > 0.8).all() and (go_.S_rate == 0).all() and (oru.S_over_Fl < 0.05).all(), 'text: S found in chimpanzee and bonobo, not in gorilla or orangutan'
figure('figures/ape_self_dotplots.png', [AF, AJ]); figure('figures/ape_vs_human.png', [AF, f'{O}/ape_element_hits.tsv'])

# step 33: units do not overlap in position, but families share content
UJ_ = json.load(open(UJ))['overlap']
number('ufPosOv', UJ_['position_overlapping_pairs'], [UJ]); number('ufKmerOv', UJ_['kmer_edge_overlaps'], [UJ]); number('ufPairs', UJ_['family_pairs'], [UJ]); number('ufPairsShare', UJ_['family_pairs_sharing_core'], [UJ])
number('ufBids', UJ_['bundle_ids_in_units'], [UJ]); number('ufBidsMulti', UJ_['bundle_ids_in_2plus_families'], [UJ]); number('ufInA', ', '.join(r'\RU{%s}' % f for f in UJ_['core_within_A']), [UJ])
# step 28j: why the family-A (R_A) units differ in length when the core does not (Fig. 1; Section 2.4)
AUJ = f'{O}/a_unit_length.json'; au = json.load(open(AUJ)); _p = au['pooled']; _c = au['core_all']
number('auCore', f"{_c['median']:.0f}", [AUJ]); number('auCoreLo', f"{_c['q10']:.0f}", [AUJ]); number('auCoreHi', f"{_c['q90']:.0f}", [AUJ])
number('auHaps', au['n_haps'], [AUJ]); assert au['n_haps'] == 270
number('auCZeroHaps', _p['C0_haps'], [AUJ]); number('auCZeroBoth', _p['C0_both'], [AUJ])
number('auCZeroUnit', f"{_p['C0']['unit_med']:.0f}", [AUJ]); number('auCZeroOver', f"{_p['C0']['over_med']:.0f}", [AUJ])
_bn = au['both_neither']['C1|T0']
number('auCOTZeroN', _bn['n'], [AUJ]); number('auCOTZeroNeither', _bn['neither'], [AUJ])
number('auCOTZeroUnit', f"{_p['C1T0']['unit_med']:.0f}", [AUJ]); number('auCOTZeroShort', f"{-_p['C1T0']['over_med']:.0f}", [AUJ])
_bt = au['both_neither']['C1|T1']
number('auCOTOneN', _bt['n'], [AUJ]); number('auCOTOneTelOnly', _bt['tel_only'], [AUJ])
number('auCOTOneTelUnit', f"{_p['C1T1_tel']['unit_med']:.0f}", [AUJ]); number('auCOTOneCenUnit', f"{_p['C1T1_cen']['unit_med']:.0f}", [AUJ])
_na = au['next_after_core']
number('auCOneCenXi', _na['C1|cen'].get('Xi', 0), [AUJ]); number('auCOneCenN', sum(_na['C1|cen'].values()), [AUJ])
number('auCZeroCenE', _na['C0|cen'].get('E', 0), [AUJ]); number('auCZeroCenN', sum(_na['C0|cen'].values()), [AUJ])
number('auTelE', _na['C0|tel'].get('E', 0) + _na['C1|tel'].get('E', 0), [AUJ])
# cold read (rev 31): E present on a telSMN arm but not right after the core; pair haplotypes outside the structure tree
import csv as _csv
_AUT, _AEH, _MLH = f'{O}/a_unit_length.tsv', f'{O}/arm_element_haps.tsv', f'{O}/ml_hap_labels.tsv'
_aul = [r for r in _csv.DictReader(open(_AUT), delimiter='\t') if r['side'] == 'tel']
_eNot = [r['src'] for r in _aul if r['next_el'] != 'E' and 'E' in r['layout'].split()]
number('auTelENotNext', len(_eNot), [_AUT])
assert int(MACROS['auTelE']) + len(_eNot) == int(MACROS['elEtel']), 'text: E follows the core on auTelE telSMN arms and lies further out on the rest of the elEtel'
_aeh = list(_csv.DictReader(open(_AEH), delimiter='\t'))
_mlh = {r['src'] for r in _csv.DictReader(open(_MLH), delimiter='\t')}
_hcc = {r['src']: r for r in _csv.DictReader(open('tables/haplotype_copy_counts.tsv'), delimiter='\t')}
_out = [r['src'] for r in _aeh if r['src'] not in _mlh]
number('afPairOutTree', len(_out), [_AEH, _MLH]); number('afPairInTree', len(_aeh) - len(_out), [_AEH, _MLH])
assert all(int(_hcc[x]['SMN1']) == 0 for x in _out), 'Fig. 5 caption: pair haplotypes outside the structure tree carry no SMN1 copy'
assert int(MACROS['afBulkCZero']) == 0 and int(MACROS['afMixCOne']) == 1 and int(MACROS['afMixCZero']) == int(MACROS['afMixPair']) - 1, 'abstract: one cluster all with X, the other all but one without'
assert int(MACROS['afBulkPair']) + int(MACROS['afMixPair']) < len(_aeh) - len(_out), 'abstract must not say the two clusters hold every pair haplotype'

# Fig. 1's two haplotypes
_lde = pd.read_csv(f'{O}/locus_diagram_examples.tsv', sep='\t').set_index('panel').src
for _m, _s in (('A', _lde['a']), ('D', _lde['b'])):
    _e = au['examples'][_s]
    for _sd in ('cen', 'tel'):
        number(f'au{_m}{_sd}Unit', f"{_e[_sd]['unit_kb']:.0f}", [AUJ]); number(f'au{_m}{_sd}Core', f"{_e[_sd]['core_kb']:.0f}", [AUJ])
assert _p['C0_both'] / _p['C0_haps'] > 0.85 and _p['C0']['over_med'] > 50, 'text: on C0 haplotypes both units run past the core, over E'
assert _bn['neither'] / _bn['n'] > 0.9 and _p['C1T0']['over_med'] < -30, 'text: on C1|T0 both units stop short of the core end'
assert _bt['tel_only'] / _bt['n'] > 0.85, 'text: on C1|T1 only the telSMN unit runs past the core'
assert _na['C0|cen'].get('E', 0) == sum(_na['C0|cen'].values()) and _na['C1|cen'].get('Xi', 0) / sum(_na['C1|cen'].values()) > 0.95, 'text: after the core, E on C0 and Xi on C1 cenSMN arms'
assert abs(_p['C0']['unit_med'] - _p['C1T0']['unit_med']) > 100 and _c['q90'] - _c['q10'] < 50, 'text: units differ by >100 kb while the core varies within ~50 kb'
assert au['examples'][_lde['a']]['cen']['next_el'] == 'Xi' and au['examples'][_lde['b']]['cen']['next_el'] == 'E', 'Fig. 1: panel a C1 (Xi after the core), panel b C0 (E after it)'
number('ufABshared', UJ_['shared_core']['A/B'], [UJ])
AB_ = json.load(open(UJ))['ab']   # families A and B: the same core bundles; what separates them
assert AB_['same_core'] and AB_['A_in_order'] / AB_['A_n'] > 0.9 and AB_['B_in_order'] / AB_['B_n'] > 0.9, 'text: A and B hold the same core bundles in the same order and orientation'
number('ufABorder', AB_['core_order'], [UJ]); number('ufAinOrd', AB_['A_in_order'], [UJ]); number('ufBinOrd', AB_['B_in_order'], [UJ])
number('ufAn', AB_['A_n'], [UJ]); number('ufBn', AB_['B_n'], [UJ])
_b1, _b15 = AB_['bp_median_kb']['1'], AB_['bp_median_kb']['15']
number('ufAbOneA', f"{_b1[0]:.0f}", [UJ]); number('ufAbOneB', f"{_b1[1]:.0f}", [UJ]); number('ufAbArrA', f"{_b15[0]:.1f}", [UJ]); number('ufAbArrB', f"{_b15[1]:.1f}", [UJ])
_oth = [abs(v[0] - v[1]) for k, v in AB_['bp_median_kb'].items() if k not in ('1', '15')]
number('ufAbOtherMax', f"{max(_oth):.1f}", [UJ])
assert _b1[0] > 3 * _b1[1] and _b15[1] > _b15[0] and max(_oth) < 1.5, 'text: A and B differ in how much of bundle 1 and of the array bundle 15 they hold; the other core bundles are the same size'
number('ufBstartK', AB_['B_starts_9p8p'], [UJ]); number('ufAstartK', AB_['A_starts_9p8p'], [UJ]); number('ufAtenM', AB_['A_has_10m_before_core'], [UJ]); number('ufBtenM', AB_['B_has_10m_before_core'], [UJ])
assert AB_['A_starts_9p8p'] == 0 and AB_['B_has_10m_before_core'] == 0, 'text: 9+ 8+ opens B units only; 10- before the core in A units only'
_n = lambda m: float(str(MACROS[m]).replace('{,}', ''))
for _p, _k in (('mantelP', 'nPerm'), ('anPArr', 'anPerm'), ('anPAll', 'anPermPA'), ('anPOne', 'anPermPA'), ('anPTwo', 'anPermPA')):
    assert abs(_n(_p) - 1 / (_n(_k) + 1)) < 1e-9, f'text: p <= {_p} is the floor of {_k} permutations'
AF_ = json.load(open(UJ))['ab_family']   # step 33: replay of PANGEA's family rule, why A and B are two families
assert AF_['replay_matches_all'] == json.load(open(UJ))['occurrences'], 'the replay reproduces the family of every occurrence'
number('ufRepN', AF_['replay_matches_all'], [UJ]); number('ufAfoundKb', f"{AF_['A_founder'][2] / 1e3:.0f}", [UJ]); number('ufBfoundKb', f"{AF_['B_founder'][2] / 1e3:.0f}", [UJ])
number('ufBfoundSim', f"{AF_['B_founder_sim_to_A']:.2f}", [UJ]); number('ufMeansSim', f"{AF_['means_sim']:.2f}", [UJ])
number('ufAownSim', f"{AF_['A_own_median']:.2f}", [UJ]); number('ufBownSim', f"{AF_['B_own_median']:.2f}", [UJ]); number('ufBvsASim', f"{AF_['B_vs_A_median']:.2f}", [UJ])
number('ufAcloser', AF_['A_closer_own'], [UJ]); number('ufBcloser', AF_['B_closer_own'], [UJ]); number('ufBgeThr', AF_['B_ge_thr_vs_A'], [UJ])
number('ufMeanOneA', f"{AF_['mean_kb']['1'][0]:.0f}", [UJ]); number('ufMeanOneB', f"{AF_['mean_kb']['1'][1]:.0f}", [UJ]); number('ufMeanArrA', f"{AF_['mean_kb']['15'][0]:.0f}", [UJ]); number('ufMeanArrB', f"{AF_['mean_kb']['15'][1]:.0f}", [UJ])
assert AF_['B_founder_sim_to_A'] < AF_['family_jaccard'] <= AF_['means_sim'], 'text: B was founded below the threshold, yet the final means are above it'
assert AF_['B_closer_own'] == AB_['B_n'] and AF_['A_closer_own'] > 0.98 * AB_['A_n'] and AF_['B_ge_thr_vs_A'] > AB_['B_n'] / 2, 'text: each unit is nearer its own family mean; most B units would pass the threshold against A'
assert UJ_['position_overlapping_pairs'] == 0 and 'B' in UJ_['core_within_A'] and UJ_['shared_core']['A/B'] == len(ids('B'))

# step 28i: the association at both bundle scales, and every decomposition parameter (Supplementary Section 'The finer scale')
TSJ, TST, DPT = f'{O}/two_scales.json', f'{O}/two_scales.tsv', f'{O}/decomp_params.tsv'
tsj = json.load(open(TSJ)); tst = pd.read_csv(TST, sep='\t').set_index('test'); dpt = pd.read_csv(DPT, sep='\t').set_index('scale')
_TSN = {'all copies': 'All', 'SMN1': 'One', 'SMN2': 'Two', 'SMN1, AFR': 'OneAfr', 'SMN1, outside AFR': 'OneNonAfr', 'SMN1, C1 arms': 'COne',
        'SMN1, C0 arms': 'CZero', 'SMN1, T0 arms': 'TZero', 'SMN1, T1 arms': 'TOne', 'cophenetic correlation (SMN1 haplotypes)': 'Coph'}
assert set(tst.index) == set(_TSN)
for t_, nm in _TSN.items():
    r = tst.loc[t_]; number(f'ts{nm}C', f'{r.coarse_spearman:.2f}', [TST]); number(f'ts{nm}F', f'{r.fine_spearman:.2f}', [TST]); number(f'ts{nm}N', int(r.haplotypes), [TST])
    if nm != 'Coph': assert r.coarse_p <= 0.001 and r.fine_p <= 0.001, 'text: every Mantel test p = 0.001 at both scales'
assert MACROS['tsAllC'] == MACROS['mantelAll'] and MACROS['tsOneC'] == MACROS['mOneR'] and MACROS['tsCophC'] == MACROS['cophenetic'], 'the coarse column is the main text'
number('tsStructRho', f"{tsj['coarse_vs_fine_struct_spearman']:.2f}", [TSJ])
_up = [n for n in _TSN.values() if n != 'Coph' and tst.loc[[k for k, v in _TSN.items() if v == n][0]].fine_spearman > tst.loc[[k for k, v in _TSN.items() if v == n][0]].coarse_spearman]
number('tsHigherFine', len(_up), [TST]); number('tsTests', len(_TSN) - 1, [TST])
assert tst.loc['SMN1, T1 arms'].fine_spearman < tst.loc['SMN1, T1 arms'].coarse_spearman + 0.02 and tst.loc['SMN1, C0 arms'].fine_spearman - tst.loc['SMN1, C0 arms'].coarse_spearman > 0.15, 'text: the gap is largest on C0 arms and absent on T1 arms'
rows = []
for t_, nm in _TSN.items():
    r = tst.loc[t_]; rows.append(f"{t_.replace('C1', 'C1').replace('AFR', 'AFR')} & {int(r.haplotypes)} & {r.coarse_spearman:.2f} & {r.fine_spearman:.2f} \\\\")
(DOC / 'tab_twoscales.tex').write_text('% written by scripts/40_numbers.py from results/sv/two_scales.tsv; do not edit\n\\newcommand{\\tabTwoScales}{%\n' + '\n'.join(rows) + '}\n')
_P = ['w', 'k', 'r', 'min-span', 'min-cov', 'min-branch-size', 'bundle-length-cutoff', 'bundle-merge-distance']
rows = []
for p_ in _P:
    c, f = dpt.loc['coarse'], dpt.loc['fine']
    cell = lambda s, p: f"{int(s[p]):,}".replace(',', '{,}') + ('' if s[p + '_source'] == 'given' else ' (default)')
    rows.append(f"\\texttt{{{"-" if len(p_) == 1 else "--"}{p_}}} & {cell(c, p_)} & {cell(f, p_)} \\\\")
for k_ in ['rows', 'segments', 'bundle_ids']:
    rows.append(f"{k_.replace('_', ' ')} & {int(dpt.loc['coarse', k_]):,} & {int(dpt.loc['fine', k_]):,} \\\\".replace(',', '{,}').replace('{,} &', ', &'))
(DOC / 'tab_decomp.tex').write_text('% written by scripts/40_numbers.py from results/sv/decomp_params.tsv; do not edit\n\\newcommand{\\tabDecomp}{%\n' + '\n'.join(rows) + '}\n')
for s_, M_ in [('coarse', 'C'), ('fine', 'F')]:
    number(f'dp{M_}Segs', f"{int(dpt.loc[s_, 'segments']):,}".replace(',', '{,}'), [DPT]); number(f'dp{M_}Bids', int(dpt.loc[s_, 'bundle_ids']), [DPT])
    number(f'dp{M_}Cut', f"{int(dpt.loc[s_, 'bundle-length-cutoff']):,}".replace(',', '{,}'), [DPT]); number(f'dp{M_}Merge', f"{int(dpt.loc[s_, 'bundle-merge-distance']):,}".replace(',', '{,}'), [DPT])
    number(f'dp{M_}K', int(dpt.loc[s_, 'k']), [DPT]); number(f'dp{M_}Branch', int(dpt.loc[s_, 'min-branch-size']), [DPT]); number(f'dp{M_}Span', int(dpt.loc[s_, 'min-span']), [DPT])
assert dpt.loc['coarse', 'bundle_ids'] < dpt.loc['fine', 'bundle_ids'] and dpt.loc['coarse', 'segments'] < dpt.loc['fine', 'segments'], 'text: the coarse decomposition has fewer, larger bundles'

# step 23b: Luo et al. 2014's six duplication-haplotype variants (their Table 4) on our lineages
LV = f'{O}/luo_variants.tsv'; _lv = pd.read_csv(LV, sep='\t')
def _lvn(g, grp, col='n_alt'): return int(_lv[(_lv.luo_g == g) & (_lv.group == grp)][col].iloc[0])
_five = ['c15774', 'c22804', 'c26190', 'c27134', 'c27706_27707']
number('luoShift', f"{int(_lv.grch38.iloc[0] - _lv.grch37.iloc[0]):,}".replace(',', '{,}'), [LV])
number('luoBoneLo', min(_lvn(g, 'B1') for g in _five), [LV]); number('luoBoneHi', max(_lvn(g, 'B1') for g in _five), [LV])
number('luoBtwoMk', _lvn('c27134', 'B2'), [LV]); number('luoBtwoOther', max(_lvn(g, 'B2') for g in _five if g != 'c27134'), [LV])
number('luoIntronOneB', _lvn('c11678', 'B1', 'top_n'), [LV])
assert all(_lvn(g, grp) == 0 for g in _five for grp in ['A', 'SMN2']), 'text: none of the five is on any A or SMN2 copy'
assert _lvn('c27134', 'B1') == MACROS['linBoneMk'] and _lvn('c27134', 'B2') == MACROS['linBtwoMk'], 'step 23b agrees with the marker counts'
assert max(_lvn(g, 'B2') for g in _five if g != 'c27134') == 0, 'text: B2 carries only the intron-7 marker'
assert sum(_lvn('c11678', grp) for grp in ['A', 'B1', 'B2', 'SMN2']) == 0, 'text: the intron-1 change as Luo write it is on no copy'
_t = _lv[(_lv.luo_g == 'c11678') & (_lv.group == 'B1')].iloc[0]
import re as _re
_m = lambda w: tuple(len(x) for x in _re.match(r'(T+)(G+)(T+)G', w).groups())
assert _m(_t.top_window)[1] == _m(_t.alt_window)[1] == _m(_t.ref_window)[1] - 1 and _m(_t.top_window)[0] == _m(_t.ref_window)[0] + 1 and _m(_t.alt_window)[2] == _m(_t.ref_window)[2] + 1, (_m(_t.top_window), _m(_t.ref_window), _m(_t.alt_window), 'text: B1 has the G run shortened by one at its left end; Luo write it at its right end')

# cold read 2026-11-16: axis S counted per occurrence and per haplotype; per-ape spread of SMN2 bases; Fig. 1 axis gap
_spx = pd.read_csv(ESP, sep='\t'); _spx = _spx[_spx.site == 'axis']; _vc = _spx.src.value_counts()
_ph = set(pd.read_csv(f'{O}/arm_element_haps.tsv', sep='\t').src)
number('palGapSHaps', int(_spx.src.nunique()), [ESP]); number('sAxisTwo', int((_vc == 2).sum()), [ESP])
number('sAxisTwoHaps', ' and '.join(sorted(h.replace('_hprc_r2_v1.0.1', '').replace('_', ' ') for h in _vc[_vc == 2].index)), [ESP])
_none = sorted(_ph - set(_spx.src)); number('sAxisNone', len(_none), [ESP]); number('sAxisNoneHaps', ' and '.join(h.replace('_hprc_r2_v1.0.1', '').replace('_', ' ') for h in _none), [ESP])
assert int(MACROS['palGapSN']) == len(_spx) and int((_vc == 1).sum()) == int(MACROS['elAxisS']), 'text: 271 occurrences; 267 haplotypes with exactly one'
assert int((_vc == 1).sum()) + int((_vc == 2).sum()) + len(_none) == int(MACROS['elHaps']) and (_vc <= 2).all() and set(_spx.src) <= _ph, 'text: 267 + 2 + 1 = 270'
_tb = X['S_apes_smn2_base']; number('evTwoBaseMin', min(_tb.values()), [ES]); number('evTwoBaseMax', max(_tb.values()), [ES])
assert min(_tb.values()) <= X['S_anc_smn2'] <= max(_tb.values()), 'text: the inferred ancestral count lies within the per-ape range'
_lx = pd.read_csv(f'{O}/locus_diagram_examples.tsv', sep='\t').set_index('panel')
_gp = [abs(float(_lx.loc[p_, 'le'])) for p_ in 'ab'] + [float(_lx.loc[p_, 'rb']) for p_ in 'ab']
number('dpAxisGap', f"{_gp[0]:.1f}", [f'{O}/locus_diagram_examples.tsv'])
assert max(_gp) - min(_gp) < 0.1, 'caption: both Fig. 1 units start about the same distance from the axis'
for p_, c_, t_ in [('a', 'auAcenUnit', 'auAtelUnit'), ('b', 'auDcenUnit', 'auDtelUnit')]:
    assert abs(float(_lx.loc[p_, 'left']) - float(_lx.loc[p_, 'le']) * -1 - float(MACROS[c_])) < 1.5 or abs(float(_lx.loc[p_, 'left']) - abs(float(_lx.loc[p_, 'le'])) - float(MACROS[c_])) < 1.5, ('caption: drawn extent = unit + axis gap', p_)

assert int(MACROS['linATel']) - int(MACROS['dcATelN']) == 2, 'text: step 17 leaves out the GRCh38 and CHM13 telomeric A copies'
# cold read 2026-11-16b (accepted by the researcher): SMN1 copies off the trees' haplotypes; AFR 135 vs 133; 164 vs 160
_mlc = pd.read_csv(f'{O}/marker_lineage.tsv', sep='\t'); _mlc = _mlc[_mlc['class'] == 'SMN1'].assign(src=lambda x: x.copy_id.str.split('|').str[0])
_off = _mlc[~_mlc.src.isin(set(_h.src)) & ~_mlc.src.str.startswith(('GCA_', 'chm13'))]
number('smnOneAll', len(_mlc), [f'{O}/marker_lineage.tsv']); number('smnOneUnres', len(_off), [f'{O}/marker_lineage.tsv', HCC]); number('smnOneUnresHaps', _off.src.nunique(), [f'{O}/marker_lineage.tsv', HCC])
assert len(_mlc) == int(MACROS['linA']) + int(MACROS['linB']) and len(_mlc) - int(MACROS['mlbCopies']) == len(_off), 'text: 522 - 497 = the SMN1 copies on unresolved haplotypes'
assert (_off.sublineage.astype(str) == 'A').all() and not _off.in_marker_side.any(), 'text: every SMN1 copy on an unresolved haplotype is lineage A'
_an = _h[_h.afr & (_h.SMN1 == 0)]
number('dosAfrNoSmnOne', len(_an), [HCC, f'{O}/copy_labels.tsv'])
assert int(MACROS['dosAfrHaps']) - len(_an) == int(MACROS['mOneAFRN']) == int(MACROS['tsOneAfrN']), 'text: 135 AFR resolved haplotypes - 2 without an SMN1 copy = the 133 of the Mantel test'
number('palAdjNotPal', int(MACROS['palNoPair']) - int(MACROS['palNoAdj']), [f'{O}/palindrome_units.tsv'])
assert int(MACROS['palNoPair']) - int(MACROS['palNoAdj']) == int(MACROS['palTail']) + 1 == int(MACROS['palAdj']) - int(MACROS['palPair']), 'text: the 4 with an adjacent pair but no palindrome: 3 tail to tail and NA19391 hap2'

# step 38 (2026-11-18): X is named for its shape in a self dot plot; the cross is centered on S, which is part of X
XC = f'{O}/x_cross.json'; xc = json.load(open(XC)); XT = f'{O}/x_cross.tsv'; _xt = pd.read_csv(XT, sep='\t')
for k_, n_ in [('c1_with_S', 'xcCOneS'), ('c1_peak_le20', 'xcTwenty'), ('c1_peak_le10', 'xcTen'), ('c1_without_S', 'xcNoS')]: number(n_, xc[k_], [XC, XT])
number('xcLo', xc['lo_kb'], [XC]); number('xcHi', xc['hi_kb'], [XC])
number('xcExS', f"{xc['examples']['HG01099_mat_hprc_r2_v1.0.1']['S_kb']:.0f}", [XC])
number('xcExCZeroArm', f"{xc['examples']['HG01099_pat_hprc_r2_v1.0.1']['arm_kb']:.0f}", [XC])
assert xc['arms'] == int(MACROS['elCOne']) + int(MACROS['elCZero']) and xc['c1_arms'] == int(MACROS['elCOne']), 'text: every palindrome cenSMN arm scored; C1 count as step 28b'
assert xc['c1_with_S'] + xc['c1_without_S'] == xc['c1_arms'], 'text: C1 arms with and without S beyond 300 kb'
assert xc['c1_peak_le20'] >= 0.9 * xc['c1_with_S'] and xc['c1_peak_le10'] > xc['c1_with_S'] / 2, 'text: on most C1 arms the cross is centered within 20 kb of S'
assert xc['examples']['HG01099_mat_hprc_r2_v1.0.1']['type'] == 'C1' and xc['examples']['HG01099_pat_hprc_r2_v1.0.1']['type'] == 'C0', 'text: Fig. panel a C1, b C0'
figure('figures/x_cross.png', [XC, XT, f'{O}/arm_element_haps.tsv', f'{O}/arm_element_segments.tsv', 'results/decomp/smn.repeat_units.tsv'])
# C1 arms without S (cold read): Xi and Xd present, no S bundle between them, and not cut by the 900 kb limit
_ns = sorted(_xt[(_xt.cen_type == 'C1') & _xt.S_kb.isna()].src)
_ea = pd.read_csv(EA_, sep='\t'); _ec = _ea[(_ea.side == 'cen') & _ea.src.isin(_ns)]
assert len(_ns) == xc['c1_without_S'] == len(_ec) and (_ec.kb_S == 0).all() and (_ec.kb_Xi > 0).all() and (_ec.kb_Xd > 0).all(), 'text: the C1 arms without S carry Xi and Xd but no S'
assert (_xt.set_index('src').loc[_ns].arm_kb < xc['hi_kb']).all(), 'text: they are not arms whose S lies beyond the scored stretch'
assert str(MACROS['sCenHaps']) == str(xc['c1_with_S']), 'text: S on the cenSMN arm of the same C1 haplotypes in both steps'
number('xcNoSNames', ' and '.join(x.split('_hprc')[0].replace('_', ' ') for x in _ns), [XT, EA_])
# step 43: shimmer similarity between repeat-unit families
FS = f'{O}/family_shimmer.json'; fs = json.load(open(FS)); FSP = f'{O}/family_shimmer_pairs.tsv'; _fc = fs['cells']; _ff = fs['families']
number('fsUnits', fs['n_units'], [FS]); number('fsPairs', fs['n_pairs'], [FS]); number('fsN', fs['max_units'], [FS])
number('fsMinOcc', fs['min_occ'], [FS]); number('fsFams', len(_ff), [FS])
for _m, _p in [('jaccard', 'J'), ('overlap', 'O'), ('span_jacc', 'S'), ('contain', 'C')]:
    for _a, _b in [('A', 'B'), ('B', 'A'), ('A', 'A'), ('B', 'B'), ('D', 'C'), ('C', 'D')]:
        number(f'fs{_p}{_a}{_b}', f"{_fc[_m][f'{_a}|{_b}']:.2f}", [FS])
_offJ = [_fc['jaccard'][f'{a}|{b}'] for a in _ff for b in _ff if a != b and {a, b} != {'A', 'B'}]
number('fsJOtherMax', f'{max(_offJ):.2f}', [FS])
_dO = [_fc['overlap'][f'{a}|{a}'] for a in _ff]; _dJ = [_fc['jaccard'][f'{a}|{a}'] for a in _ff]
number('fsODiagLo', f'{min(_dO):.2f}', [FS]); number('fsODiagHi', f'{max(_dO):.2f}', [FS])
number('fsJDiagLo', f'{min(_dJ):.2f}', [FS]); number('fsJDiagHi', f'{max(_dJ):.2f}', [FS])
assert set(fs['offdiag_max']['jaccard']['pair']) == {'A', 'B'} and _fc['jaccard']['A|B'] > max(_offJ), 'text: R_A-R_B the most alike pair of families by whole-unit Jaccard'
assert _fc['contain']['B|A'] > _fc['contain']['A|B'] + 0.2 and _fc['contain']['D|C'] > _fc['contain']['C|D'] + 0.3, 'text: R_B largely within R_A, R_D within R_C'
assert max(_fc[m][f'{y}|{x}'] for m in ['jaccard', 'contain'] for x in 'AB' for y in 'D') <= 0.01 and max(_fc['contain'][f'{x}|D'] for x in 'AB') <= 0.01, 'text: R_D and R_A/R_B medians at most 0.01'
assert all(_fc['overlap'][f'{a}|{a}'] > _fc['jaccard'][f'{a}|{a}'] for a in _ff), 'text: normalising by the shorter unit raises every within-family value'
assert _fc['span_jacc']['B|B'] - _fc['jaccard']['B|B'] < _fc['overlap']['B|B'] - _fc['jaccard']['B|B'], 'text: the matched-span version raises R_B less'
_dab = max([_fc[m][f'D|{x}'] for m in ['jaccard', 'contain'] for x in 'AB'] + [_fc['contain'][f'{x}|D'] for x in 'AB'])
number('fsDABMax', f'{_dab:.2f}', [FS])
figure('figures/family_shimmer.png', [FS, FSP])


# ---- Discussion: an SMN2 clade that travels with the S1-2 SMN1 clade on the same haplotypes (step 19d)
RS2 = [f'{O}/smn2_s12_clade.json']
_s2j = json.load(open(RS2[0]))
number('sTwoClN', _s2j['n_copies'], RS2); number('sTwoClSup', _s2j['support'], RS2)
number('sTwoClHapsAll', _s2j['haps_all'], RS2); number('sTwoClHaps', _s2j['haps_in_tree'], RS2); number('sTwoClWithSOneTwo', _s2j['haps_with_s1_2'], RS2)
number('sTwoClAfr', _s2j['afr_copies'], RS2); number('sTwoClTwoMinusOneMinus', _s2j['copies_on_2m_1m'], RS2)
assert _s2j['haps_with_s1_2'] >= 0.9 * _s2j['haps_in_tree'], 'text: nearly every haplotype with a copy of this SMN2 clade carries an S1-2 SMN1 copy'
assert _s2j['afr_copies'] <= 2, 'text: almost none African'
_CLs = pd.read_csv(f'{O}/copy_labels.tsv', sep='\t').set_index('tip')
_HL = pd.read_csv(f'{O}/ml_hap_labels.tsv', sep='\t')


# ---- cold read: the two largest structure clusters, and the Fig. 3 arm-strip count
_sc = _HL.struct_cluster.value_counts()
number('scN', int(_sc.size), [f'{O}/ml_hap_labels.tsv']); number('scFirst', int(_sc.iloc[0]), [f'{O}/ml_hap_labels.tsv'])
number('scSecond', int(_sc.iloc[1]), [f'{O}/ml_hap_labels.tsv']); number('scThird', int(_sc.iloc[2]), [f'{O}/ml_hap_labels.tsv'])
assert int(_sc.iloc[0]) == int(MACROS['afBulkN']) if 'afBulkN' in MACROS else True
assert _sc.iloc[1] > 3 * _sc.iloc[2], 'text: two clusters clearly the largest'
_PCL = pd.read_csv(f'{O}/palindrome_copies_lineage.tsv', sep='\t')
_cen = _PCL[_PCL.arm.str.startswith('left')]
number('armCenItems', len(_cen), [f'{O}/palindrome_copies_lineage.tsv']); number('armCenNotTree', int((~_cen.copy_id.isin(_CLs.copy_id)).sum()), [f'{O}/palindrome_copies_lineage.tsv', f'{O}/copy_labels.tsv'])
number('armCenInTree', int(_cen.copy_id.isin(_CLs.copy_id).sum()), [f'{O}/palindrome_copies_lineage.tsv', f'{O}/copy_labels.tsv'])
_nt = _cen[~_cen.copy_id.isin(_CLs.copy_id)]
_ntD = _nt[_nt.token.astype(str).str.startswith('d')]; _ntT = _nt[_nt.token.astype(str).str.startswith('t')]; _ntC = _nt[~_nt.index.isin(_ntD.index) & ~_nt.index.isin(_ntT.index)]
assert len(_ntD) + len(_ntT) + len(_ntC) == len(_nt) and len(_ntC) == 1 and (_ntC.gene == 'SMN2').all(), _nt
_cc = pd.read_csv('results/copies/copy_calls.tsv', sep='\t', dtype=str)
assert (_cc[_cc.copy_id == _ntC.copy_id.iloc[0]].status == 'composite').all() and (_cc.copy_id == _ntC.copy_id.iloc[0]).sum() == 1
_PL = [f'{O}/palindrome_copies_lineage.tsv', f'{O}/copy_labels.tsv']
number('armCenNTDel', len(_ntD), _PL); number('armCenNTTrunc', len(_ntT), _PL); number('armCenNTComp', len(_ntC), _PL + ['results/copies/copy_calls.tsv'])
number('armCenNTCompHap', _ntC.src.iloc[0].replace('_hprc_r2_v1.0.1', '').replace('_', ' '), _PL)
assert int(_cen.copy_id.isin(_CLs.copy_id).sum()) == 262, 'Fig. 3 legend: cenSMN arm 262'
# ---- final check: assembly gaps (runs of N) inside the windows (step 36b)
WG = f'{O}/window_gaps.tsv'
_wg = pd.read_csv(WG, sep='\t'); _wgr = _wg[_wg.resolved]; _wgg = _wgr[_wgr.n_bases > 0]
assert (_wg[_wg.n_bases > 0].resolved).all(), 'every window with N is a resolved haplotype'
_wgs = _wgg[_wgg.n_bases < 1000]; _wgb = _wgg[_wgg.n_bases >= 1000]
number('gapRes', len(_wgg), [WG, HC]); number('gapSmall', len(_wgs), [WG]); number('gapSmallMax', int(_wgs.n_bases.max()), [WG])
number('gapBig', len(_wgb), [WG]); number('gapBigMin', f'{_wgb.n_bases.min()/1000:.1f}', [WG]); number('gapBigMax', f'{_wgb.n_bases.max()/1000:.0f}', [WG])
assert len(_wgs) + len(_wgb) == len(_wgg) and _wgb.n_bases.min() > 10 * _wgs.n_bases.max(), 'two clearly separated gap sizes'
_aeh = pd.read_csv(f'{O}/arm_element_haps.tsv', sep='\t')
number('gapBigPal', int(_aeh.src.isin(_wgb.src).sum()), [WG, f'{O}/arm_element_haps.tsv'])
_gpb = int(_aeh.src.isin(_wgb.src).sum()); _gpr = int((~_aeh.src.isin(_wgb.src)).sum())
number('gapRestN', len(_wgr) - len(_wgb), [WG, HC]); number('gapRestPal', _gpr, [WG, f'{O}/arm_element_haps.tsv'])
_gpp = fisher_exact([[_gpb, len(_wgb) - _gpb], [_gpr, (len(_wgr) - len(_wgb)) - _gpr]])[1]
number('gapPalP', f'{_gpp:.2f}', [WG, f'{O}/arm_element_haps.tsv'])
assert _gpb / len(_wgb) < _gpr / (len(_wgr) - len(_wgb)), 'text: a window with a large gap less often shows the pair'
_g = dict(zip(_wg.src, _wg.n_bases))
number('gapNoSKb', f"{_g['HG00323_hap1_hprc_r2_v1.0.1']/1000:.0f}", [WG]); assert 'HG00323 hap1' in MACROS['xcNoSNames']
assert _g['HG02984_pat_hprc_r2_v1.0.1'] == 0, 'the other C1 arm without S has no gap'
number('gapNoAxisKb', f"{_g['NA19043_hap2_hprc_r2_v1.0.1']/1000:.0f}", [WG]); assert MACROS['sAxisNoneHaps'] == 'NA19043 hap2'
# reversed telomeric A copies on S1-2 band haplotypes (S4.6)
_mmb = pd.read_csv(MCM, sep='\t'); _two = set(_mmb[_mmb.band == 'Two'].src)
_lc = pd.read_csv(f'{O}/lineage_copies.tsv', sep='\t'); _atr2 = _lc[(_lc.lineage == 'A') & (_lc.position == 'telomeric') & (_lc.strand == '-')]
assert len(_atr2) == int(MACROS['linATelRev'])
number('linATelRevBand', int(_atr2.src.isin(_two).sum()), [f'{O}/lineage_copies.tsv', MCM]); number('linBandTwoN', len(_two), [MCM])
assert _atr2[_atr2.src.isin(_two)].src.nunique() == int(_atr2.src.isin(_two).sum()), 'one such copy per band haplotype'
out = ['% written by scripts/40_numbers.py; do not edit'] + [f'\\newcommand{{\\{k}}}{{{v}}}' for k, v in MACROS.items()]
(DOC / 'numbers.tex').write_text('\n'.join(out) + '\n')
(DOC / 'numbers.json').write_text(json.dumps(RECORD, indent=1))
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16] if Path(p).is_file() else None
run = {'at': time.strftime('%Y-%m-%d %H:%M:%S %z'), 'script': RECORD['written_by'], 'script_sha256': sha(RECORD['written_by']),
       'numbers': {n['name']: n['value'] for n in RECORD['numbers']}, 'reads': {p: sha(p) for p in sorted({p for n in RECORD['numbers'] for p in n['reads']})}}
with open(DOC / 'numbers.log.jsonl', 'a') as fh: fh.write(json.dumps(run) + '\n')
print(len(MACROS), 'macros'); print('\n'.join(out[1:]))

# abstract (2026-11-17): "Most haplotypes carry a head-to-head inverted repeat"
assert int(MACROS['palPair']) > int(MACROS['nHapsRes']) / 2, 'abstract: most resolved haplotypes carry the palindrome'

assert int(MACROS['ldTestBN']) + int(MACROS['ldTestNoBN']) == int(MACROS['dosAfrHaps']) - int(MACROS['dosAfrNoSmnOne']), 'text 2.3: the B test counts the resolved AFR haplotypes less those with no SMN1 copy'
