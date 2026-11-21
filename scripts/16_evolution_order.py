# Pipeline step 16 (evolution_order); ported from explore/structure_vs_sequence/11_evolution_order.py (sha256 f2ac82dbce84, commit 5e2b3bb),
# paths rewritten to results/sv/ and figures/; run from the project root after the steps before it.
"""Exploration step 11: the proposed order of events for SMN1 lineage A, lineage B and SMN2, drawn with the
evidence each event rests on.

Recomputes, from out/copies_with_outgroups.aln.fa (step 09's alignment with four great-ape SMN gene bodies):
  * S sites: gene-body columns where SMN1 (A+B) and SMN2 are each >=90% fixed for different bases;
  * AB sites: gene-body columns where A and B are each >=90% fixed for different bases (as step 09).
The ancestral base at a site is the base carried by at least 3 of the 4 apes (else unresolved).
Reads results/sv/evolution_sites.tsv (step 13b) and ape_site_reads.tsv (step 13c); writes evolution_sites_checked.tsv. Writes results/sv/evolution_summary.json, figures/evolution_order.png/.pdf.
Literature statements in the figure are quoted (Rochette 2001, Suzuki 2026 abstracts) and marked as such.
Ape copy numbers come from step 12 (out/ape_copies.tsv, primary haplotypes) and step 15 (out/ape_copies_hap2.tsv,
other haplotypes; out/ape_smn_gene_records.tsv, NCBI Gene records): event 1 is drawn as the human-ape split, with the
timing of the SMN duplication left open (one SMN copy in each of the 8 ape haplotypes searched).
2026-11-14: lineages A/B are the ML-tree halves (step 05). The order drawn follows the ML tree (B's half against the
half holding A and SMN2) and the site evidence; the evidence against it is drawn beside it."""
import json, numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
O = 'results/sv'; F = 'figures'
THR = 0.9
names, seqs = [], []
for l in open(f'{O}/copies_with_outgroups.aln.fa'):
    l = l.strip()
    if l[0] == '>': names.append(l[1:].split()[0]); seqs.append([])
    else: seqs[-1].append(l.upper())
M = np.array([np.frombuffer(''.join(s).encode(), dtype=np.uint8) for s in seqs]); del seqs
idx = {n: i for i, n in enumerate(names)}
t = pd.read_csv(f'{O}/copy_tips.tsv', sep='\t')
ml = pd.read_csv(f'{O}/marker_lineage.tsv', sep='\t').set_index('tip')
lab = pd.read_csv(f'{O}/copy_labels.tsv', sep='\t').set_index('tip')
ref = idx[t[(t.copy_class == 'SMN1') & t.src.str.startswith('GCA_000001405')].tip.iloc[0]]
rpos = np.cumsum(M[ref] != ord('-')) - 1; L = rpos[-1] + 1
gene = (rpos >= 2000) & (rpos < L - 2000) & (M[ref] != ord('-'))
grp = {'A': [], 'B': [], 'SMN2': []}
for r in t.itertuples():
    if r.copy_class == 'SMN1' and r.tip in ml.index: grp['B' if ml.loc[r.tip, 'in_marker_subclade'] else 'A'].append(idx[r.tip])
    elif r.copy_class == 'SMN2': grp['SMN2'].append(idx[r.tip])
grp['SMN1'] = grp['A'] + grp['B']
APES = ['og_chimp', 'og_bonobo', 'og_gorilla', 'og_orangutan']
APE_LABEL = {'og_chimp': 'chimpanzee', 'og_bonobo': 'bonobo', 'og_gorilla': 'gorilla', 'og_orangutan': 'orangutan'}
assert all(a in idx for a in APES)
def maj(rows, c):
    col = M[rows, c]; col = col[np.isin(col, list(b'ACGT'))]
    if len(col) < 0.8 * len(rows): return None, 0.0
    v, n = np.unique(col, return_counts=True); k = n.argmax(); return chr(v[k]), n[k] / len(col)
def ancestral(c):
    b = [chr(M[idx[a], c]) for a in APES]; b = [x for x in b if x in 'ACGT']
    if not b: return None, 0
    v, n = np.unique(b, return_counts=True); k = n.argmax()
    return (str(v[k]), int(n[k])) if n[k] >= 3 else (None, int(n[k]))
S = pd.read_csv(f'{O}/evolution_sites.tsv', sep='\t')   # step 13b
# 2026-11-15: the ape bases are the read-checked ones of step 13c (raw primary-assembly sequence), not the multiple
# alignment's: an accepted read is used; a read that is not accepted keeps the alignment's base only when it agrees with
# it; otherwise the ape's base at that site is unresolved ('?').
_rd = pd.read_csv(f'{O}/ape_site_reads.tsv', sep='\t', keep_default_na=False)
_fin = {}
for r in _rd.itertuples():
    _fin[(r.kind, int(r.grch38), r.ape)] = r.read if r.accepted else (r.msa if r.read == r.msa else '?')
for ap in APES:
    S[ap] = [_fin[(k, int(g), ap)] for k, g in zip(S.kind, S.grch38)]
_anc = []
for _, r in S.iterrows():
    b = [r[a] for a in APES if r[a] in ('A', 'C', 'G', 'T')]
    v, n = np.unique(b, return_counts=True) if b else (np.array([]), np.array([0]))
    k = int(n.argmax()) if b else 0
    _anc.append((str(v[k]), int(n[k])) if b and n[k] >= 3 else (None, int(n[k]) if b else 0))
S['ancestral'] = [a for a, _ in _anc]; S['apes_agreeing'] = [n for _, n in _anc]
S.to_csv(f'{O}/evolution_sites_checked.tsv', sep='\t', index=False)
s = S[S.kind == 'S']; ab = S[S.kind == 'AB']
C840, CS239 = 70951946, 70952674
assert {C840, CS239} <= set(s.grch38), 'c.840 / c.*239 not among the fixed SMN1/SMN2 sites'
def ape_state(df, ap, which):
    return int((df[ap] == df[which]).sum())
X = {}
X['n_S'] = len(s); X['n_AB'] = len(ab)
X['S_span'] = [int(s.grch38.min()), int(s.grch38.max())]
X['S_anc_smn1'] = int((s.ancestral == s.SMN1).sum()); X['S_anc_smn2'] = int((s.ancestral == s.SMN2).sum())
X['S_unres'] = int(s.ancestral.isna().sum())
X['S_AB_both_smn1'] = int(((s.A == s.SMN1) & (s.B == s.SMN1)).sum())
X['S_shared_derived_AB'] = int(((s.ancestral == s.SMN2) & (s.A == s.SMN1) & (s.B == s.SMN1)).sum())
X['S_apes_unanimous'] = int((s.apes_agreeing == 4).sum())
X['c840_apes'] = ''.join(s[s.grch38 == C840][APES].iloc[0]); X['cs239_apes'] = ''.join(s[s.grch38 == CS239][APES].iloc[0])
X['c840'] = s[s.grch38 == C840][['SMN1', 'SMN2']].iloc[0].tolist(); X['cs239'] = s[s.grch38 == CS239][['SMN1', 'SMN2']].iloc[0].tolist()
X['S_apes_smn2_base'] = {APE_LABEL[a]: ape_state(s, a, 'SMN2') for a in APES}
X['AB_span'] = [int(ab.grch38.min()), int(ab.grch38.max())]
X['AB_anc_A'] = int((ab.ancestral == ab.A).sum()); X['AB_anc_B'] = int((ab.ancestral == ab.B).sum())
X['AB_unres'] = int(ab.ancestral.isna().sum())
X['AB_smn2_eq_A'] = int((ab.SMN2 == ab.A).sum())
X['AB_smn2_eq_A_anc'] = int(((ab.SMN2 == ab.A) & (ab.ancestral == ab.A)).sum())
X['AB_shared_derived_A_SMN2'] = int(((ab.SMN2 == ab.A) & (ab.ancestral == ab.B)).sum())
X['AB_apes_A_base'] = {APE_LABEL[a]: ape_state(ab, a, 'A') for a in APES}
g = pd.read_csv(f'{O}/group_pair_distances.tsv', sep='\t').set_index(['a', 'b'])
for k, (p, q) in {'d_AA': ('A', 'A'), 'd_BB': ('B', 'B'), 'd_22': ('SMN2', 'SMN2'), 'd_AB': ('A', 'B'),
                  'd_A2': ('A', 'SMN2'), 'd_B2': ('B', 'SMN2')}.items(): X[k] = float(g.loc[(p, q), 'mean_per_kb'])
btips = ml.index[ml.in_marker_subclade.astype(bool)]
X['nA'], X['nB'], X['nSMN2'] = len(grp['A']), len(grp['B']), len(grp['SMN2'])
X['B_afr'] = int(lab.loc[btips, 'afr'].astype(bool).sum())
X['n_h12'] = int((lab.cls == 'h12').sum()); X['n_h21'] = int((lab.cls == 'h21').sum())
P1 = pd.read_csv(f'{O}/ape_copies.tsv', sep='\t'); P1 = P1[P1.ape != 'human GRCh38']
P2 = pd.read_csv(f'{O}/ape_copies_hap2.tsv', sep='\t'); P2 = P2[P2.rule == 'strict']
GR = pd.read_csv(f'{O}/ape_smn_gene_records.tsv', sep='\t'); GR = GR[GR.live & GR.description.str.startswith('survival') & ~GR.symbol.str.startswith('SMNDC')]
X['ape_haps'] = int(len(P1.ape.unique()) + len(P2.ape.unique()))
X['ape_copies'] = int(len(P1) + len(P2))
X['ape_haps_smn1_tests'] = int(((P1.c840 == 'C') & (P1.c239 == 'G')).sum() + ((P2.c840 == 'C') & (P2.c239 == 'G')).sum())
X['gene_full'] = int((GR.length > 20000).sum()); X['gene_short'] = int((GR.length <= 20000).sum())
X['gene_short_range'] = [int(GR[GR.length <= 20000].length.min()), int(GR[GR.length <= 20000].length.max())]
assert X['ape_haps'] == 8 and X['ape_copies'] == 8 and X['ape_haps_smn1_tests'] == 8, 'not one SMN1-type copy per ape haplotype'
assert X['gene_full'] == 4, 'not one full-length SMN gene record per ape'
assert X['S_AB_both_smn1'] == X['n_S'], 'A and B do not both carry the SMN1 base at every S site'
assert X['c840_apes'] == X['c840'][0] * 4 and X['cs239_apes'] == X['cs239'][0] * 4, 'apes not all SMN1-type at the test sites'
# ML-tree halves and their support (step 04a tree, step 05 lineages)
import sys; sys.path.insert(0, 'scripts'); from smnlib import parse_newick; sys.setrecursionlimit(100000)
_tr = parse_newick(open(f'{O}/mltree.treefile').read().strip().rstrip(';'))
def _ts(n):
    n['tips'] = {n['name']} if not n['kids'] else set().union(*[_ts(k) for k in n['kids']]); return n['tips']
_ts(_tr); _hum = set(ml.index)
def _find(n):
    if n['tips'] == _hum: return n
    return _find(next(k for k in n['kids'] if _hum <= k['tips']))
_H = _find(_tr); _bside = set(ml.index[ml.in_marker_side.astype(bool)])
_hb = next(k for k in _H['kids'] if k['tips'] == _bside); _ha = next(k for k in _H['kids'] if k is not _hb)
X['ml_sup_B_half'], X['ml_sup_A_half'] = _hb['sup'], _ha['sup']
X['ml_SMN2_in_A_half'] = int(len(set(lab.index[lab.cls == 'SMN2']) & _ha['tips']))
assert X['ml_SMN2_in_A_half'] == X['nSMN2'], 'not every SMN2 copy in A\'s half'
assert X['AB_shared_derived_A_SMN2'] == X['AB_anc_B'], 'where B has the old base, A and SMN2 do not always share the new one'
assert X['AB_smn2_eq_A'] == X['n_AB'], 'SMN2 does not carry A\'s base at every A/B site'
assert X['AB_anc_B'] > X['AB_anc_A'], 'B no longer carries the old base at most A/B sites'
assert X['S_shared_derived_AB'] > 0, 'no new SMN1 base shared by A and B'
assert X['d_B2'] > X['d_A2']
json.dump(X, open(f'{O}/evolution_summary.json', 'w'), indent=1)
print(json.dumps(X, indent=1))
print(S.to_string())

# ---------------- figure ----------------
# Lineage colours from smnlib.COPYCOL, as in Figs. 1-5: A red, B black (B1's colour; B2 is brown), SMN2 blue.
# The shared human trunk is drawn dark grey so that it is not read as lineage B.
from smnlib import COPYCOL
CA, CB, C2, CAPE, CT = COPYCOL['A'], COPYCOL['B1'], COPYCOL['SMN2'], '#7f7f7f', '0.45'
fig = plt.figure(figsize=(14, 11.5))
ax = fig.add_axes([0.02, 0.40, 0.40, 0.57]); ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis('off')
ax.text(0, 9.95, 'a  Proposed order of events (a hypothesis)', fontsize=12, weight='bold', va='top')
def ln(x0, y0, x1, y1, c=CT, ls='-', lw=2.5): ax.plot([x0, x1], [y0, y1], color=c, ls=ls, lw=lw, solid_capstyle='round')
def ev(x, y, n): ax.text(x, y, n, ha='center', va='center', fontsize=10, weight='bold', color='white',
                         bbox=dict(boxstyle='circle,pad=0.25', fc='k', ec='k'))
ax.annotate('', xy=(0.25, 0.7), xytext=(0.25, 9.0), arrowprops=dict(arrowstyle='->', lw=1.2, color='0.4'))
ax.text(0.45, 4.8, 'older  →  younger', rotation=-90, va='center', fontsize=9, color='0.4')
ln(5, 9.0, 5, 8.0)
ax.text(5.15, 9.15, 'ancestral SMN gene (one copy?)', fontsize=9.5, va='bottom')
ln(5, 8.0, 2.0, 6.6, CAPE); ln(2.0, 6.6, 2.0, 1.0, CAPE)
ev(5, 8.0, '1')
ax.text(2.0, 0.9, 'great apes:\none SMN copy per\nhaplotype searched\n(SMN1-type at\nc.840, c.*239)', ha='center', va='top', fontsize=8.5, color=CAPE)
ev(2.0, 4.2, '2')
ln(5, 8.0, 6.0, 6.8)
ax.text(5.45, 7.55, 'human\nlineage', fontsize=8.5, ha='left')
ax.text(6.25, 6.95, 'duplication?', fontsize=8, ha='left', style='italic', color='0.3')
ln(6.0, 6.8, 4.0, 5.2, CB); ln(4.0, 5.2, 4.0, 1.0, CB)
ln(6.0, 6.8, 7.2, 5.4, CT); ev(6.0, 6.8, '3')
ln(7.2, 5.4, 6.0, 4.0, CA); ln(6.0, 4.0, 6.0, 1.0, CA)
ln(7.2, 5.4, 8.6, 4.0, C2); ln(8.6, 4.0, 8.6, 1.0, C2)
ev(7.2, 5.4, '4')
ax.text(4.0, 0.9, 'SMN1\nlineage B\n(B1, B2)', color=CB, ha='center', va='top', fontsize=10, weight='bold')
ax.text(6.0, 0.9, 'SMN1\nlineage A', color=CA, ha='center', va='top', fontsize=10, weight='bold')
ax.text(8.6, 0.9, 'SMN2', color=C2, ha='center', va='top', fontsize=10, weight='bold')
for x0, x1, y in [(4.0, 6.0, 2.2), (6.0, 8.6, 1.6)]:
    ax.annotate('', xy=(x1 - 0.08, y), xytext=(x0 + 0.08, y), arrowprops=dict(arrowstyle='<->', ls='--', color='0.35', lw=1.1))
ev(6.0, 2.75, '5')
ex = fig.add_axes([0.43, 0.40, 0.56, 0.57]); ex.axis('off'); ex.set_xlim(0, 1); ex.set_ylim(0, 1)
kb = lambda a, b: f'{(b - a + 1) / 1000:.1f} kb'
apes2 = ', '.join(f'{k} {v}' for k, v in X['S_apes_smn2_base'].items())
apesA = ', '.join(f'{k} {v}' for k, v in X['AB_apes_A_base'].items())
E = [
 ('1', 'Human and ape lineages split; when SMN duplicated is open',
  [f"[here] T2T ape assemblies: {X['ape_copies']} SMN copies in {X['ape_haps']} haplotypes (2 per ape), one each; NCBI Gene: {X['gene_full']} full-length",
   f"   SMN genes (one per ape), {X['gene_short']} other SMN-like records of {X['gene_short_range'][0]:,}-{X['gene_short_range'][1]:,} bp.",
   '[quoted] Rochette 2001 (abstract, qPCR): chimpanzees carry 2-7 SMN copies per diploid genome.']),
 ('2', 'The single ape copy: a mix at the SMN1/SMN2 sites',
  [f"[here] All {X['ape_haps_smn1_tests']} ape copies carry SMN1's base at c.840 ({X['c840'][0]}) and at c.*239 ({X['cs239'][0]}).",
   f"[here] At the {X['n_S']} fixed SMN1/SMN2 sites the primary ape copy has SMN2's base in {apes2}."]),
 ('3', 'Human copies split in two: lineage B against the A + SMN2 side',
  [f"[here] ML tree: B's half (support {X['ml_sup_B_half']:.0f}) against the half holding A and all {X['ml_SMN2_in_A_half']} SMN2 copies (support {X['ml_sup_A_half']:.0f}).",
   f"[here] {X['n_AB']} sites fixed between A and B, chr5:{X['AB_span'][0]:,}–{X['AB_span'][1]:,} (exon 1 - intron 2). Ape base = B's at {X['AB_anc_B']}, A's at {X['AB_anc_A']}" + (f", unresolved {X['AB_unres']}." if X['AB_unres'] else '.'),
   f"[here] At all {X['AB_shared_derived_A_SMN2']} sites where B has the ape base, A and SMN2 share the new base."]),
 ('4', 'SMN2 arises on the A side; diverges in one 2.5-kb stretch',
  [f"[here] All {X['n_S']} fixed SMN1/SMN2 sites lie in chr5:{X['S_span'][0]:,}–{X['S_span'][1]:,} ({kb(*X['S_span'])}, intron 6 – exon 8).",
   f"[here] Ape base = SMN1's at {X['S_anc_smn1']} (c.840, c.*239 among them: SMN2's bases there are new), = SMN2's at {X['S_anc_smn2']}" + (f", unresolved {X['S_unres']}." if X['S_unres'] else '.'),
   '[quoted] Rochette 2001: no non-human primate carries the SMN2 exon 7-8 bases.']),
 ('5', 'Conversion and crossing-over continue',
  [f"[here] Hybrid copies in the tree: {X['n_h12']} h12, {X['n_h21']} h21. Mean distance A-SMN2 {X['d_A2']:.2f}, A-B {X['d_AB']:.2f}, B-SMN2 {X['d_B2']:.2f} per kb.",
   '[quoted] Suzuki 2026 (abstract): putatively biased, one-way gene conversion around SMN.']),
 ('vs', 'Against this order: the intron 6 - exon 8 stretch',
  [f"[here] A and B both carry SMN1's base at all {X['n_S']} SMN1/SMN2 sites, and share the new SMN1 base at {X['S_shared_derived_AB']}: there, A groups with B.",
   f"[here] Upstream {X['AB_shared_derived_A_SMN2']} new bases group A with SMN2; here {X['S_shared_derived_AB']} group A with B: different parts of the gene",
   '   tell different histories, as conversion between the copies would make them.']),
]
y = 0.995
for n, head, lines in E:
    ex.text(0.0, y, n if n != 'vs' else '≠', fontsize=10, weight='bold', color='white', va='top',
            bbox=dict(boxstyle='circle,pad=0.25', fc='k' if n != 'vs' else '0.45', ec='none'))
    ex.text(0.045, y, head, fontsize=10, weight='bold', va='top')
    y -= 0.045
    for s_ in lines:
        ex.text(0.045, y, s_, fontsize=8.3, va='top'); y -= 0.037
    y -= 0.018
bx = fig.add_axes([0.10, 0.155, 0.87, 0.19])
cols = pd.concat([s, ab]).reset_index(drop=True)
rowsdef = [('lineage A (%d)' % X['nA'], 'A'), ('lineage B (%d)' % X['nB'], 'B'), ('SMN2 (%d)' % X['nSMN2'], 'SMN2')] + [(APE_LABEL[a], a) for a in APES]
# Site states in panel b: teal / ochre, so that no state colour is a lineage colour.
C1, C2b, CG = '#5ab4ac', '#d8b365', '#dddddd'
for j, r in cols.iterrows():
    for i, (lab_, key) in enumerate(rowsdef):
        base = r[key]
        col = C1 if base == r.state1 else C2b if base == r.state2 else CG
        bx.add_patch(plt.Rectangle((j, i), 0.92, 0.88, fc=col, ec='none'))
        if i < 3 and isinstance(r.ancestral, str) and base in 'ACGT' and base != r.ancestral:
            bx.plot(j + 0.46, i + 0.44, 'o', ms=5, mfc='k', mec='white', mew=0.6)
        if i >= 3: bx.text(j + 0.46, i + 0.44, base if base in 'ACGT' else ('?' if base == '?' else '–'), ha='center', va='center', fontsize=6.5)
bx.set_xlim(-0.2, len(cols)); bx.set_ylim(len(rowsdef) + 0.1, -0.2)
bx.set_yticks([i + 0.44 for i in range(len(rowsdef))]); bx.set_yticklabels([r[0] for r in rowsdef], fontsize=8.5)
lbl = []
for r in cols.itertuples():
    tag = '\nc.840' if r.grch38 == C840 else '\nc.*239' if r.grch38 == CS239 else ''
    lbl.append(f'{r.grch38:,}{tag}')
bx.set_xticks([j + 0.46 for j in range(len(cols))]); bx.set_xticklabels(lbl, rotation=90, fontsize=7, ha='center')
bx.set_xlabel('site (GRCh38 chr5 position)', fontsize=9)
for sp in bx.spines.values(): sp.set_visible(False)
bx.tick_params(length=0)
ns = len(s)
bx.axvline(ns - 0.04, color='k', lw=1)
bx.text(ns / 2, -0.35, f'{ns} sites fixed between SMN1 and SMN2 (event 4): teal = SMN1 base, ochre = SMN2 base', ha='center', va='bottom', fontsize=8.5)
bx.text(ns + len(ab) / 2, -0.35, f'{len(ab)} sites fixed between A and B (event 3): teal = A base, ochre = B base', ha='center', va='bottom', fontsize=8.5)
fig.text(0.02, 0.375, 'b  Evidence at each site: majority base of each human group, and the base of each ape copy', fontsize=12, weight='bold')
fig.legend(handles=[Patch(fc=C1, label='first state (SMN1 / A)'), Patch(fc=C2b, label='second state (SMN2 / B)'),
                   Patch(fc=CG, label='other base or gap; ? = ape base unresolved in its raw sequence'),
                   Line2D([], [], marker='o', ls='', mfc='k', mec='white', label='new (derived) base: differs from the base of ≥3 apes')],
          loc='lower center', bbox_to_anchor=(0.5, 0.005), ncol=4, fontsize=8.5, frameon=False)
for ext in ('png', 'pdf'):
    fig.savefig(f'{F}/evolution_order.{ext}', dpi=200, **({'metadata': {'CreationDate': None}} if ext == 'pdf' else {}))
print('wrote', f'{F}/evolution_order.png')
