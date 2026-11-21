"""Pipeline step 05: SMN1 lineages A and B, cut from the gene-body maximum-likelihood tree (step 04a), as the
researcher decided on 2026-11-14 ("(b) redefine the lineages from the ML halves, so that S1-8 goes with B").
The tree is rooted on the orangutan copy; the human copies form one clade, whose two child clades are the two halves.
Lineage B = the SMN1 copies of the half that holds the c.*3+80T>G copies (results/copies/markers.tsv); lineage A = the
SMN1 copies of the other half. Columns: in_marker_side = the copy is in B's half (SMN1 or h21);
in_marker_subclade = the copy is an SMN1 copy of lineage B (name kept for the steps that read it).
Sub-lineages of B (researcher, 2026-11-14: "maybe we need to subdivide B to B and B'"): the clades inside B's half with
ultrafast bootstrap support >= 95 that lie in no larger such clade are exactly two; B1 = the SMN1 copies of the one that
holds h21 copies, B2 = the SMN1 copies of the other; B copies in neither are 'B-' (column sublineage: A, B1, B2 or B-).
(Until 2026-11-19 B1 was found from a UPGMA lineage; the researcher asked that every sequence tree be the ML tree, and the
rule above gives the same B1 and B2.)
Writes results/sv/marker_lineage.tsv."""
import numpy as np, pandas as pd, collections, sys
sys.path.insert(0, 'scripts'); from smnlib import parse_newick
O = 'results/sv'
t = pd.read_csv(f'{O}/copy_tips.tsv', sep='\t')
t = t[t.copy_class != 'del7_8_SMN2like'].reset_index(drop=True)
mrk = pd.read_csv('results/copies/markers.tsv', sep='\t'); alt = set(mrk.copy_id[mrk.c3p80 == 'alt'])
is_mk = ((t.copy_class == 'SMN1') & t.copy_id.isin(alt)).to_numpy()
# ---- the ML tree: halves of the human copies
sys.setrecursionlimit(100000)
tree = parse_newick(open(f'{O}/mltree.treefile').read().strip().rstrip(';'))
def tipset(n):
    n['tips'] = {n['name']} if not n['kids'] else set().union(*[tipset(k) for k in n['kids']]); return n['tips']
tipset(tree)
human = set(t.tip)
def find(n):
    if n['tips'] == human: return n
    for k in n['kids']:
        if human <= k['tips']: return find(k)
    raise SystemExit('human copies are not one clade')
H = find(tree)
assert len(H['kids']) == 2, 'the human clade splits in two'
halves = sorted(H['kids'], key=lambda k: -len(set(t.tip[is_mk]) & k['tips']))
Bh, Ah = halves
assert set(t.tip[is_mk]) <= Bh['tips'], 'every marker copy is in one half'
assert Bh['sup'] >= 90 and Ah['sup'] >= 90, (Bh['sup'], Ah['sup'])
inB = t.tip.isin(Bh['tips']).to_numpy()
cb = collections.Counter(t.copy_class[inB])
assert set(cb) == {'SMN1', 'hyb_2_1'}, cb   # no SMN2 or h12 in B's half
linB = inB & (t.copy_class == 'SMN1').to_numpy()
# sub-lineages of B
def clades(n, out):
    if n['kids']: out.append(n); [clades(k, out) for k in n['kids']]
    return out
inner = [c for c in clades(Bh, []) if c is not Bh and (c['sup'] or 0) >= 95]
maxi = [c for c in inner if not any(c['tips'] < o['tips'] for o in inner)]
assert len(maxi) == 2, f'{len(maxi)} maximal supported clades in B half'
hyb = set(t.tip[t.copy_class == 'hyb_2_1'])
c1 = max(maxi, key=lambda c: len(c['tips'] & hyb)); c2 = next(c for c in maxi if c is not c1)
assert c1['tips'] & hyb and not (c2['tips'] & hyb), 'B1 holds the h21 copies, B2 none'
assert c1['sup'] >= 90 and c2['sup'] >= 90, (c1['sup'], c2['sup'])
assert set(t.copy_class[t.tip.isin(c2['tips'])]) == {'SMN1'}, 'B2 clade is SMN1 only'
smn1 = (t.copy_class == 'SMN1').to_numpy()
subl = np.where(~smn1, '', np.where(~linB, 'A', np.where(t.tip.isin(c1['tips']), 'B1', np.where(t.tip.isin(c2['tips']), 'B2', 'B-'))))
print('B1', int((subl == 'B1').sum()), 'support', c1['sup'], '; B2', int((subl == 'B2').sum()), 'support', c2['sup'], '; B-', int((subl == 'B-').sum()))
pd.DataFrame({'tip': t.tip, 'copy_id': t.copy_id, 'class': t.copy_class, 'c3p80_alt': is_mk, 'in_marker_side': inB,
              'in_marker_subclade': linB, 'sublineage': subl}).to_csv(f'{O}/marker_lineage.tsv', sep='\t', index=False)
print('ML halves: B side', len(Bh['tips']), 'support', Bh['sup'], '; A side', len(Ah['tips']), 'support', Ah['sup'])
print('lineage B', int(linB.sum()), 'SMN1 copies; marker copies', int(is_mk.sum()), 'all in B')
