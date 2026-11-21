"""Pipeline step 29: which copy groups form clades in the gene-body ML tree (step 28), with ultrafast-bootstrap
support. The tree is rooted on the orangutan copy. For each group: the smallest clade holding all its copies, that
clade's size, how many copies in it are outside the group, and its support; and the largest clade made only of the
group's copies, with its size and support. Groups: lineage A, lineage B (step 05), SMN1 (A+B), SMN2, h21, h12, the
four apes, and Chen et al.'s SMN1 haplogroups with at least 10 copies (step 11).
Writes results/sv/mltree_clades.tsv."""
import pandas as pd, re, sys
O = 'results/sv'

def parse(s):
    """Newick -> nested dicts {'name','sup','len','kids'}."""
    i = 0
    def node():
        nonlocal i
        n = {'kids': [], 'name': '', 'sup': None, 'len': 0.0}
        if s[i] == '(':
            i += 1
            while True:
                n['kids'].append(node())
                if s[i] == ',': i += 1; continue
                if s[i] == ')': i += 1; break
        m = re.match(r"[^:,();]*", s[i:]); lab = m.group(0); i += len(lab)
        if n['kids']:
            n['sup'] = float(lab) if lab else None
        else: n['name'] = lab
        if i < len(s) and s[i] == ':':
            m = re.match(r":[0-9.eE+-]+", s[i:]); n['len'] = float(m.group(0)[1:]); i += len(m.group(0))
        return n
    return node()

sys.setrecursionlimit(100000)
tree = parse(open(f'{O}/mltree.treefile').read().strip().rstrip(';'))
# reroot on og_orangutan: the treefile is written with the outgroup as a child of the top node
top = [k['name'] for k in tree['kids'] if not k['kids']]
assert 'og_orangutan' in top, top
nodes = []
def walk(n):
    n['tips'] = {n['name']} if not n['kids'] else set().union(*[walk(k) for k in n['kids']])
    nodes.append(n); return n['tips']
walk(tree)
all_tips = tree['tips'] - {'og_orangutan'}
clades = [n for n in nodes if n['kids'] and n is not tree]
clades.append({'tips': all_tips, 'sup': None})   # the root's ingroup side

ml = pd.read_csv(f'{O}/marker_lineage.tsv', sep='\t')
tips = pd.read_csv(f'{O}/copy_tips.tsv', sep='\t')
hg = pd.read_csv(f'{O}/chen_hg_copies.tsv', sep='\t')
G = {'lineage A': set(ml.tip[(ml['class'] == 'SMN1') & ~ml.in_marker_subclade]),
     'lineage B': set(ml.tip[(ml['class'] == 'SMN1') & ml.in_marker_subclade]),
     'SMN1': set(ml.tip[ml['class'] == 'SMN1']), 'SMN2': set(ml.tip[ml['class'] == 'SMN2']),
     'h21': set(ml.tip[ml['class'] == 'hyb_2_1']), 'h12': set(ml.tip[ml['class'] == 'hyb_1_2']),
     'apes (bonobo, chimp, gorilla)': {'og_bonobo', 'og_chimp', 'og_gorilla'},
     'human copies': set(ml.tip)}
hcol = [c for c in hg.columns if c in ('haplogroup', 'hg', 'best')]
assert hcol, hg.columns
for h, sub in hg.groupby(hcol[0]):
    s = set(sub.tip) if 'tip' in sub.columns else set(tips.set_index('copy_id').loc[sub.copy_id].tip)
    if len(s) >= 10 and isinstance(h, str) and h.startswith('S1-'): G[f'haplogroup {h}'] = s
# Also each haplogroup FAMILY that spans more than one label (e.g. S1-6 with S1-6c), as the text counts by family
if 'family' in hg.columns:
    for f, sub in hg.groupby('family'):
        if isinstance(f, str) and f.startswith('S1-') and sub[hcol[0]].nunique() > 1:
            s = set(sub.tip) if 'tip' in sub.columns else set(tips.set_index('copy_id').loc[sub.copy_id].tip)
            if len(s) >= 10: G[f'haplogroup family {f} ({"/".join(sorted(sub[hcol[0]].unique()))})'] = s
rows = []
for g, s in G.items():
    s = s & all_tips; assert s, g
    mrca = min((c for c in clades if s <= c['tips']), key=lambda c: len(c['tips']))
    pure = [c for c in clades if c['tips'] <= s]
    big = max(pure, key=lambda c: (len(c['tips']), c['sup'] or 0)) if pure else None
    rows.append(dict(group=g, copies=len(s), mrca_clade_size=len(mrca['tips']), mrca_outside=len(mrca['tips'] - s),
                     mrca_support=mrca['sup'], monophyletic=len(mrca['tips']) == len(s),
                     largest_pure_clade=len(big['tips']) if big else 0, largest_pure_support=big['sup'] if big else None))
R = pd.DataFrame(rows); R.to_csv(f'{O}/mltree_clades.tsv', sep='\t', index=False)
print(R.to_string())

# backbone: every clade of 15+ copies down to depth 7 below the root, with its make-up and support
import collections
mlx = ml.set_index('tip'); hgx = hg.set_index('tip'); cl = pd.read_csv(f'{O}/copy_labels.tsv', sep='\t').set_index('tip')
def lab(x):
    if x not in mlx.index: return x
    r = mlx.loc[x]
    if r['class'] == 'SMN1': return 'B' if r.in_marker_subclade else 'A'
    return {'SMN2': 'SMN2', 'hyb_2_1': 'h21', 'hyb_1_2': 'h12'}[r['class']]
def fmt(c): return ';'.join(f'{k}:{v}' for k, v in sorted(c.items(), key=lambda kv: (-kv[1], str(kv[0]))))
BR = []
def back(n, d, path):
    t = n['tips']
    BR.append(dict(node=path, depth=d, copies=len(t), support=n['sup'], make_up=fmt(collections.Counter(lab(x) for x in t)),
                   smn1_haplogroups=fmt(collections.Counter(('tie' if hgx.family.get(x) == 'tie' else hgx.hg.get(x)) for x in t if lab(x) in 'AB')),
                   ancestry=fmt(collections.Counter(str(cl.superpop.get(x)) for x in t if x in cl.index))))
    if d < 7:
        ks = sorted(n['kids'], key=lambda k: (-len(k['tips']), min(k['tips'])))
        for i, k in enumerate(ks):
            if len(k['tips']) >= 15 or any(x.startswith('og_') for x in k['tips']): back(k, d + 1, f'{path}.{i}')
back(tree, 0, 'r')
pd.DataFrame(BR).to_csv(f'{O}/mltree_backbone.tsv', sep='\t', index=False)
print(pd.DataFrame(BR)[['node', 'copies', 'support', 'make_up']].to_string())
