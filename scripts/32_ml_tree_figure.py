"""Pipeline step 32: the gene-body ML tree (step 28) drawn as a phylogram, rooted on the orangutan copy, with strips
for each copy: class and SMN1 lineage (step 05), Chen et al.'s SMN1 haplogroup (step 11), the arm of the head-to-head R_A pair that holds
it, as in Fig. 1 (step 22) and the donor's continental group. Ultrafast-bootstrap support is printed on the
backbone clades of 15 or more copies (results/sv/mltree_backbone.tsv, step 29). Writes figures/ml_tree.png/.pdf."""
import sys, collections, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, 'scripts'); from smnlib import parse_newick
sys.setrecursionlimit(100000)
O = 'results/sv'
tree = parse_newick(open(f'{O}/mltree.treefile').read().strip().rstrip(';'))
ml = pd.read_csv(f'{O}/marker_lineage.tsv', sep='\t').set_index('tip')
hg = pd.read_csv(f'{O}/chen_hg_copies.tsv', sep='\t').set_index('tip')
cl = pd.read_csv(f'{O}/copy_labels.tsv', sep='\t').set_index('tip')
# arm of the head-to-head R_A pair that holds the copy, as in Fig. 1 (results/sv/palindrome_copies_lineage.tsv, step 22)
_pl = pd.read_csv(f'{O}/palindrome_copies_lineage.tsv', sep='\t', keep_default_na=False)
ARM = {c: ('cenSMN arm' if a.startswith('left') else 'telSMN arm') for c, a in zip(_pl.copy_id, _pl.arm)}
assert set(_pl.arm.str.split(' ').str[0]) == {'left', 'right'}
def slot(t):
    if t not in cl.index or pd.isna(cl.position.get(t)): return None  # reference rows and unresolved haplotypes: not assessed
    return ARM.get(cl.copy_id[t], 'no R_A pair')
# arrangement type of the copy's haplotype (Fig. 1, Fig. 2 groups); types on fewer than 10 resolved haplotypes pooled
_hc = pd.read_csv('tables/haplotype_copy_counts.tsv', sep='\t')
_vc = _hc[_hc.resolved].arrangement.value_counts()
TOPARR = [a for a, n in _vc.items() if n >= 10]
_pal = ['#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02', '#a6761d', '#1f78b4']
assert len(TOPARR) <= len(_pal)
ARRC = dict(zip(TOPARR, _pal)); ARRC['other types'] = '#d9d9d9'
def arr(t):
    if t not in cl.index or pd.isna(cl.position.get(t)): return None
    a = cl.arrangement[t]
    return a if a in ARRC else 'other types'
# ladderize, then place leaves
def size(n): n['n'] = 1 if not n['kids'] else sum(size(k) for k in n['kids']); return n['n']
size(tree)
leaves = []
def place(n, x):
    n['x'] = x + n['len'] if n is not tree else 0.0
    if not n['kids']: n['y'] = len(leaves); leaves.append(n); return n['y']
    ys = [place(k, n['x']) for k in sorted(n['kids'], key=lambda k: k['n'])]
    n['y'] = (min(ys) + max(ys)) / 2; return n['y']
place(tree, 0.0)
N = len(leaves)
XH = max(l['x'] for l in leaves if not l['name'].startswith('og_'))   # deepest human copy
GAP = XH * 0.25
def X(x):   # the ape branches run far past the human copies: beyond XH they are compressed 50-fold (broken axis)
    return x if x <= XH else XH + min(GAP, (x - XH) / 50)
fig = plt.figure(figsize=(8.5, 11)); ax = fig.add_axes([0.04, 0.1, 0.5, 0.82])
def draw(n):
    for k in n['kids']:
        ax.plot([X(n['x']), X(n['x'])], [n['y'], k['y']], color='k', lw=.35)
        ax.plot([X(n['x']), X(k['x'])], [k['y'], k['y']], color='k', lw=.35)
        draw(k)
    if n['kids'] and n['n'] >= 40 and n['sup'] is not None and n is not tree:
        ax.text(X(n['x']) - XH * .004, n['y'] + 2, f"{int(n['sup'])}", fontsize=5.5, ha='right', va='bottom', color='#b2182b')
draw(tree)
apes = sorted([l for l in leaves if l['name'].startswith('og_')], key=lambda l: l['y'])
for k, l in enumerate(apes):
    ax.annotate(l['name'][3:].replace('chimp', 'chimpanzee'), xy=(X(l['x']), l['y']), xytext=(XH * 1.08, l['y'] + 9 + 11 * k), zorder=5,
                fontsize=6, va='center', arrowprops=dict(arrowstyle='-', lw=.3, color='#777'))
ax.axvspan(XH * 1.01, XH * 1.03, color='white', zorder=3); ax.text(XH * 1.02, -1, '//', ha='center', va='center', fontsize=8, zorder=4)
ax.set_ylim(-2, N + 2); ax.set_yticks([]); ax.set_xlim(0, XH + GAP * 1.6)
ax.set_xticks([x for x in np.arange(0, XH, 0.001)])
for sp in ('left', 'right', 'top'): ax.spines[sp].set_visible(False)
ax.set_xlabel('substitutions per site (GTR+F+G4)', fontsize=8)
def cls(t):
    if t not in ml.index: return 'ape'
    r = ml.loc[t]
    if r['class'] == 'SMN1': return {'A': 'SMN1 A', 'B1': 'SMN1 B1', 'B2': 'SMN1 B2', 'B-': 'SMN1 B other'}[r.sublineage]
    return {'SMN2': 'SMN2', 'hyb_2_1': 'h21', 'hyb_1_2': 'h12'}[r['class']]
HGF = {'S1-1': 'S1-1/1c', 'S1-1c': 'S1-1/1c', 'S1-2': 'S1-2', 'S1-2c': 'S1-2', 'S1-2|S1-2c': 'S1-2', 'S1-3': 'S1-3', 'S1-4': 'S1-4', 'S1-4c': 'S1-4', 'S1-5': 'S1-5',
       'S1-6': 'S1-6', 'S1-6c': 'S1-6', 'S1-8': 'S1-8', 'S1-9': 'S1-9 family', 'S1-9c': 'S1-9 family', 'S1-9d': 'S1-9 family', 'S1-10': 'S1-10'}
def hgf(t):
    if t not in hg.index: return None
    if hg.family.get(t) == 'tie': return 'tie'
    return HGF[hg.hg[t]]  # every non-tie label is mapped
S = [('class', cls, {'SMN1 A': '#e41a1c', 'SMN1 B1': '#000000', 'SMN1 B2': '#a65628', 'SMN1 B other': '#666666', 'SMN2': '#377eb8', 'h21': '#984ea3', 'h12': '#c2a5cf', 'ape': '#ffffff'}),
     ('haplogroup', hgf, {'S1-1/1c': '#a6cee3', 'S1-2': '#1f78b4', 'S1-3': '#b2df8a', 'S1-4': '#fdbf6f', 'S1-5': '#33a02c', 'S1-6': '#fb9a99',
                                'S1-8': '#e31a1c', 'S1-9 family': '#3d3d3d', 'S1-10': '#ff7f00', 'tie': '#cab2d6'}),
     ('arm (Fig. 1)', slot, {'cenSMN arm': '#8073ac', 'telSMN arm': '#e08214', 'no R_A pair': '#bdbdbd'}),
     ('arrangement', arr, ARRC),
     ('ancestry', lambda t: cl.superpop.get(t) if t in cl.index else None, {'AFR': '#8c510a', 'AMR': '#d73027', 'EAS': '#4575b4', 'EUR': '#1a9850', 'SAS': '#984ea3', 'other': '#bdbdbd'})]
for j, (title, f, cmap) in enumerate(S):
    a = fig.add_axes([0.555 + 0.04 * j, 0.1, 0.035, 0.82], sharey=ax)
    for l in leaves:
        v = f(l['name']); c = cmap.get(v) if v is not None else None
        if c: a.add_patch(plt.Rectangle((0, l['y'] - .5), 1, 1, color=c, lw=0))
    a.set_xlim(0, 1); a.set_xticks([]); a.set_yticks([]); a.set_title(title, fontsize=7, rotation=90, ha='center', va='bottom')
    for sp in a.spines.values(): sp.set_visible(False)
handles = []
for title, f, cmap in S:
    cnt = collections.Counter(f(l['name']) for l in leaves)
    handles += [plt.Rectangle((0, 0), 1, 1, color=c, ec='#999', lw=.3, label=f'{title.replace("haplogroup", "SMN1 haplogroup")}: {v} ({cnt.get(v, 0)})') for v, c in cmap.items() if v != 'ape']
fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(0.765, 0.96), fontsize=5.5, frameon=False, ncol=1)
fig.text(0.775, 0.30, 'Red numbers: ultrafast-bootstrap\nsupport (%) of clades with\n40 or more copies.\nRooted on the orangutan copy;\nape branches compressed\n50-fold past the break (//).\nWhite strip cells: not assessed.', fontsize=6.5, va='top')
fig.savefig('figures/ml_tree.png', dpi=220); fig.savefig('figures/ml_tree.pdf', metadata={'CreationDate': None})
print(N, 'leaves')
