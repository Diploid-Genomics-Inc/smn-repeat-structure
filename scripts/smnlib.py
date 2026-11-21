"""Small helpers shared by the pipeline steps.
read_genbank: a minimal GenBank flat-file reader for data/refs_seq/smn1_smn2_refseq.gb (NCBI RefSeq records
NM_000344.4 and NM_017411.4). Returns {VERSION: {'seq': str (upper case), 'exons': [(start, end), ...] 1-based
inclusive in transcript order, 'cds': (start, end)}}. Written for this project (the clinical study's smnlib.py, which
exploration step 20 imported, is not in this project); it reads only LOCUS/VERSION, exon, CDS and ORIGIN lines."""
import re

def read_genbank(path):
    out, cur, in_seq = {}, None, False
    for line in open(path):
        if line.startswith('VERSION'):
            cur = {'seq': [], 'exons': [], 'cds': None}; out[line.split()[1]] = cur; in_seq = False
        elif line.startswith('ORIGIN'):
            in_seq = True
        elif line.startswith('//'):
            in_seq = False
        elif in_seq and cur is not None:
            cur['seq'].append(re.sub(r'[^a-zA-Z]', '', line))
        elif cur is not None:
            m = re.match(r'^     (exon|CDS)\s+(\d+)\.\.(\d+)\s*$', line)
            if m:
                a, b = int(m.group(2)), int(m.group(3))
                if m.group(1) == 'exon': cur['exons'].append((a, b))
                elif cur['cds'] is None: cur['cds'] = (a, b)
    for v in out.values():
        v['seq'] = ''.join(v['seq']).upper()
        assert v['cds'] and v['exons'], 'record without CDS or exons'
    return out


def parse_newick(s):
    """Newick string (no trailing ';') -> nested dicts {'name','sup','len','kids'}; internal labels read as support."""
    import re as _re
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
        m = _re.match(r"[^:,();]*", s[i:]); lab = m.group(0); i += len(lab)
        if n['kids']: n['sup'] = float(lab) if lab else None
        else: n['name'] = lab
        if i < len(s) and s[i] == ':':
            m = _re.match(r":[0-9.eE+-]+", s[i:]); n['len'] = float(m.group(0)[1:]); i += len(m.group(0))
        return n
    return node()


# One colour per repeat-unit family, shared by every figure that draws units (Figs. 1, 2, 4), so a family has the same
# colour everywhere. No family colour is a copy-tick colour: D was purple, the hybrid-copy colour, and is now tan.
# G was green (#27ae60), too close to B's teal; it is now violet (2026-11-14, at the researcher's request).
FAMCOL = {'A': '#e67e22', 'B': '#16a085', 'E': '#f1c40f', 'D': '#bfa27a', 'F': '#34495e', 'C': '#95a5a6', 'G': '#6a3d9a', 'J': '#c0392b'}
FAM_OTHER = '#d9d9d9'


def famlab(f, strand=None):
    """Label of a repeat-unit family as drawn: R with the family letter as subscript (R_A), so it is not read as an
    SMN1 lineage (A, B) or as the core; strand '+' or '-' appended (minus sign)."""
    s = '' if strand is None else ('+' if strand == '+' else '−')
    return f'R$_{{{f}}}${s}'
# Copy-tick colours (Figs. 1, 2, 4, ML tree) and nested-repeat bar colours (Figs. 1 and 4), shared for the same reason.
COPYCOL = {'A': '#e41a1c', 'B1': '#000000', 'B2': '#a65628', 'B-': '#666666', 'SMN2': '#377eb8', 'hybrid': '#984ea3', 'other': '#999999'}
NESTCOL = {'spans': '#8dc63f', 'beyond': '#e7298a', 'inside': '#999999'}   # nested inverted repeat; direct NAIP duplication
# Arm elements (step 28b), shared by Figs. 1 and 4 and the self dot plots: A = axis to the end of the core; Xi = inverted
# copy of the core's inner part; S = palindrome spacer; Xd = Xi's outer half again, opposite to Xi (direct, as the core); E = the GTF2H2-bearing stretch.
ELCOL = {'A': '#c9a978', 'Xi': '#2b6cb0', 'S': '#c51b7d', 'Xd': '#66a61e', 'E': '#e6ab02', 'Fl': '#9e9e9e', 'K': '#c9a978', 'other': '#e0e0e0'}
ELLAB = {'A': 'core, to the end of the core (also later copies of its end)', 'Xi': 'Xi′ and Xi, the inner and outer parts of an inverted copy of the core',
         'S': 'S, palindrome spacer', 'Xd': 'Xd, Xi again, opposite to Xi (direct, as the core)', 'E': 'E, carries GTF2H2',
         'Fl': 'unique flank', 'other': 'other bundles'}


# Xi' and Xi (researcher 2026-11-18): step 28b's element 'Xi' is one inverted copy of the core, split by its tandem array
# into an inner part (bundles 1-, 5-; nearer the core), written Xi' (prime = the part outside X), and an outer part
# (13- 16- 17- 21- 24- 3-), written Xi, which with S and Xd makes X. Every Xi block of arm_element_blocks.tsv holds
# bundles of one part only (checked in 40_numbers.py), so a block is labelled by the part its bundles belong to.
XI_INNER = {'1-', '5-'}
XI_PRIME = 'Xi\u2032'
# Xi' is drawn in a paler blue than Xi (researcher 2026-11-20), so X's boundary can be read off the element strip
# of Figs. 1, 4 and the self dot plots without reading the labels.
ELCOL[XI_PRIME] = '#a6cee3'
ELLAB['Xi'] = 'Xi, the outer part of an inverted copy of the core (with S and Xd it makes X)'
ELLAB[XI_PRIME] = 'Xi\u2032, the inner part of that copy, nearer the core (outside X)'
EL_DARK = {'Xi', 'S', 'Xd'}   # blocks dark enough to take white text; Xi' and the others take black
def eltxt(lab):
    """Text color for a label drawn inside an element block."""
    return 'white' if lab in EL_DARK else 'black'
def elcol(el, lab=None):
    """Block color: by the display label where that has a color of its own (Xi\u2032), else by the element."""
    return ELCOL.get(lab, ELCOL.get(el, '#d9d9d9'))
def xi_part_labels(blk, seg_path='results/sv/arm_element_segments.tsv'):
    """Return a Series (index of blk) of display labels: el, with Xi blocks of inner bundles labelled Xi'."""
    import pandas as _pd
    sg = _pd.read_csv(seg_path, sep='\t'); sg = sg[sg.el == 'Xi']
    lab = blk.el.astype(str).copy()
    for k, r in blk[blk.el == 'Xi'].iterrows():
        t = sg[(sg.src == r.src) & (sg.side == r.side) & (sg.a0 >= r.a0 - 0.5) & (sg.a1 <= r.a1 + 0.5)]
        if len(t) and t.tok.isin(XI_INNER).all(): lab[k] = XI_PRIME
    return lab


def ml_smn1_hap(hs, O='results/sv'):
    """The gene-body maximum-likelihood tree (step 04a, results/sv/mltree.treefile) at the haplotype level, for the
    haplotypes hs (each must carry at least one SMN1 tip): returns (H, order).
    H[i, j] = patristic distance between the SMN1 copy sets of hs[i] and hs[j], each copy to its nearest copy on the
    other haplotype, both ways, averaged (as step 07 combines p-distances).
    order = indices of hs in the tree's tip order (each haplotype at its first SMN1 tip), for runs of leaves."""
    import sys as _sys, numpy as _np, pandas as _pd
    _sys.setrecursionlimit(100000)
    tree = parse_newick(open(f'{O}/mltree.treefile').read().strip().rstrip(';'))
    path, tiporder = {}, []
    def walk(nd, anc, d):
        me = anc + [(id(nd), d)]
        if not nd['kids']:
            path[nd['name']] = me; tiporder.append(nd['name']); return
        for k in nd['kids']: walk(k, me, d + k['len'])
    walk(tree, [], 0.0)
    tips = _pd.read_csv(f'{O}/copy_tips.tsv', sep='\t'); s1 = tips[(tips.copy_class == 'SMN1') & tips.tip.isin(path)]
    byh = {h: list(s1[s1.src == h].tip) for h in hs}
    assert all(byh[h] for h in hs), 'a haplotype has no SMN1 tip in the ML tree'
    allt = sorted({x for h in hs for x in byh[h]}); ti = {x: i for i, x in enumerate(allt)}
    anc = [{nid for nid, _ in path[x]} for x in allt]
    P = _np.zeros((len(allt), len(allt)))
    for i, a in enumerate(allt):
        for j in range(i + 1, len(allt)):
            lca = 0.0
            for nid, d in path[a]:
                if nid in anc[j]: lca = d
                else: break
            P[i, j] = P[j, i] = path[a][-1][1] + path[allt[j]][-1][1] - 2 * lca
    m = len(hs); H = _np.zeros((m, m))
    for i in range(m):
        for j in range(i + 1, m):
            sub = P[_np.ix_([ti[x] for x in byh[hs[i]]], [ti[x] for x in byh[hs[j]]])]
            H[i, j] = H[j, i] = (sub.min(1).mean() + sub.min(0).mean()) / 2
    pos = {x: k for k, x in enumerate(tiporder)}
    first = [min(pos[x] for x in byh[h]) for h in hs]
    order = list(_np.argsort(first, kind='stable'))
    return H, order
