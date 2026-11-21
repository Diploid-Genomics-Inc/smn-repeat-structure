"""Pipeline step 34 (supplement): (1) self-against-self shimmer dot plots of one haplotype per palindrome form (step 28b:
cenSMN arm C1/C0 x telSMN arm T0/T1; C1|T0 and C0|T1 are Fig. 1's haplotypes, results/sv/locus_diagram_examples.tsv),
with the arm elements (results/sv/arm_element_blocks.tsv), S at the axis, the ~20 kb arrays and the SMN copies along
the axes; (2) for illustration only, PANGEA's principal MAP graph (`pangea cohort bundle decomp`, the PGR-TK graph)
built from five fixed windows: HG03816 pat, HG01099 mat, GRCh38, CHM13 and the alphabetically first '1- 1+' haplotype
of the B1+B2 band (results/sv/ml_hap_labels.tsv label S1-8+S1-9). Window sequences come from the HPRC release-2
archive ($STRATA_DATA_HPRC/hprc_r2_full.agcrs) at results/copies/windows.tsv, oriented as the study's windows, in
data/derived/windows/ (not in git). Writes figures/self_dotplots.png/.pdf, figures/window_bundles.png/.pdf,
results/sv/pangenome_graph.gexf (for Gephi) and results/sv/dotplot_graph_summary.json."""
import os, re, subprocess, json, sys
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D
import pangeatk as ptk, networkx as nx
sys.path.insert(0, 'scripts'); from smnlib import FAMCOL, famlab
O, WD = 'results/sv', 'data/derived/windows'; os.makedirs(WD, exist_ok=True)
ARC = os.path.join(os.environ['STRATA_DATA_HPRC'], 'hprc_r2_full.agcrs')
W = pd.read_csv('results/copies/windows.tsv', sep='\t'); W = W[W.n_copies > 0].set_index('src')
EX = pd.read_csv(f'{O}/locus_diagram_examples.tsv', sep='\t').set_index('panel')
RU = pd.read_csv('results/decomp/smn.repeat_units.tsv', sep='\t', comment='#')
CU = pd.read_csv('results/decomp/copy_units.tsv', sep='\t')
HL = pd.read_csv(f'{O}/ml_hap_labels.tsv', sep='\t')
COMP = str.maketrans('ACGTN', 'TGCAN')
def window(src):
    f = f'{WD}/{src}.fa'
    if not os.path.exists(f):
        r = W.loc[src]
        out = subprocess.run(['agc-rs', 'get', ARC, f'{src}/{r.ctg}:{int(r.ctg_start)}-{int(r.ctg_end)}'],
                             capture_output=True, text=True, check=True).stdout
        s = ''.join(l.strip() for l in out.splitlines() if not l.startswith('>')).upper()
        assert len(s) == int(r.length), (src, len(s), r.length)
        if int(r.orientation) == 1: s = s.translate(COMP)[::-1]
        open(f, 'w').write(f'>{src}\n{s}\n')
    return ''.join(l.strip() for l in open(f) if not l.startswith('>'))
summ = {'shimmer': dict(w=24, k=48, r=2, min_span=12)}
# ---- (1) self dot plots, one haplotype per palindrome form (step 28b). C1|T0 and C0|T1 are Fig. 1's two haplotypes;
# for C1|T1 and C0|T0, the haplotype with both arms in their type's commonest layout whose window length is nearest
# the form's median (ties by name).
EH = pd.read_csv(f'{O}/arm_element_haps.tsv', sep='\t').set_index('src'); EA = pd.read_csv(f'{O}/arm_elements.tsv', sep='\t')
BLK = pd.read_csv(f'{O}/arm_element_blocks.tsv', sep='\t'); SPC = pd.read_csv(f'{O}/arm_element_spacers.tsv', sep='\t')
ARR = pd.read_csv(f'{O}/arm_element_arrays.tsv', sep='\t'); AXR = pd.read_csv(f'{O}/palindrome_diagram_rows.tsv', sep='\t').set_index('src')
CPL = pd.read_csv(f'{O}/palindrome_copies_lineage.tsv', sep='\t', keep_default_na=False)
from smnlib import ELCOL, ELLAB, COPYCOL, elcol, eltxt, xi_part_labels, XI_PRIME
BLK['lab'] = xi_part_labels(BLK)   # Xi' (inner part) or Xi (outer part), display only
canon = EA.groupby('src').canonical.all(); EH['form'] = EH.cen_type + '|' + EH.tel_type
reps = {'C1|T0': EX.loc['a', 'src'], 'C0|T1': EX.loc['b', 'src']}
for f in ('C1|T1', 'C0|T0'):
    m = EH[(EH.form == f) & EH.index.map(canon)].copy(); m['L'] = [W.loc[s, 'length'] for s in m.index]
    reps[f] = (m.L - m.L.median()).abs().reset_index().sort_values(['L', 'src'], kind='stable').src.iloc[0]
assert EH.loc[reps['C1|T0'], 'form'] == 'C1|T0' and EH.loc[reps['C0|T1'], 'form'] == 'C0|T1'
GCOL = {'SMN1 A': COPYCOL['A'], 'SMN1 B1': COPYCOL['B1'], 'SMN1 B2': COPYCOL['B2'], 'SMN2': COPYCOL['SMN2']}
fig, axs = plt.subplots(2, 2, figsize=(15, 16.2)); axs = axs.ravel()
for i, (ax, f) in enumerate(zip(axs, ['C1|T0', 'C1|T1', 'C0|T0', 'C0|T1'])):
    src = reps[f]; s = window(src); L = len(s) / 1e3; ax_ = AXR.axis[src] / 1e3
    x, y = ptk.get_shmmr_dots(s.encode(), s.encode(), 24, 48, 2, 12)
    x, y = np.asarray(x) / 1e3, np.asarray(y) / 1e3
    off = np.abs(x - y) > 5
    ax.scatter(x, y, s=0.05, c='k', lw=0, rasterized=True)
    def strip(lo, hi, col, o=.045, w=.03):
        ax.add_patch(Rectangle((lo, -L * o), hi - lo, L * w, color=col, lw=0, clip_on=False))
        ax.add_patch(Rectangle((-L * o, lo), L * w, hi - lo, color=col, lw=0, clip_on=False))
    for b in BLK[BLK.src == src].itertuples():
        lo, hi = (ax_ - b.a1, ax_ - b.a0) if b.side == 'cen' else (ax_ + b.a0, ax_ + b.a1)
        strip(lo, hi, elcol(b.el, b.lab))
        if b.el in ('Xi', 'S', 'Xd'): ax.add_patch(Rectangle((lo, 0), hi - lo, L, color=elcol(b.el, b.lab), alpha=.06, lw=0))
        if b.el in ('Xi', 'Xd', 'E') and hi - lo > 40: ax.text((lo + hi) / 2, -L * .03, b.lab, ha='center', va='center', fontsize=8, color=eltxt(b.lab), clip_on=False)
    for q in SPC[(SPC.src == src) & (SPC.site == 'axis')].itertuples(): strip(q.bgn / 1e3, q.end / 1e3, ELCOL['S'])
    for r in ARR[ARR.src == src].itertuples():
        lo, hi = (ax_ - r.d1, ax_ - r.d0) if r.side == 'cen' else (ax_ + r.d0, ax_ + r.d1)
        strip(lo, hi, '#333', o=.01, w=.008)
    for c in CPL[CPL.src == src].itertuples():
        p_ = ax_ + c.d_kb; ax.plot(p_, p_, 'o', ms=8, mfc='none', mec=GCOL.get(c.mark, COPYCOL['other']), mew=2)
    ax.axvline(ax_, color='k', lw=.6, ls='--'); ax.axhline(ax_, color='k', lw=.6, ls='--')
    ax.tick_params(axis='y', pad=24); ax.set_xlim(0, L); ax.set_ylim(L, 0); ax.set_aspect('equal'); ax.set_xlabel(f'{src.split("_hprc")[0]} window (kb)'); ax.set_ylabel('same window (kb)')
    nf = int((EH.form == f).sum())
    ax.set_title(f'({"abcd"[i]}) cenSMN arm {f[:2]}, telSMN arm {f[3:]}: {src.split("_hprc")[0].replace("_", " ")} ({EH.loc[src, "superpop"]})\n'
                 f'{nf} haplotypes have this form', loc='left', fontsize=10.5, pad=34)
    ae = {sd: EA[(EA.src == src) & (EA.side == sd)].a_end_kb.iloc[0] for sd in ('cen', 'tel')}
    beyond_cen = off & (x < ax_ - ae['cen']) & (y < ax_ - ae['cen'])   # both ends beyond A on the cenSMN arm
    beyond_tel = off & (x > ax_ + ae['tel']) & (y > ax_ + ae['tel'])
    summ['abcd'[i]] = dict(form=f, src=src, length=len(s), dots=int(len(x)), off_diagonal=int(off.sum()), axis_kb=float(ax_),
                           beyond_A_cen_dots=int(beyond_cen.sum()), beyond_A_tel_dots=int(beyond_tel.sum()))
hd = [Patch(color=ELCOL[e], label=ELLAB[e]) for e in ['A', XI_PRIME, 'Xi', 'S', 'Xd', 'E', 'Fl', 'other']]
hd += [Line2D([], [], color='#333', lw=3, label='~20 kb tandem array (inner strip)'), Line2D([], [], color='k', ls='--', lw=.8, label='R$_A$+ R$_A$− axis')]
hd += [Line2D([], [], marker='o', ls='', mfc='none', mec=c, mew=2, label=k) for k, c in GCOL.items()]
fig.legend(handles=hd, loc='lower center', ncol=4, fontsize=9.5, frameon=False, bbox_to_anchor=(.5, .005))
fig.tight_layout(rect=(0, .065, 1, .985), h_pad=4)
fig.savefig('figures/self_dotplots.png', dpi=150); fig.savefig('figures/self_dotplots.pdf', dpi=150, metadata={'CreationDate': None}); plt.close(fig)
# ---- (2) principal MAP graph of five windows (illustration)
band = sorted(HL[(HL.seq_label == 'S1-8+S1-9') & (HL.arrangement == '1- 1+')].src)[0]
refs = [s for s in W.index if s.startswith('GCA_')] + [s for s in W.index if s.startswith('chm13')]
srcs = ['HG03816_pat_hprc_r2_v1.0.1', 'HG01099_mat_hprc_r2_v1.0.1'] + refs + [band]   # fixed: the five windows of Fig. S1 (chosen 2026-11-14, before Fig. 1 changed)
os.makedirs('data/derived/graph', exist_ok=True); fa = 'data/derived/graph/five_windows.fa'
with open(fa, 'w') as f:
    for s in srcs: f.write(f'>{s}\n{window(s)}\n')
pre = 'data/derived/graph/five'
subprocess.run(['pangea', 'cohort', 'bundle', 'decomp', '-f', fa, '-o', pre, '-w', '48', '-r', '4',
                '--bundle-merge-distance', '10000', '--bundle-length-cutoff', '2500'], check=True, capture_output=True)
G = nx.Graph(); bn = {}; seglen = {}
for l in open(f'{pre}.pmapg.gfa'):
    t = l.rstrip('\n').split('\t')
    if t[0] == 'S':
        tags = {x.split(':', 2)[0]: x.split(':', 2)[2] for x in t[3:]}; G.add_node(t[1]); bn[t[1]] = int(tags.get('BN', -1)); seglen[t[1]] = int(tags.get('LN', 0))
    elif t[0] == 'L':
        G.add_edge(t[1], t[3])
BED = pd.read_csv(f'{pre}.bed', sep='\t', comment='#', header=None, names=['ctg', 'bgn', 'end', 'tag'])
BED['bid'] = BED.tag.str.split(':').str[0].astype(int); BED['dir'] = BED.tag.str.split(':').str[2].astype(int)
top = list(BED.assign(bp=BED.end - BED.bgn).groupby('bid').bp.sum().sort_values(ascending=False, kind='stable').index[:20])
cmap = plt.get_cmap('tab20'); bcol = lambda b: cmap(top.index(b)) if b in top else '#cccccc'   # graph export (GEXF) colours only
# Fig. S1 draws the STUDY's decomposition (researcher 2026-11-15: the figure must use the coarse bundles of the paper),
# so its bundle ids are those of Table S1; the five-window decomposition above is kept only for the GFA/GEXF export.
SBED = pd.read_csv('results/decomp/smn.bed', sep='\t', comment='#', header=None, names=['ctg', 'bgn', 'end', 'tag'])
SBED = SBED[SBED.ctg.isin(srcs)].copy()
SBED['bid'] = SBED.tag.str.split(':').str[0].astype(int); SBED['dir'] = SBED.tag.str.split(':').str[2].astype(int)
assert set(SBED.ctg) == set(srcs), 'every Fig. S1 window must be in the study decomposition'
sids = list(SBED.assign(bp=SBED.end - SBED.bgn).groupby('bid').bp.sum().sort_values(ascending=False, kind='stable').index)
assert len(sids) <= 40
scol = lambda b: (plt.get_cmap('tab20')(sids.index(b)) if sids.index(b) < 20 else plt.get_cmap('tab20b')(sids.index(b) - 20))
ids = sorted({b for b in bn.values() if b >= 0})
# researcher 2026-11-14: 'remove the graph, not very useful' -- the graph is no longer drawn; it is still exported
# (GFA, GEXF below). The figure keeps the windows as bundles, repeat units and arm elements.
fig = plt.figure(figsize=(14, 5.2)); bx = fig.add_axes([.11, .11, .86, .8])
nodes = list(G.nodes)
YS = 1.35   # row spacing in panel b (bundles, units, elements per window)
for i, s in enumerate(srcs):
    sub = SBED[SBED.ctg == s]
    for _, r in sub.iterrows():
        bx.add_patch(Rectangle((r.bgn / 1e3, -YS * i - .05), (r.end - r.bgn) / 1e3, .4, color=scol(r.bid), lw=0, alpha=1 if r.dir == 0 else .55))
    if s in AXR.index:   # arm elements (step 28b) under the units: the same elements whatever the unit call
        ax0 = AXR.axis[s]
        for b in BLK[BLK.src == s].itertuples():
            x0, x1 = (ax0 - b.a1 * 1e3, ax0 - b.a0 * 1e3) if b.side == 'cen' else (ax0 + b.a0 * 1e3, ax0 + b.a1 * 1e3)
            bx.add_patch(Rectangle((x0 / 1e3, -YS * i - .72), (x1 - x0) / 1e3, .24, color=elcol(b.el, b.lab), lw=0))
            if x1 - x0 > 60000: bx.text((x0 + x1) / 2e3, -YS * i - .60, 'core' if b.el == 'A' else b.lab, ha='center', va='center', fontsize=6, color='w' if (b.lab in ('Xi', 'Xd') or b.el == 'A') else 'k')
        for q in ARR[ARR.src == s].itertuples():
            x0, x1 = (ax0 - q.d1 * 1e3, ax0 - q.d0 * 1e3) if q.side == 'cen' else (ax0 + q.d0 * 1e3, ax0 + q.d1 * 1e3)
            bx.plot([x0 / 1e3, x1 / 1e3], [-YS * i - .76] * 2, color='k', lw=2, solid_capstyle='butt')
            summ.setdefault('graph_arrays', []).append(dict(src=s, side=q.side, where=q.where, unit_kb=q.unit_kb, copies=q.copies, kb=round(q.d1 - q.d0, 1)))
    for _, u in RU[RU.row == s].iterrows():   # the study's repeat units (from its own decomposition, smn.bed) on the same window
        assert u.end <= len(window(s)) + 1
        bx.add_patch(Rectangle((u.bgn / 1e3, -YS * i - .42), (u.end - u.bgn) / 1e3, .3, color=FAMCOL.get(u.family, '#d9d9d9'), alpha=1 if u.strand == '+' else .45, lw=0))
        if u.end - u.bgn > 80000: bx.text((u.bgn + u.end) / 2e3, -YS * i - .27, famlab(u.family, u.strand), ha='center', va='center', fontsize=7)
bx.set_yticks(-YS * np.arange(len(srcs)))
lab = {s: ('GRCh38' if s.startswith('GCA_') else 'CHM13' if s.startswith('chm13') else s.split('_hprc')[0]) for s in srcs}
bx.set_yticklabels([lab[s] for s in srcs]); bx.set_ylim(-YS * (len(srcs) - 1) - .9, .5)
bx.set_xlim(0, max(len(window(s)) for s in srcs) / 1e3); bx.set_xlabel('window (kb)')
bx.set_title('Five windows as bundle segments (top bar), repeat units (middle) and arm elements (bottom)', loc='left')
for sp in ['top', 'right']: bx.spines[sp].set_visible(False)
fig.savefig('figures/window_bundles.png', dpi=200); fig.savefig('figures/window_bundles.pdf', dpi=200, metadata={'CreationDate': None}); plt.close(fig)
for n in nodes:   # GEXF for Gephi: bundle id, colour as drawn, vertex length
    c = matplotlib.colors.to_rgb(bcol(bn[n])); G.nodes[n]['bundle'] = bn[n]; G.nodes[n]['top_bundle'] = bn[n] in top
    G.nodes[n]['viz'] = {'color': {'r': int(c[0] * 255), 'g': int(c[1] * 255), 'b': int(c[2] * 255), 'a': 1.0}, 'size': 1.0}
for n, l in seglen.items(): G.nodes[n]['length'] = l
# haplotype paths (2026-11-14, researcher: GFA and GEXF of the MAP graph): each window's shimmer pairs, with the
# decomp's settings, matched to the vertex names (SN tag = the two shimmer hashes, smaller first)
sn2v = {}
for l in open(f'{pre}.pmapg.gfa'):
    t = l.rstrip('\n').split('\t')
    if t[0] == 'S': sn2v[[x for x in t[3:] if x.startswith('SN:Z:')][0][5:]] = t[1]
LS = set()
for l in open(f'{pre}.pmapg.gfa'):
    t = l.split('\t')
    if t[0] == 'L': LS.add((t[1], t[2], t[3], t[4])); LS.add((t[3], '-' if t[4] == '+' else '+', t[1], '-' if t[2] == '+' else '+'))
paths, pstat = {}, {}
for s_ in srcs:
    steps = []
    for h0, h1, p0, p1, o in ptk.get_shmmr_pairs_from_seq(window(s_).encode(), 48, 56, 4, 12):
        v = sn2v.get(f'{min(h0, h1):016x}_{max(h0, h1):016x}')
        if v is not None and (not steps or steps[-1][0] != v): steps.append((v, '+' if o == 0 else '-'))
    ok = sum((a_[0], a_[1], b_[0], b_[1]) in LS for a_, b_ in zip(steps, steps[1:]))
    paths[s_] = steps; pstat[lab[s_]] = dict(steps=len(steps), steps_on_edges=ok)
    for v, _ in steps: G.nodes[v]['hap_' + lab[s_].replace(' ', '_')] = True
for n in G.nodes:
    for s_ in srcs: G.nodes[n].setdefault('hap_' + lab[s_].replace(' ', '_'), False)
    G.nodes[n]['n_haplotypes'] = sum(G.nodes[n]['hap_' + lab[s_].replace(' ', '_')] for s_ in srcs)
with open(f'{O}/pangenome_graph.gfa', 'w') as f:
    for l in open(f'{pre}.pmapg.gfa'):
        f.write(l)
        if l.startswith('H'): f.write('# principal MAP graph of Fig. S1a (pangea cohort bundle decomp w48 k56 r4); P lines added by step 34\n')
    for s_ in srcs: f.write(f"P\t{lab[s_].replace(' ', '_')}\t{','.join(v + o for v, o in paths[s_])}\t*\n")
summ['graph_paths'] = pstat
nx.write_gexf(G, f'{O}/pangenome_graph.gexf')
# networkx stamps the run's date in <meta lastmodifieddate>; fix it so the export is byte-identical on every run
_gp = f'{O}/pangenome_graph.gexf'
_gt = open(_gp).read()
open(_gp, 'w').write(re.sub(r'lastmodifieddate="[0-9-]+"', 'lastmodifieddate="2026-11-15"', _gt, count=1))
summ['graph'] = dict(haplotypes=srcs, nodes=G.number_of_nodes(), edges=G.number_of_edges(),
                     components=nx.number_connected_components(G), bundles=len(ids), bundles_coloured=len(top), bed_segments=int(len(BED)), study_segments=int(len(SBED)), study_bundles=len(sids),
                     nodes_in_bundles=int(sum(b >= 0 for b in bn.values())),
                     decomp=dict(w=48, k=56, r=4, bundle_merge_distance=10000, bundle_length_cutoff=2500))
json.dump(summ, open(f'{O}/dotplot_graph_summary.json', 'w'), indent=1)
print(json.dumps(summ, indent=1))
