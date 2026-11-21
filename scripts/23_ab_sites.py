# Pipeline step 23 (ab_sites); ported from explore/structure_vs_sequence/20_ab_sites_along_gene.py (sha256 5f4d30ca02f4, commit 5e2b3bb),
# paths rewritten to results/sv/ and figures/; run from the project root after the steps before it.
"""Exploration step 20: every alignment column of the SMN1 gene body where the majority base (a gap counts as a state)
differs between lineage A and lineage B (ML halves, step 05), drawn along the gene with the exons of NM_000344.4.
Reads out/copies.aln.fa, out/copy_tips.tsv, doc/marker_lineage.tsv, out/outgroup_polarized_0.9.tsv,
data/refs_seq/smn1_smn2_refseq.gb. Writes out/ab_sites.tsv, out/ab_site_events.tsv, fig/ab_sites.png/.pdf.
Note: the lineages were cut from a tree built on this same alignment, so how cleanly a site separates them is partly built in."""
import sys, collections, numpy as np, pandas as pd, matplotlib
sys.path.insert(0, 'scripts'); from smnlib import read_genbank
matplotlib.use('Agg'); import matplotlib.pyplot as plt
O = 'results/sv'
names, seqs = [], []
for l in open(f'{O}/copies.aln.fa'):
    l = l.strip()
    if l.startswith('>'): names.append(l[1:].split()[0]); seqs.append([])
    else: seqs[-1].append(l.upper())
M = np.array([np.frombuffer(''.join(s).encode(), dtype=np.uint8) for s in seqs]); idx = {n: i for i, n in enumerate(names)}
t = pd.read_csv(f'{O}/copy_tips.tsv', sep='\t'); ml = pd.read_csv('results/sv/marker_lineage.tsv', sep='\t').set_index('tip')
ref = idx[t[(t.copy_class == 'SMN1') & t.src.str.startswith('GCA_000001405')].tip.iloc[0]]
rseq = bytes(M[ref][M[ref] != ord('-')]).decode()
rpos = np.cumsum(M[ref] != ord('-')) - 1; L = rpos[-1] + 1
OFF = 70923086                     # GRCh38 position of the reference copy's first base minus 1 (as in step 09)
gene = (rpos >= 2000) & (rpos < L - 2000)
G = {'A': [], 'B': [], 'SMN2': []}
for r in t.itertuples():
    if r.tip not in idx: continue
    if r.copy_class == 'SMN1' and r.tip in ml.index: G['B' if ml.loc[r.tip, 'in_marker_subclade'] else 'A'].append(idx[r.tip])
    elif r.copy_class == 'SMN2': G['SMN2'].append(idx[r.tip])
nA, nB = len(G['A']), len(G['B'])
ml_ = pd.read_csv(f'{O}/marker_lineage.tsv', sep='\t'); ml_ = ml_[ml_['class'] == 'SMN1']
assert (nA, nB) == (int((~ml_.in_marker_subclade).sum()), int(ml_.in_marker_subclade.sum())), {k: len(v) for k, v in G.items()}
# exons of NM_000344.4 located in the reference copy
g = read_genbank('data/refs_seq/smn1_smn2_refseq.gb')['NM_000344.4']
EX = []
for name, (a, b) in zip(['1', '2a', '2b', '3', '4', '5', '6', '7', '8'], g['exons']):
    p = rseq.find(g['seq'][a - 1:b]); assert p >= 0, name
    EX.append((name, OFF + p + 1, OFF + p + (b - a + 1)))
# c.840 must fall in exon 7 at chr5:70,951,946 (the note's position) - checks the coordinate frame
ex7 = [e for e in EX if e[0] == '7'][0]; cds_a = g['cds'][0]
c840 = ex7[1] + (cds_a + 840 - 1) - g['exons'][7][0]
assert c840 == 70951946, c840
rows = []
for c in np.where(gene)[0]:
    ca = collections.Counter(M[G['A'], c].tobytes().decode()); cb = collections.Counter(M[G['B'], c].tobytes().decode())
    a = ca.most_common(1)[0][0]; b = cb.most_common(1)[0][0]
    if a == b: continue
    s2 = collections.Counter(M[G['SMN2'], c].tobytes().decode())
    rows.append(dict(col=int(c), grch38=int(OFF + rpos[c] + 1), ref_gap=bool(M[ref, c] == ord('-')), A=a, fA=ca[a] / nA, B=b, fB=cb[b] / nB,
                     B_with_A=cb[a], A_with_B=ca[b], SMN2_as_A=s2[a] / len(G['SMN2']), SMN2_as_B=s2[b] / len(G['SMN2'])))
S = pd.DataFrame(rows)
S['kind'] = np.where((S.A == '-') | (S.B == '-'), 'indel', 'SNV')
S['sep'] = np.minimum(S.fA, S.fB)            # how cleanly the site splits the lineages (both sides)
S['class'] = np.where(S.sep >= 0.9, 'near-fixed (both >=90%)', 'partial')
S['SMN2_shares'] = np.where(S.SMN2_as_A >= 0.5, 'A', np.where(S.SMN2_as_B >= 0.5, 'B', 'neither'))
def region(p):
    for i, (n, a, b) in enumerate(EX):
        if a <= p <= b: return 'exon ' + n
    for i in range(len(EX) - 1):
        if EX[i][2] < p < EX[i + 1][1]: return 'intron ' + EX[i][0]
    return 'outside'
S['region'] = S.grch38.map(region)
# events: adjacent columns of one indel (or one gap run) merged
S = S.sort_values('col').reset_index(drop=True)
ev = (S.col.diff() != 1) | (S.kind != S.kind.shift()); S['event'] = ev.cumsum()
E = S.groupby('event').agg(grch38=('grch38', 'min'), cols=('col', 'size'), kind=('kind', 'first'), sep=('sep', 'max'), cls=('class', lambda x: x.iloc[0] if x.nunique() == 1 else 'near-fixed (both >=90%)' if (x == 'near-fixed (both >=90%)').any() else 'partial'),
                           SMN2_shares=('SMN2_shares', 'first'), region=('region', 'first')).reset_index()
P9 = pd.read_csv(f'{O}/outgroup_polarized_0.9.tsv', sep='\t')
nf = S[(S['class'].str.startswith('near')) & (S.kind == 'SNV')]
assert set(P9.col) <= set(nf.col), 'step 09 sites missing'
S.to_csv(f'{O}/ab_sites.tsv', sep='\t', index=False); E.to_csv(f'{O}/ab_site_events.tsv', sep='\t', index=False)
perfect = S[(S.fA == 1) & (S.fB == 1)]   # no exception: every A copy carries A's state and every B copy B's
print('columns', len(S), 'events', len(E), E.kind.value_counts().to_dict(), 'near-fixed events', int(E.cls.str.startswith('near').sum()),
      'SNV near-fixed cols', len(nf), 'perfect cols', perfect[['col', 'grch38', 'A', 'B', 'region']].to_dict('records'))
print(E.groupby(['cls', 'SMN2_shares']).size()); print(E.region.value_counts().to_dict())
# marker c.*3+80T>G
mk = S[S.grch38 == 70952074]; print('marker row', mk[['A', 'fA', 'B', 'fB', 'sep']].to_dict('records'))
# --- figure
fig, (ax0, ax) = plt.subplots(2, 1, figsize=(12, 5.2), sharex=True, gridspec_kw=dict(height_ratios=[1, 4], hspace=0.05))
x0, x1 = OFF + 2000, OFF + L - 2000
ax0.plot([x0, x1], [0, 0], color='k', lw=1)
for n, a, b in EX:
    ax0.add_patch(plt.Rectangle((a, -0.35), b - a, 0.7, color='#333')); ax0.text((a + b) / 2, 0.55, n, ha='center', fontsize=8)
ax0.axvspan(70950189, 70952674, color='#f5b041', alpha=0.25, lw=0); ax0.text(70950100, -0.9, 'shaded: SMN1/SMN2 fixed differences (intron 6 - exon 8)', ha='right', va='center', fontsize=7, color='#9c640c')
ax0.set_ylim(-1.6, 1.1); ax0.set_yticks([]); ax0.set_title('SMN1 gene body (GRCh38 copy, NM_000344.4 exons) and every site where lineage A and B differ in their majority base', loc='left', fontsize=10)
for s in ['top', 'right', 'left']: ax0.spines[s].set_visible(False)
col = {'A': '#2471a3', 'B': '#c0392b', 'neither': '#7f8c8d'}
for (k, sh), q in E.groupby(['kind', 'SMN2_shares']):
    ax.vlines(q.grch38, 0, q.sep, color=col[sh], lw=0.7, alpha=0.7)
    ax.scatter(q.grch38, q.sep, s=22, marker='o' if k == 'SNV' else 's', color=col[sh], edgecolor='k', lw=0.3, zorder=3,
               label=f'{k}, SMN2 carries {"A" if sh == "A" else "B" if sh == "B" else "neither"}\'s state ({len(q)})')
ax.axhline(0.9, color='grey', ls=':', lw=0.8); ax.text(x1, 0.905, '90%', fontsize=7, color='grey', ha='right', va='bottom')
ax.axvspan(70950189, 70952674, color='#f5b041', alpha=0.15, lw=0)
m = E[E.grch38 == 70952074]
if len(m): ax.annotate('c.*3+80T>G', (70952074, m.sep.iloc[0]), xytext=(70947500, 0.35), fontsize=8, arrowprops=dict(arrowstyle='->', lw=0.6))
d = E[(E.grch38 >= 70952640) & (E.grch38 <= 70952650) & (E.kind == 'indel')]
assert len(d) == 1
ax.annotate('c.*211_*212del', (d.grch38.iloc[0], d.sep.iloc[0]), xytext=(70948800, 0.18), fontsize=8, arrowprops=dict(arrowstyle='->', lw=0.6))
for r in E[E.sep >= 0.995].itertuples():
    ax.annotate(f'{r.grch38:,}', (r.grch38, r.sep), xytext=(0, 6), textcoords='offset points', fontsize=7, ha='center', rotation=90, va='bottom')
ax.set_ylim(0, 1.18); ax.set_xlim(x0 - 200, x1 + 200)
ax.set_ylabel('separation: share of the less-fixed\nlineage carrying its own state')
ax.set_xlabel('GRCh38 chr5 position (SMN1 copy)'); ax.ticklabel_format(axis='x', style='plain', useOffset=False)
ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f'{v/1e6:.3f} Mb'))
ax.legend(fontsize=7, loc='upper left', bbox_to_anchor=(0, -0.2), ncol=3, frameon=False)
fig.savefig('figures/ab_sites.png', dpi=200, bbox_inches='tight'); fig.savefig('figures/ab_sites.pdf', bbox_inches='tight', metadata={'CreationDate': None})
