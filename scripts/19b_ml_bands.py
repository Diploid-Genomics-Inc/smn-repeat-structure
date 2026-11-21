"""Pipeline step 19b: the bands that link large-scale structure to SMN1 sequence, with the sequence side taken from
the gene-body maximum-likelihood tree (step 04a: IQ-TREE GTR+F+G4, 1,000 ultrafast bootstraps) instead of a cut UPGMA
tree. Promoted from explore/ml_bands/ml_bands.py (commit 1e48190) by the researcher's decision of 2026-11-14 ("ok, use
the ML ... it may be ok to leave the UPGAMA on the supplemenatry as another check"); step 19 (UPGMA on both sides)
stays as that check.

Structure side: UPGMA (average linkage, optimal leaf order) on PANGEA bundle distance over the haplotypes with an SMN1
copy, cut into K clusters, exactly as step 19 (checked against step 19's member table).
Sequence side: each SMN1 copy goes to the smallest clade of the ML tree with UFBoot >= T that holds at least M SMN1
copies of these haplotypes; copies in no such clade go to the rest of the human clade ('rest'). A cluster is named by
the commonest Chen et al. haplogroup family of its own copies. A haplotype's sequence label is its SMN1 copy's
cluster, or for several copies the sorted '+'-joined clusters.
Test: every structure-cluster x sequence-label pair of >= 5 haplotypes, hypergeometric, Bonferroni over all
structure-cluster x sequence-label pairs. A band is a pair of >= MIN_OV haplotypes with Bonferroni p < ALPHA.
Sensitivity: the same at T in {80, 90, 95} and M in {5, 10, 20}; and each band matched to step 19's UPGMA band.
Writes results/sv/ml_seq_clusters.tsv, ml_hap_labels.tsv, ml_cluster_pairs.tsv, ml_cluster_pairs_members.tsv,
ml_band_sensitivity.tsv, ml_vs_upgma_bands.tsv, ml_bands_summary.json. Fig. 5 (figures/tanglegram_ml.png/.pdf) is drawn by step
28j, which reruns this script after step 28b so the figure can carry the arm types C0/C1, T0/T1."""
import sys, json, collections, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, 'scripts'); from smnlib import parse_newick, COPYCOL
from scipy.cluster.hierarchy import average, optimal_leaf_ordering, dendrogram, fcluster
from scipy.spatial.distance import squareform
from scipy.stats import hypergeom
sys.setrecursionlimit(100000)
O, F = 'results/sv', 'figures'
K, T, M, MIN_OV, ALPHA = 16, 95, 20, 10, 1e-3
DRAW_FIG5 = globals().get('DRAW_FIG5', False)   # True only when run by step 28j (needs step 28b's arm types)
import tempfile
OUT = tempfile.mkdtemp() if DRAW_FIG5 else O   # step 28j's rerun leaves step 19b's tables to step 19b (one writer per file)
# --- structure side, as step 19
P = pd.read_csv(f'{O}/pair_distances.tsv.gz', sep='\t')
rows = sorted(set(P.a) | set(P.b)); ri = {r: i for i, r in enumerate(rows)}; n0 = len(rows)
def mat(col):
    X = np.full((n0, n0), np.nan); np.fill_diagonal(X, 0)
    ia, ib = P.a.map(ri).to_numpy(), P.b.map(ri).to_numpy(); X[ia, ib] = X[ib, ia] = P[col].to_numpy(); return X
S, G1 = mat('struct'), mat('seq_SMN1')
k1 = np.where((~np.isnan(G1)).sum(1) > 1)[0]; rk = [rows[i] for i in k1]; n = len(rk)
S1 = S[np.ix_(k1, k1)]
Zs = optimal_leaf_ordering(average(squareform(S1)), squareform(S1))
cs = fcluster(Zs, K, 'maxclust')
leafS = dendrogram(Zs, no_plot=True)['leaves']; posS = np.empty(n); posS[leafS] = np.arange(n) / (n - 1)
cm19 = pd.read_csv(f'{O}/cluster_pairs_members.tsv', sep='\t'); cm19 = cm19[cm19.k == K]
c19 = cm19.drop_duplicates('src').set_index('src').struct_cluster
assert all(c19[r] == cs[i] for i, r in enumerate(rk) if r in c19.index), 'structure clusters differ from step 19'
# --- annotations
L = pd.read_csv(f'{O}/lineage_haplotypes.tsv', sep='\t').set_index('src')
arr = np.array([L.loc[r, 'arrangement'] if r in L.index else 'NA' for r in rk], dtype=object)
sp = np.array([L.loc[r, 'superpop'] if r in L.index else 'reference' for r in rk], dtype=object)
tips = pd.read_csv(f'{O}/copy_tips.tsv', sep='\t')
ml = pd.read_csv(f'{O}/marker_lineage.tsv', sep='\t').set_index('tip')
hg = pd.read_csv(f'{O}/chen_hg_copies.tsv', sep='\t').set_index('tip')
sm1 = tips[(tips.copy_class == 'SMN1') & tips.src.isin(rk)]
sm1tips = set(sm1.tip); tip2src = dict(zip(sm1.tip, sm1.src))
sub = ml.sublineage.fillna('').to_dict(); fam = hg.family.to_dict()
# --- ML tree
tree = parse_newick(open(f'{O}/mltree.treefile').read().strip().rstrip(';'))
def walk(nd, out):
    nd['tips'] = {nd['name']} if not nd['kids'] else set().union(*[walk(k, out) for k in nd['kids']])
    nd['n1'] = len(nd['tips'] & sm1tips)
    if nd['kids']: out.append(nd)
    return nd['tips']
nodes = []; walk(tree, nodes)
human = set(ml.index)
HC = min((nd for nd in nodes if human <= nd['tips']), key=lambda nd: len(nd['tips']))
def top(v, m=3):
    s = pd.Series(list(v)).value_counts(); return '; '.join(f'{a} {b}' for a, b in s.head(m).items())
def partition(t, m):
    Q = [nd for nd in nodes if nd is HC or (nd['sup'] is not None and nd['sup'] >= t and nd['n1'] >= m and nd['tips'] <= HC['tips'])]
    cid = {}
    for nd in sorted(Q, key=lambda nd: -len(nd['tips'])):   # larger first, smaller overwrite
        for x in nd['tips'] & sm1tips: cid[x] = id(nd)
    names, info, used = {}, [], collections.Counter()
    for nd in sorted(Q, key=lambda nd: (-len(nd['tips']), min(nd['tips']))):
        own = sorted(x for x in nd['tips'] & sm1tips if cid[x] == id(nd))
        if not own: continue
        if nd is HC: nm = 'rest'
        else:
            f = pd.Series([fam.get(x, '?') for x in own]).value_counts().index[0]
            used[f] += 1; nm = f if used[f] == 1 else f'{f}/{used[f]}'
        names[id(nd)] = nm
        info.append({'T': t, 'M': m, 'cluster': nm, 'support': nd['sup'] if nd is not HC else None, 'clade_SMN1': nd['n1'],
                     'own_SMN1': len(own), 'sublineages': top([sub.get(x, '') for x in own], 5),
                     'haplogroups': top([fam.get(x, '?') for x in own], 4), 'tips': own})
    hl = {}
    for x, s in tip2src.items(): hl.setdefault(s, []).append(names[cid[x]])
    lab = np.array(['+'.join(sorted(hl[r])) for r in rk], dtype=object)
    return info, lab
def pairs(lab):
    ul = sorted(set(lab)); npairs = len(set(cs)) * len(ul); out = []
    for i in sorted(set(cs)):
        for j in ul:
            mk = (cs == i) & (lab == j); o = int(mk.sum())
            if o < 5: continue
            a, b = int((cs == i).sum()), int((lab == j).sum()); p = float(hypergeom.sf(o - 1, n, a, b))
            out.append({'struct_cluster': i, 'seq_label': j, 'struct_size': a, 'seq_size': b, 'overlap': o,
                        'expected': round(a * b / n, 1), 'p': p, 'p_bonf': min(1.0, p * npairs),
                        'struct_pos_median': round(float(np.median(posS[mk])), 3),
                        'arrangements': top(arr[mk]), 'superpop': top(sp[mk], 5), 'members': sorted(np.array(rk)[mk])})
    return pd.DataFrame(out).sort_values('p').reset_index(drop=True), len(ul)
# --- the UPGMA bands of step 19
cp19 = pd.read_csv(f'{O}/cluster_pairs.tsv', sep='\t'); cp19 = cp19[(cp19.k == K) & (cp19.p_bonf < ALPHA) & (cp19.overlap >= MIN_OV)]
old = {(r.struct_cluster, r.seq_cluster): set(cm19[(cm19.struct_cluster == r.struct_cluster) & (cm19.seq_cluster == r.seq_cluster)].src) for r in cp19.itertuples()}
BANDNAME = {'rest': 'Bulk', 'S1-2': 'Two', 'S1-3': 'Three', 'S1-8+S1-9': 'AB'}
# --- sensitivity over T x M
sens = []
for t in (80, 90, 95):
    for m in (5, 10, 20):
        _, lab = partition(t, m); PR, nl = pairs(lab)
        B = PR[(PR.p_bonf < ALPHA) & (PR.overlap >= MIN_OV)]
        for ok, ov in old.items():
            best = max(PR.itertuples(), key=lambda d: len(set(d.members) & ov)); sh = len(set(best.members) & ov)
            sens.append({'T': t, 'M': m, 'n_seq_labels': nl, 'n_bands': len(B), 'upgma_band': f'{ok[0]}x{ok[1]}', 'upgma_n': len(ov),
                         'ml_seq_label': best.seq_label, 'ml_n': best.overlap, 'shared': sh, 'ml_p_bonf': best.p_bonf, 'ml_is_band': bool(best.overlap >= MIN_OV and best.p_bonf < ALPHA)})
SE = pd.DataFrame(sens); SE.to_csv(f'{OUT}/ml_band_sensitivity.tsv', sep='\t', index=False)
# --- headline setting
info, lab = partition(T, M); PR, nl = pairs(lab)
pd.DataFrame([{k: v for k, v in d.items() if k != 'tips'} for d in info]).to_csv(f'{OUT}/ml_seq_clusters.tsv', sep='\t', index=False)
pd.DataFrame({'src': rk, 'struct_cluster': cs, 'seq_label': lab, 'arrangement': arr, 'superpop': sp}).to_csv(f'{OUT}/ml_hap_labels.tsv', sep='\t', index=False)
BD = PR[(PR.p_bonf < ALPHA) & (PR.overlap >= MIN_OV)].copy()
BD['band'] = BD.seq_label.map(BANDNAME)
assert BD.band.notna().all() and BD.band.is_unique, BD[['seq_label', 'overlap']]   # which bands appear is reported, and the text's claims are checked in step 40
print('bands', BD[['band', 'struct_cluster', 'seq_label', 'overlap', 'expected', 'p_bonf']].to_string() if 'expected' in BD else BD[['band', 'struct_cluster', 'seq_label', 'overlap', 'p_bonf']].to_string())
PR = PR.merge(BD[['struct_cluster', 'seq_label', 'band']], how='left')
PR.drop(columns='members').to_csv(f'{OUT}/ml_cluster_pairs.tsv', sep='\t', index=False)
MEM = pd.DataFrame([{'band': r.band, 'struct_cluster': r.struct_cluster, 'seq_label': r.seq_label, 'src': s,
                     'arrangement': L.loc[s, 'arrangement'] if s in L.index else 'NA', 'superpop': L.loc[s, 'superpop'] if s in L.index else 'reference',
                     'SMN1_lineages': L.loc[s, 'composition'] if s in L.index else 'NA'} for r in BD.itertuples() for s in r.members])
MEM.to_csv(f'{OUT}/ml_cluster_pairs_members.tsv', sep='\t', index=False)
# --- figure: ML tree pruned to these SMN1 copies (left) against the structure UPGMA (right), both trees' leaves facing the middle; the pruned tree keeps
# the topology and branch lengths of the orangutan-rooted ML tree inside the human clade
def prune(nd):
    if not nd['kids']: return dict(nd, kids=[]) if nd['name'] in sm1tips else None
    ks = [k for k in (prune(k) for k in nd['kids']) if k is not None]
    if not ks: return None
    if len(ks) == 1: k = ks[0]; k['len'] += nd['len']; return k
    return dict(nd, kids=ks)
PT = prune(HC); PT['len'] = 0.0
def size(t): t['n'] = 1 if not t['kids'] else sum(size(k) for k in t['kids']); return t['n']
size(PT)
leaves = []
def place(t, x):
    t['x'] = x + t['len']
    if not t['kids']: t['y'] = len(leaves); leaves.append(t); return t['y']
    ys = [place(k, t['x']) for k in sorted(t['kids'], key=lambda k: k['n'])]
    t['y'] = (min(ys) + max(ys)) / 2; return t['y']
place(PT, 0.0)
NT = len(leaves); ty = {l['name']: l['y'] for l in leaves}; tx = {l['name']: l['x'] for l in leaves}
assert NT == len(sm1tips)
src2tips = collections.defaultdict(list)
for x, s in tip2src.items(): src2tips[s].append(x)
# where each band sits in the drawn trees (0 = bottom, 1 = top): structure leaf, and the haplotypes' SMN1 copies
BD['seq_pos_median'] = [round(float(np.median([ty[x] / (NT - 1) for s in r.members for x in src2tips[s]])), 3) for r in BD.itertuples()]
# --- ML band against UPGMA band
u19 = cp19.set_index(['struct_cluster', 'seq_cluster'])
cmp = []
for r in BD.itertuples():
    ok = max(old, key=lambda k: len(old[k] & set(r.members))); ov = old[ok]; sh = len(ov & set(r.members)); u = u19.loc[ok]
    cmp.append({'band': r.band, 'ml_struct_cluster': r.struct_cluster, 'ml_seq_label': r.seq_label, 'ml_n': r.overlap, 'ml_expected': r.expected, 'ml_p_bonf': r.p_bonf,
                'struct_pos_median': r.struct_pos_median, 'seq_pos_median': r.seq_pos_median,
                'upgma_struct_cluster': ok[0], 'upgma_seq_cluster': ok[1], 'upgma_n': len(ov), 'upgma_expected': u.expected, 'upgma_p_bonf': u.p_bonf,
                'shared': sh, 'only_ml': len(set(r.members) - ov), 'only_upgma': len(ov - set(r.members)),
                'only_ml_arr': top([L.loc[s, 'arrangement'] for s in sorted(set(r.members) - ov)]), 'only_upgma_arr': top([L.loc[s, 'arrangement'] for s in sorted(ov - set(r.members))])})
CMP = pd.DataFrame(cmp); CMP.to_csv(f'{OUT}/ml_vs_upgma_bands.tsv', sep='\t', index=False)
BCOL = {'Bulk': '#d73027', 'Two': '#4575b4', 'Three': '#1a9850', 'AB': '#984ea3'}
BLAB = {'Bulk': 'bulk', 'Two': 'S1-2 band', 'Three': 'S1-3 band', 'AB': 'B1 + B2 block'}
SUPER_COL = {'AFR': '#8c510a', 'AMR': '#d73027', 'EAS': '#4575b4', 'EUR': '#1a9850', 'SAS': '#984ea3', 'other': '#bdbdbd', 'reference': '#000000'}
TOPARR = [a for a, _ in collections.Counter(arr).most_common(7) if a != 'NA'][:6]
ARRCOL = dict(zip(TOPARR, ['#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02']))
HGF = {'S1-1': 'S1-1', 'S1-1c': 'S1-1', 'S1-2': 'S1-2', 'S1-2c': 'S1-2', 'S1-3': 'S1-3', 'S1-4': 'S1-4', 'S1-4c': 'S1-4', 'S1-5': 'S1-5',
       'S1-6': 'S1-6', 'S1-6c': 'S1-6', 'S1-8': 'S1-8', 'S1-9': 'S1-9', 'S1-9c': 'S1-9', 'S1-9d': 'S1-9', 'S1-10': 'S1-10'}
HGCOL = {'S1-1': '#a6cee3', 'S1-2': '#1f78b4', 'S1-3': '#b2df8a', 'S1-4': '#fdbf6f', 'S1-5': '#33a02c', 'S1-6': '#fb9a99', 'S1-8': '#01665e',
         'S1-9': '#3d3d3d', 'S1-10': '#ff7f00', 'tie or other': '#cab2d6'}
hgf = lambda x: HGF.get(fam.get(x), 'tie or other')
if DRAW_FIG5:
    # arm types (step 28b) for the two inner right strips; the figure is drawn only when step 28j runs this script after 28b
    AE = pd.read_csv(f'{O}/arm_element_haps.tsv', sep='\t').set_index('src')
    ctype = {s: AE.loc[s, 'cen_type'] for s in rk if s in AE.index}; ttype = {s: AE.loc[s, 'tel_type'] for s in rk if s in AE.index}
    ATCOL = {'C1': '#525252', 'T1': '#525252', '1': '#525252', 'C0': '#fde68a', 'T0': '#fde68a', '0': '#fde68a', 'none': '#f2f2f2'}
    nct = collections.Counter(ctype.values()); ntt = collections.Counter(ttype.values())
    fig = plt.figure(figsize=(11, 13.5))
    Y0, YH = 0.15, 0.79
    aG = fig.add_axes([0.02, Y0, 0.35, YH]); aM = fig.add_axes([0.37, Y0, 0.36, YH]); aS = fig.add_axes([0.73, Y0, 0.25, YH])
    plt.rcParams['lines.linewidth'] = 0.5
    dendrogram(Zs, ax=aS, orientation='right', no_labels=True, color_threshold=0, above_threshold_color='k')
    aS.set_ylim(0, 10 * n); aS.set_yticks([]); aS.set_xticks([]); [s_.set_visible(False) for s_ in aS.spines.values()]
    aS.set_title(f'large-scale structure, {n} haplotypes\n(PANGEA bundle distance, UPGMA)', fontsize=10)
    def tdraw(t):
        for k in t['kids']:
            aG.plot([t['x'], t['x']], [t['y'], k['y']], color='k', lw=.35); aG.plot([t['x'], k['x']], [k['y'], k['y']], color='k', lw=.35); tdraw(k)
    tdraw(PT)
    XM = max(tx.values())
    for d in info:   # the supported clades that define the sequence clusters (5 or more copies of their own)
        if d['cluster'] == 'rest' or len(d['tips']) < 5: continue
        ys = sorted(ty[x] for x in d['tips'])
        aG.plot([XM * 1.04] * 2, [min(ys), max(ys)], color='#b2182b', lw=1.2)
        aG.text(XM * 1.06, (min(ys) + max(ys)) / 2, f"{d['cluster']} ({int(d['support'])})", fontsize=7.5, va='center', color='#b2182b')
    aG.set_ylim(-.5, NT - .5); aG.set_yticks([]); [aG.spines[s_].set_visible(False) for s_ in ('left', 'right', 'top')]
    aG.set_xlim(0, XM * 1.45)
    from matplotlib.ticker import FuncFormatter, MaxNLocator
    aG.xaxis.set_major_locator(MaxNLocator(3)); aG.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v*1e3:.2g}'))
    aG.set_xlabel('substitutions per kb (GTR+F+G4)', fontsize=8)
    aG.set_title(f'SMN1 copies, {NT} (gene-body ML tree, pruned;\nred: clades with UFBoot support $\\geq${T})', fontsize=10)
    memb = {s: r.band for r in BD.itertuples() for s in r.members}
    fs = lambda i: (posS[i] * (n - 1) + .5) / n; fg = lambda x: (ty[x] + .5) / NT
    for layer in (0, 1):
        for i, s in enumerate(rk):
            b = memb.get(s)
            if (b is None) != (layer == 0): continue
            for x in src2tips[s]:
                aM.plot([.15, .71], [fg(x), fs(i)], lw=.3 if b is None else .7, color=(.85, .85, .85) if b is None else BCOL[b], zorder=1 + layer)
    for i, s in enumerate(rk):
        aM.add_patch(plt.Rectangle((.94, fs(i) - .5 / n), .06, 1 / n, color=SUPER_COL.get(sp[i], '#bdbdbd'), lw=0))
        aM.add_patch(plt.Rectangle((.87, fs(i) - .5 / n), .06, 1 / n, color=ARRCOL.get(arr[i], '#e0e0e0'), lw=0))
        aM.add_patch(plt.Rectangle((.80, fs(i) - .5 / n), .06, 1 / n, color=ATCOL[ttype.get(s, 'none')], lw=0))
        aM.add_patch(plt.Rectangle((.73, fs(i) - .5 / n), .06, 1 / n, color=ATCOL[ctype.get(s, 'none')], lw=0))
    for x in sorted(sm1tips):
        aM.add_patch(plt.Rectangle((.07, fg(x) - .5 / NT), .06, 1 / NT, color=COPYCOL.get(sub.get(x, ''), '#999999'), lw=0))
        aM.add_patch(plt.Rectangle((0, fg(x) - .5 / NT), .06, 1 / NT, color=HGCOL[hgf(x)], lw=0))
    for xx, t_ in [(.97, 'ancestry'), (.90, 'arrangement'), (.83, 'telSMN arm'), (.76, 'cenSMN arm'), (.10, 'SMN1 lineage'), (.03, 'haplogroup')]:
        aM.text(xx, 1.004, t_, rotation=90, ha='center', va='bottom', fontsize=7.5)
    aM.set_xlim(0, 1); aM.set_ylim(0, 1); aM.axis('off')
    aM.text(.5, 1.004, 'each SMN1 copy (left)\nis joined to the haplotype\nwhose contig carries it (right)', ha='center', va='bottom', fontsize=8.5)
    h1 = [plt.Line2D([], [], color=BCOL[r.band], lw=2, label=f'{BLAB[r.band]}: {r.overlap} haplotypes ({r.expected} expected)') for r in BD.itertuples()]
    h1 += [plt.Line2D([], [], color=(.75, .75, .75), lw=2, label='in no band')]
    fig.legend(handles=h1, loc='upper left', bbox_to_anchor=(0.40, Y0 - 0.025), fontsize=7.5, frameon=False, title='line: one per SMN1 copy,\ncopy (left) to its haplotype (right);\ncolor: band', title_fontsize=8)
    nsp = collections.Counter(sp); narr = collections.Counter(arr[sp != 'reference'])   # arrangement counts: donor haplotypes only
    h2 = [plt.Rectangle((0, 0), 1, 1, color=SUPER_COL[v], label=f'{"ASL or no label" if v == "other" else "GRCh38/CHM13" if v == "reference" else v} ({nsp[v]})') for v in SUPER_COL if nsp.get(v)]
    h2 += [plt.Rectangle((0, 0), 1, 1, color=ARRCOL[a], label=f'"{a}" ({narr[a]})') for a in TOPARR] + [plt.Rectangle((0, 0), 1, 1, color='#e0e0e0', label='other types')]
    h2 += [plt.Rectangle((0, 0), 1, 1, color=ATCOL[k_], label=l_) for k_, l_ in [('1', f'Xd: C1 ({nct["C1"]}), T1 ({ntt["T1"]})'), ('0', f'no Xd: C0 ({nct["C0"]}), T0 ({ntt["T0"]})'), ('none', f'no arm type ({n - len(ctype)})')]]
    fig.legend(handles=h2, loc='upper left', bbox_to_anchor=(0.73, Y0 - 0.025), fontsize=7, frameon=False, ncol=2, title='right strips (haplotypes): cenSMN and\ntelSMN arm types, arrangement, ancestry', title_fontsize=8)
    nsub = collections.Counter(sub.get(x, '') for x in sorted(sm1tips)); nhg = collections.Counter(hgf(x) for x in sorted(sm1tips))
    h3 = [plt.Rectangle((0, 0), 1, 1, color=COPYCOL[v], label=f"{'B, neither clade' if v == 'B-' else v} ({nsub[v]})") for v in ['A', 'B1', 'B2', 'B-'] if nsub.get(v)]
    h3 += [plt.Rectangle((0, 0), 1, 1, color=c_, label=f'{k_} ({nhg[k_]})') for k_, c_ in HGCOL.items() if nhg.get(k_)]
    fig.legend(handles=h3, loc='upper left', bbox_to_anchor=(0.02, Y0 - 0.025), fontsize=7, frameon=False, ncol=3, title='left strips: haplogroup, lineage (copies)', title_fontsize=8)
    fig.savefig(f'{F}/tanglegram_ml.png', dpi=200); fig.savefig(f'{F}/tanglegram_ml.pdf', metadata={'CreationDate': None}); plt.close(fig)

G_ = SE.groupby(['T', 'M'])
summ = {'K': K, 'T': T, 'M': M, 'min_overlap': MIN_OV, 'alpha': ALPHA, 'haplotypes': n, 'SMN1_copies': NT, 'seq_labels': nl,
        'seq_clusters': len(info), 'haplotypes_by_SMN1_copies': {str(k): int(v) for k, v in sorted(collections.Counter(len(src2tips[s]) for s in rk).items())},
 'bands': {r.band: {'n': int(r.overlap), 'expected': r.expected, 'p_bonf': r.p_bonf} for r in BD.itertuples()},
        'sensitivity_settings': int(G_.ngroups), 'settings_all_bands': int(G_.ml_is_band.all().sum()),
        'settings_n_bands': {f'{t}/{m}': int(v) for (t, m), v in G_.n_bands.first().items()}}
json.dump(summ, open(f'{OUT}/ml_bands_summary.json', 'w'), indent=1, sort_keys=True)
pd.set_option('display.width', 250)
print(CMP.drop(columns=['only_ml_arr', 'only_upgma_arr']).to_string(index=False)); print(CMP[['band', 'only_ml_arr', 'only_upgma_arr']].to_string(index=False)); print(json.dumps(summ, indent=1))
