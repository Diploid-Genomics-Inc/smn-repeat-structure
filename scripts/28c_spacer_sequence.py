"""Pipeline step 28c: the palindrome spacer S checked by sequence, not by bundles.
Each S occurrence of step 28b (results/sv/arm_element_spacers.tsv: a run of bundle-15/33 segments; site = axis, cen, tel)
is cut from its window (data/derived/windows, step 00), 21.3 kb centred on its bundle.
(1) Orientation by sequence: 21-mer containment of a reference S (the axis S of REF) in the occurrence and in its
    reverse complement; the larger decides. Compared with the bundle order along the window (15..33 or 33..15).
(2) Variant: each non-axis occurrence's edit distance per kb to its own haplotype's axis S; 'Sa' (like the axis S) below
    VAR_CUT per kb, 'St' above (the distances fall in two groups with none between 2 and 3 per kb; checked).
(3) Identity: each occurrence, put in the reference orientation, is aligned end to end to the others of the same
    haplotype (edlib, global edit distance), and to the same site of a fixed sample of other haplotypes.
Out: results/sv/spacer_sequence.tsv (per occurrence), spacer_pairs.tsv (per pair), spacer_summary.json."""
import json, random, itertools, numpy as np, pandas as pd, edlib
O = 'results/sv/'; WD = 'data/derived/windows/'; K = 21; REF = 'HG01258_pat_hprc_r2_v1.0.1'
SO = pd.read_csv(O + 'arm_element_spacers.tsv', sep='\t')
COMP = str.maketrans('ACGTN', 'TGCAN'); rc = lambda s: s.translate(COMP)[::-1]
_seq = {}
def window(src):
    if src not in _seq: _seq[src] = ''.join(l.strip() for l in open(WD + src + '.fa') if not l.startswith('>'))
    return _seq[src]
km = lambda s: {s[i:i + K] for i in range(0, len(s) - K + 1)}
FLANK = 5000                        # S's core = the span of its bundle (bundle 20, ~5.8 kb); each core is aligned inside the other
                                    # occurrence's core plus FLANK bp each side (infix alignment), so bundle ends that differ by
                                    # a few kb between copies do not count as differences
SO['seq'] = [window(r.src)[r.bgn:r.end] for r in SO.itertuples()]
SO['ext'] = [window(r.src)[max(0, r.bgn - FLANK):r.end + FLANK] for r in SO.itertuples()]
ref = SO[(SO.src == REF) & (SO.site == 'axis')]; assert len(ref) == 1
RK = km(ref.seq.iloc[0])
fw, rv = [], []
for s in SO.seq:
    k = km(s); fw.append(len(RK & k) / len(RK)); rv.append(len(RK & km(rc(s))) / len(RK))
SO['cont_fwd'] = np.round(fw, 3); SO['cont_rev'] = np.round(rv, 3)
SO['seq_orient'] = np.where(SO.cont_fwd >= SO.cont_rev, 'ref', 'reverse')
SO['oseq'] = [s if o == 'ref' else rc(s) for s, o in zip(SO.seq, SO.seq_orient)]
SO['oext'] = [s if o == 'ref' else rc(s) for s, o in zip(SO.ext, SO.seq_orient)]
SO['len_bp'] = SO.seq.str.len()
def ident(a, b):
    r = edlib.align(a.oseq, b.oext, mode='HW', task='distance')
    return r['editDistance'], len(a.oseq)
pairs = []
for src, g in SO.groupby('src'):
    for (i, a), (j, b) in itertools.combinations(g.iterrows(), 2):
        ed, L = ident(a, b)
        pairs.append(dict(kind='within', src_a=src, site_a=a.site, src_b=src, site_b=b.site, edits=ed, len=L, per_kb=round(1e3 * ed / L, 3)))
random.seed(1); srcs = sorted(SO.src.unique()); samp = sorted(random.sample(srcs, 40))
for site in ('axis', 'cen', 'tel'):
    g = SO[(SO.site == site) & SO.src.isin(samp)].drop_duplicates('src')
    for (i, a), (j, b) in itertools.combinations(g.iterrows(), 2):
        ed, L = ident(a, b)
        pairs.append(dict(kind='between', src_a=a.src, site_a=site, src_b=b.src, site_b=site, edits=ed, len=L, per_kb=round(1e3 * ed / L, 3)))
VAR_CUT = 2.5
ax = SO[SO.site == 'axis'].drop_duplicates('src').set_index('src')   # 2 haplotypes have two S at the axis: the first is used
SO['to_axis_per_kb'] = [round(1e3 * ident(r, ax.loc[r.src])[0] / len(r.oseq), 3) if (r.site != 'axis' and r.src in ax.index) else np.nan for r in SO.itertuples()]
SO['variant'] = np.where(SO.site == 'axis', 'Sa', np.where(SO.to_axis_per_kb.isna(), '', np.where(SO.to_axis_per_kb < VAR_CUT, 'Sa', 'St')))
print('S to axis S per kb, histogram:', np.histogram(SO.to_axis_per_kb.dropna(), bins=[0, 1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 50, 1000])[0].tolist(), 'lengths', SO.len_bp.describe().round(0).to_dict())
assert not SO.to_axis_per_kb.between(2, 3, inclusive='left').any(), 'S distances to the axis S are not in two groups'
AT = pd.read_csv(O + 'arm_elements.tsv', sep='\t').set_index(['src', 'side']).type
SO['arm_type'] = [AT.get((r.src, r.site), 'axis') for r in SO.itertuples()]
SO['rank_out'] = SO.assign(ad=SO.d_kb.abs()).groupby(['src', 'site']).ad.rank().astype(int)
P = pd.DataFrame(pairs); P.to_csv(O + 'spacer_pairs.tsv', sep='\t', index=False)
SO.drop(columns=['seq', 'oseq', 'ext', 'oext']).to_csv(O + 'spacer_sequence.tsv', sep='\t', index=False)
S = {'n_occ': len(SO), 'by_site': SO.site.value_counts().to_dict(), 'len_bp_median': float(SO.len_bp.median()),
     'len_bp_q05_q95': SO.len_bp.quantile([.05, .95]).tolist(),
     'cont_best_min': float(SO[['cont_fwd', 'cont_rev']].max(axis=1).min()), 'cont_best_median': float(SO[['cont_fwd', 'cont_rev']].max(axis=1).median()),
     'cont_other_max': float(SO[['cont_fwd', 'cont_rev']].min(axis=1).max()),
     'orient_vs_bundle': {f'{a}|{b}|{c}': int(n) for (a, b, c), n in SO.groupby(['site', 'window_order', 'seq_orient']).size().items()},
     'variant': {f'{a}|{b}|{c}': int(n) for (a, b, c), n in SO[SO.variant != ''].groupby(['site', 'arm_type', 'variant']).size().items()},
     'to_axis_by_variant': {v: dict(n=len(g), median=float(g.to_axis_per_kb.median()), q05=float(g.to_axis_per_kb.quantile(.05)), q95=float(g.to_axis_per_kb.quantile(.95))) for v, g in SO[(SO.site != 'axis') & SO.to_axis_per_kb.notna()].groupby('variant')},
     'pairs': {f'{k}:{a}-{b}': dict(n=len(g), median=float(g.per_kb.median()), q95=float(g.per_kb.quantile(.95)))
               for (k, a, b), g in P.assign(a=P[['site_a', 'site_b']].min(axis=1), b=P[['site_a', 'site_b']].max(axis=1)).groupby(['kind', 'a', 'b'])}}
json.dump(S, open(O + 'spacer_summary.json', 'w'), indent=1, default=str)
print(json.dumps(S, indent=1, default=str))
print(pd.crosstab([SO.site, SO.window_order], SO.seq_orient))
