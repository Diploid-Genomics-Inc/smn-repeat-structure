# Pipeline step 27 (arm_contents); ported from explore/structure_vs_sequence/24_arm_contents.py (sha256 aa12cb39c13f, commit 5e2b3bb),
# paths rewritten to results/sv/ and figures/; run from the project root after the steps before it.
"""Step 24: what the two arms of the inverted family-A pair hold (resolved haplotypes, Fig. 10 set).
(1) SMN copies per arm against arm length; (2) which copy class / SMN1 lineage fills each arm's SMN slot
(left arm = A+, centromeric = cenSMN position; right arm = A-, telomeric = telSMN position);
(3) nested repeats inside each A unit: pangea bundles occurring more than once inside the same unit
(repeat-units partitions a window, so a nested repeat shows only at the bundle level);
(4) genes (results/copies/features.tsv) between the axis and the copy, and beyond the copy;
(5) each recurring bundle pair: spanning the SMN copy (the nested inverted repeat) or beyond it (a direct duplication),
with the genes overlapping each occurrence.
Inputs share one frame: window coordinates in GRCh38 orientation (checked below on copy/SERF1 placement)."""
import pandas as pd, numpy as np, json
from scipy.stats import spearmanr, fisher_exact, mannwhitneyu
O = 'results/sv/'
RU = pd.read_csv('results/decomp/smn.repeat_units.tsv', sep='\t', comment='#')
B = pd.read_csv('results/decomp/smn.bed', sep='\t', comment='#', header=None, names=['row', 'bgn', 'end', 'name'])
f = B.name.str.split(':', expand=True); B['bid'] = f[0].astype(int); B['dir'] = f[2].astype(int)
CU = pd.read_csv('results/decomp/copy_units.tsv', sep='\t')
F = pd.read_csv('results/copies/features.tsv', sep='\t')
AX = pd.read_csv(O + 'palindrome_arms.tsv', sep='\t')
KEEP = set(pd.read_csv(O + 'palindrome_diagram_rows.tsv', sep='\t').src)   # step 22's pair haplotypes (a copy on each arm)
AX = AX[AX.src.isin(KEEP)]
PL = pd.read_csv(O + 'palindrome_copies_lineage.tsv', sep='\t')
CL = pd.read_csv(O + 'copy_labels.tsv', sep='\t')
sp = CL.dropna(subset=['superpop']).groupby('src').superpop.first() if 'src' in CL else None
pal = set(AX.src)
arms, units, genes, cps, pairs = [], [], [], [], []
for src, g in RU[RU.row.isin(pal)].groupby('row'):
    g = g.sort_values('bgn').reset_index(drop=True)
    for i in range(len(g) - 1):
        a, b = g.iloc[i], g.iloc[i + 1]
        if not (a.family == 'A' and b.family == 'A' and a.strand == '+' and b.strand == '-'): continue
        axis = (a.end + b.bgn) / 2
        for side, u in (('left', a), ('right', b)):
            dist = (lambda x: axis - x) if side == 'left' else (lambda x: x - axis)
            s = B[(B.row == src) & (B.bgn >= u.bgn) & (B.end <= u.end)].sort_values('bgn')
            vc = s.bid.value_counts(); rb = vc[vc > 1].index; r = s[s.bid.isin(rb)]
            fi = r.drop_duplicates('bid', keep='first').set_index('bid'); la = r.drop_duplicates('bid', keep='last').set_index('bid')
            inv = float((fi.dir != la.loc[fi.index].dir).mean()) if len(rb) else np.nan
            oc = float(np.corrcoef(fi.bgn, la.loc[fi.index].bgn)[0, 1]) if len(rb) > 2 else np.nan
            c = CU[(CU.src == src) & (CU.mid_fwd >= u.bgn) & (CU.mid_fwd <= u.end)]
            counted = c[~c.token.str[0].isin(['t', 'f'])]
            cd = [dist(m) / 1e3 for m in counted.mid_fwd]
            for m, tok in zip(c.mid_fwd, c.token):
                hit = s[(s.bgn <= m) & (s.end >= m)]
                cps.append(dict(src=src, side=side, token=tok, d_kb=dist(m) / 1e3, in_bundle=len(hit) > 0, in_nested=bool(len(hit) and hit.bid.isin(rb).any())))
            rep_d = [dist((x.bgn + x.end) / 2) / 1e3 for x in r.itertuples()]
            fg = F[(F.src == src) & (F.kind == 'gene') & (F.win_start >= u.bgn) & (F.win_end <= u.end)]
            def genes_at(lo, hi): return ','.join(sorted(set(fg[(fg.win_start < hi) & (fg.win_end > lo)].gene)))
            ccd = np.median(cd) if cd else np.nan
            for bid in vc[vc == 2].index:
                x2 = s[s.bid == bid].copy(); x2['d'] = [dist((q.bgn + q.end) / 2) / 1e3 for q in x2.itertuples()]; x2 = x2.sort_values('d')
                i0, i1 = x2.iloc[0], x2.iloc[1]
                # kind: 'spans' = the two occurrences lie either side of the arm's SMN copy; 'beyond' = both beyond it
                kind = 'spans' if i0.d < ccd < i1.d else ('beyond' if i0.d > ccd else 'inside')
                pairs.append(dict(src=src, side=side, bid=bid, inner_kb=i0.d, outer_kb=i1.d, inverted=bool(i0.dir != i1.dir), span_kb=u.span_bp / 1e3,
                                  copy_d_kb=ccd, kind=kind, occ_bp=int(i0.end - i0.bgn), genes_inner=genes_at(i0.bgn, i0.end), genes_outer=genes_at(i1.bgn, i1.end)))
            units.append(dict(src=src, side=side, span_kb=u.span_bp / 1e3, n_items=len(c), n_counted=len(counted), copy_d_kb=np.median(cd) if cd else np.nan,
                              n_bundle_ids=s.bid.nunique(), n_nested_ids=len(rb), nested_bp=int((r.end - r.bgn).sum()), nested_inverted_share=inv, nested_order_corr=oc,
                              nested_min_d_kb=min(rep_d) if rep_d else np.nan, nested_max_d_kb=max(rep_d) if rep_d else np.nan))
            ff = F[(F.src == src) & (F.win_start >= u.bgn) & (F.win_end <= u.end)]
            for x in ff.itertuples():
                d = dist((x.win_start + x.win_end) / 2) / 1e3
                genes.append(dict(src=src, side=side, span_kb=u.span_bp / 1e3, gene=x.gene, kind=x.kind, d_kb=d, beyond_copy=bool(cd) and d > max(cd)))
U = pd.DataFrame(units); G = pd.DataFrame(genes); C = pd.DataFrame(cps)
U['long'] = U.span_kb >= 500; U['nested'] = U.n_nested_ids > 0
if sp is not None: U['superpop'] = U.src.map(sp)
NP = pd.DataFrame(pairs); NP.to_csv(O + 'arm_nested_pairs.tsv', sep='\t', index=False)
U.to_csv(O + 'arm_units.tsv', sep='\t', index=False); G.to_csv(O + 'arm_genes.tsv', sep='\t', index=False); C.to_csv(O + 'arm_copies.tsv', sep='\t', index=False)
# slot occupancy: class / lineage of the counted copy on each arm
PL['slot'] = np.where(PL.gene == 'SMN1', 'SMN1 ' + PL.lineage.fillna('?'), PL.gene)
PL['side'] = np.where(PL.arm.str.startswith('left'), 'left', 'right')
if sp is not None: PL['afr'] = np.where(PL.src.map(sp) == 'AFR', 'AFR', 'other')
slot = pd.crosstab([PL.side, PL.afr] if sp is not None else PL.side, PL.slot)
slot.to_csv(O + 'arm_slots.tsv', sep='\t')
# summaries
S = {}
S['n_haps'] = int(U.src.nunique()); S['n_units'] = len(U)
S['counted_per_arm'] = {int(k): int(v) for k, v in U.n_counted.value_counts().items()}
S['copy_d_short'] = float(U[~U.long].copy_d_kb.median()); S['copy_d_long'] = float(U[U.long].copy_d_kb.median())
ct = pd.crosstab(U.long, U.nested); S['long_nested'] = ct.to_dict()
S['fisher_long_nested_p'] = fisher_exact(ct.values)[1]
S['rho_span_nestedbp'] = spearmanr(U.span_kb, U.nested_bp)[0]
N = U[U.nested]
S['nested_ids_median'] = float(N.n_nested_ids.median()); S['nested_bp_median'] = float(N.nested_bp.median())
S['nested_inverted_share_median'] = float(N.nested_inverted_share.median()); S['nested_order_corr_median'] = float(N.nested_order_corr.median())
S['nested_all_inverted'] = int((N.nested_inverted_share >= 0.9).sum())
S['nested_min_d_median'] = float(N.nested_min_d_kb.median()); S['nested_max_d_median'] = float(N.nested_max_d_kb.median())
S['bundle_ids_short'] = float(U[~U.long].n_bundle_ids.median()); S['bundle_ids_long'] = float(U[U.long].n_bundle_ids.median())
S['copies_in_nested'] = int(C.in_nested.sum()); S['copies_total'] = len(C)
S['nested_by_side'] = pd.crosstab(U.side, U.nested).to_dict()
if sp is not None: S['nested_by_afr'] = pd.crosstab(U.superpop == 'AFR', U.nested).to_dict()
Gg = G[G.kind == 'gene']
beyond = Gg[Gg.beyond_copy].merge(U[['src', 'side', 'long']], on=['src', 'side'])
S['genes_beyond_copy_per_unit'] = beyond.groupby(['long', 'gene']).size().div(U.groupby('long').size(), level=0).round(2).to_dict()
inner = Gg[~Gg.beyond_copy].merge(U[['src', 'side', 'long']], on=['src', 'side'])
S['genes_inner_per_unit'] = inner.groupby(['long', 'gene']).size().div(U.groupby('long').size(), level=0).round(2).to_dict()
json.dump({k: (v if not isinstance(v, dict) else {str(a): (b if not isinstance(b, dict) else {str(c): d for c, d in b.items()}) for a, b in v.items()}) for k, v in S.items()}, open(O + 'arm_summary.json', 'w'), indent=1, default=float)
print(json.dumps(json.load(open(O + 'arm_summary.json')), indent=1)); print(slot)
assert S['n_haps'] == len(KEEP)
# checks for the reported result
# was an assert under the clinical decomposition; reported now, text claims are checked in step 40
if not ((U.n_counted == 1).all()): print('NOTE (step 27):', "(U.n_counted == 1).all(), 'an arm holds more or fewer than one counted SMN copy'")
# was an assert under the clinical decomposition; reported now, text claims are checked in step 40
# The 'nested' checks below test the long-arm reading of the finer decomposition's unit call, which the note no longer uses
# (superseded by the arm elements of step 28b; no number reads arm_nested_pairs.tsv). On the coarse decomposition they print a NOTE only.
if not (C.in_nested.sum() == 0): print('NOTE (step 27):', "C.in_nested.sum() == 0, 'an SMN copy lies inside a nested (repeated) bundle'")
# was an assert under the clinical decomposition; reported now, text claims are checked in step 40
if not (set(NP.src + NP.side) <= set(U[U.nested].src + U[U.nested].side)): print('NOTE (step 27):', 'set(NP.src + NP.side) <= set(U[U.nested].src + U[U.nested].side)')
Gg2 = G[G.kind == 'gene']
if not (Gg2.d_kb.max() < 560 and (Gg2[Gg2.gene != 'SERF1B'].d_kb > 300).all()): print('NOTE (step 27): genes on the arms at', round(Gg2.d_kb.min(), 1), '-', round(Gg2.d_kb.max(), 1), 'kb')
# checks for the two kinds of recurring pair (Section 4.12)
# was an assert under the clinical decomposition; reported now, text claims are checked in step 40
if not ((NP[NP.kind == 'spans'].inverted).all() and (~NP[NP.kind == 'beyond'].inverted).all()): print('NOTE (step 27):', "(NP[NP.kind == 'spans'].inverted).all() and (~NP[NP.kind == 'beyond'].inverted).all()")
# was an assert under the clinical decomposition; reported now, text claims are checked in step 40
if not ((NP[NP.kind == 'spans'][['genes_inner', 'genes_outer']].fillna('') == '').all().all()): print('NOTE (step 27): a gene in the nested inverted repeat')
# was an assert under the clinical decomposition; reported now, text claims are checked in step 40
if not (NP[NP.kind == 'beyond'].genes_inner.str.contains('NAIP').all() and NP[NP.kind == 'beyond'].genes_outer.str.contains('NAIP').all()): print('NOTE (step 27):', "NP[NP.kind == 'beyond'].genes_inner.str.contains('NAIP').all() and NP[NP.kind == 'beyond'].genes_outer.str.contains('NAI")
print(NP.kind.value_counts().to_dict())
