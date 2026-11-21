"""Pipeline step 35 (from exploration p01, 2026-11-14; researcher: 'primate structure', then 'add to supplemental, both dot plot' and 'a supplement figure and a sentence in §2.7'): do the great apes carry the human 5q13 repeat
structure (the A+ A- palindrome, X = Xi + S + Xd, E)?
Inputs: data/outgroup/regions/*.fa and regions_hap2/*.fa (NCBI efetch of the ape T2T assemblies, 5-18 Mb around SMN;
data/SOURCES.md), results/sv/ape_copies.tsv / ape_copies_hap2.tsv (SMN copy position in each file, steps 14/15),
the human reference haplotype HG01099 mat (Fig. 1a, C1|T0) window data/derived/windows/, its axis
(palindrome_diagram_rows.tsv) and arm-element blocks (arm_element_blocks.tsv, step 28b).
For each of the 8 ape haplotypes: a 3 Mb stretch centred on its SMN copy (clipped to the file).
(1) self dot plot (pangeatk.get_shmmr_dots, w24 k48 r2 min_span 12, the settings of Fig. S9);
(2) large-duplication content: dots with |x-y| > 50 kb, per Mb, against the human window;
(3) human vs ape: for each human element block, the share of its shimmer dots that land in the ape stretch, and the
number of distinct ape loci (50 kb clusters) they land in.
(4) the fold: for the A core of each human arm, the densest 600-kb ape stretch of its hits, the share of hits in it,
the overlap of the two arms' stretches, and Spearman's rho of ape position against distance from the human axis;
(5) S: hits per human shimmer position in S against the unique flank (Fl) on the same ape.
Writes results/sv/ape_struct_summary.tsv, ape_element_hits.tsv, ape_fold.tsv, ape_struct.json,
figures/ape_self_dotplots.png/.pdf, figures/ape_vs_human.png/.pdf."""
import glob, os, re, json
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
import pangeatk as ptk, sys
sys.path.insert(0, 'scripts'); from smnlib import ELCOL
O = 'results/sv'; S = 'results/sv'; F = 'figures'
W, K, R, MS = 24, 48, 2, 12
HALF = 1_500_000
def readfa(p): return ''.join(l.strip() for l in open(p) if not l.startswith('>')).upper()
dots = lambda a, b: tuple(np.asarray(v) for v in ptk.get_shmmr_dots(a.encode(), b.encode(), W, K, R, MS))
c1 = pd.read_csv(f'{S}/ape_copies.tsv', sep='\t'); c1 = c1[c1.ape != 'human GRCh38']
c2 = pd.read_csv(f'{S}/ape_copies_hap2.tsv', sep='\t'); c2 = c2[c2.rule == 'strict']
apes = []
for r in c1.itertuples():
    f = glob.glob(f'data/outgroup/regions/{r.ape}_{r.chrom}_*.fa'); assert len(f) == 1, f
    apes.append((f'{r.ape} primary', f[0], int(r.start)))
for r in c2.itertuples():
    f = glob.glob(f'data/outgroup/regions_hap2/{r.ape}_{r.haplotype}_{r.chrom}_*.fa'); assert len(f) == 1, f
    apes.append((f'{r.ape} {r.haplotype}', f[0], int(r.start)))
assert len(apes) == 8
HS = 'HG01099_mat_hprc_r2_v1.0.1'
hum = readfa(f'data/derived/windows/{HS}.fa')
ax = float(pd.read_csv(f'{S}/palindrome_diagram_rows.tsv', sep='\t').set_index('src').axis[HS])
BLK = pd.read_csv(f'{S}/arm_element_blocks.tsv', sep='\t'); BLK = BLK[BLK.src == HS].copy()
BLK['h0'] = np.where(BLK.side == 'cen', ax - BLK.a1 * 1e3, ax + BLK.a0 * 1e3)
BLK['h1'] = np.where(BLK.side == 'cen', ax - BLK.a0 * 1e3, ax + BLK.a1 * 1e3)
rows, hits = [], []
hx, hy = dots(hum, hum)
hoff = int((np.abs(hx - hy) > 50_000).sum())
rows.append(dict(haplotype='human HG01099 mat', length=len(hum), copy_at=None, dots=len(hx), far_dots=hoff, far_per_mb=hoff / len(hum) * 1e6))
selfd, cross, seqs = {}, {}, {}
for name, f, cp in apes:
    s = readfa(f); a0 = max(0, cp - HALF); s = s[a0: cp + HALF]; seqs[name] = (s, cp - a0)
    x, y = dots(s, s); off = int((np.abs(x - y) > 50_000).sum())
    rows.append(dict(haplotype=name, length=len(s), copy_at=cp - a0, dots=len(x), far_dots=off, far_per_mb=off / len(s) * 1e6))
    selfd[name] = (x, y)
    qx, qy = dots(hum, s); cross[name] = (qx, qy)   # x = human window, y = ape stretch
    for b in BLK.itertuples():
        m = (qx >= b.h0) & (qx < b.h1)
        nh = int(((hx >= b.h0) & (hx < b.h1) & (hx == hy)).sum())   # human shimmer positions in the block (diagonal)
        ys = np.sort(qy[m]); loci = int((np.diff(ys) > 50_000).sum() + 1) if len(ys) else 0
        hits.append(dict(ape=name, side=b.side, el=b.el, h_kb=f'{b.h0/1e3:.0f}-{b.h1/1e3:.0f}', block_kb=round((b.h1 - b.h0) / 1e3, 1),
                         human_positions=nh, ape_hits=int(m.sum()), hits_per_position=round(m.sum() / nh, 2) if nh else None, ape_loci=loci))
    print(name, len(s), 'far dots/Mb', round(off / len(s) * 1e6), flush=True)
HT = pd.DataFrame(hits); HT.to_csv(f'{O}/ape_element_hits.tsv', sep='\t', index=False)
T = pd.DataFrame(rows); T.to_csv(f'{O}/ape_struct_summary.tsv', sep='\t', index=False)
from scipy.stats import spearmanr
WIN = 600_000
fold = []
for name, _, _ in apes:
    qx, qy = cross[name]; st = {}
    for side in ('cen', 'tel'):
        b = BLK[(BLK.side == side) & (BLK.el == 'A')].sort_values('block_kb' if 'block_kb' in BLK else 'a0').iloc[0] if False else BLK[(BLK.side == side) & (BLK.el == 'A')].assign(L=lambda d: d.h1 - d.h0).sort_values('L').iloc[-1]
        m = (qx >= b.h0) & (qx < b.h1); ys = np.sort(qy[m]); xs = qx[m]
        j = np.searchsorted(ys, ys + WIN); cnt = j - np.arange(len(ys)); i = int(np.argmax(cnt))
        lo, hi = int(ys[i]), int(ys[min(j[i], len(ys)) - 1])
        inw = (qy[m] >= lo) & (qy[m] <= hi)
        dist = np.abs(xs[inw] - ax)
        rho = float(spearmanr(dist, qy[m][inw])[0]) if inw.sum() > 10 else float('nan')
        st[side] = dict(lo=lo, hi=hi, share=float(inw.mean()), n=int(m.sum()), rho=rho)
    ov = max(0, min(st['cen']['hi'], st['tel']['hi']) - max(st['cen']['lo'], st['tel']['lo']))
    un = max(st['cen']['hi'], st['tel']['hi']) - min(st['cen']['lo'], st['tel']['lo'])
    e = HT[HT.ape == name].groupby('el').apply(lambda d: d.ape_hits.sum() / d.human_positions.sum())
    fold.append(dict(ape=name, cen_lo_kb=st['cen']['lo'] / 1e3, cen_hi_kb=st['cen']['hi'] / 1e3, cen_share=round(st['cen']['share'], 3), cen_rho=round(st['cen']['rho'], 3),
                     tel_lo_kb=st['tel']['lo'] / 1e3, tel_hi_kb=st['tel']['hi'] / 1e3, tel_share=round(st['tel']['share'], 3), tel_rho=round(st['tel']['rho'], 3),
                     overlap=round(ov / un, 3) if un else 0, copy_kb=seqs[name][1] / 1e3,
                     S_rate=round(float(e.get('S', np.nan)), 3), Fl_rate=round(float(e.get('Fl', np.nan)), 3), S_over_Fl=round(float(e.get('S', np.nan) / e.get('Fl', np.nan)), 3),
                     Xi_rate=round(float(e.get('Xi', np.nan)), 3), Xd_rate=round(float(e.get('Xd', np.nan)), 3), A_rate=round(float(e.get('A', np.nan)), 3)))
FT = pd.DataFrame(fold); FT.to_csv(f'{O}/ape_fold.tsv', sep='\t', index=False); print(FT.to_string(index=False))
json.dump(dict(human=HS, far_per_mb_human=float(T.far_per_mb.iloc[0]), far_per_mb_apes=[float(T.far_per_mb.iloc[1:].min()), float(T.far_per_mb.iloc[1:].max())], win_kb=WIN / 1e3), open(f'{O}/ape_struct.json', 'w'), indent=1)
print(T.to_string(index=False))
print(HT.groupby(['el']).hits_per_position.describe()[['min', '50%', 'max']].to_string())
fig, axs = plt.subplots(2, 4, figsize=(18, 9.5))
for a_, (name, _, _) in zip(axs.flat, apes):
    x, y = selfd[name]; L = len(seqs[name][0]) / 1e3
    a_.scatter(x / 1e3, y / 1e3, s=.2, c='k', lw=0, rasterized=True); c = seqs[name][1] / 1e3
    a_.scatter([c], [c], s=60, facecolors='none', edgecolors='#e66101', lw=1.5, zorder=4)
    a_.set_xlim(0, L); a_.set_ylim(L, 0); a_.set_aspect('equal'); a_.set_title(name, fontsize=10)
    a_.set_xlabel('kb'); a_.tick_params(labelsize=7)
fig.suptitle('Great-ape SMN regions against themselves (3 Mb centered on the one SMN copy, orange circle)', fontsize=12)
fig.tight_layout(); fig.savefig(f'{F}/ape_self_dotplots.png', dpi=150); fig.savefig(f'{F}/ape_self_dotplots.pdf', metadata={'CreationDate': None}); plt.close(fig)
fig, axs = plt.subplots(2, 4, figsize=(18, 9.5))
for a_, (name, _, _) in zip(axs.flat, apes):
    x, y = cross[name]; L = len(seqs[name][0]) / 1e3
    a_.scatter(x / 1e3, y / 1e3, s=.2, c='k', lw=0, rasterized=True)
    for b in BLK.itertuples():
        a_.axvspan(b.h0 / 1e3, b.h1 / 1e3, ymin=0, ymax=.03, color=ELCOL.get(b.el, '#ccc'), lw=0)
    a_.axvline(ax / 1e3, color='r', ls='--', lw=.7); a_.axhline(seqs[name][1] / 1e3, color='#e66101', lw=.7)
    a_.set_xlim(0, len(hum) / 1e3); a_.set_ylim(L, 0); a_.set_title(name, fontsize=10); a_.tick_params(labelsize=7)
    a_.set_xlabel('human HG01099 mat window (kb)'); a_.set_ylabel('ape (kb)', fontsize=8)
fig.suptitle('Human HG01099 mat (x; arm elements along the bottom, axis dashed red) against each ape (y; SMN copy orange)', fontsize=12)
fig.tight_layout(); fig.savefig(f'{F}/ape_vs_human.png', dpi=150); fig.savefig(f'{F}/ape_vs_human.pdf', metadata={'CreationDate': None}); plt.close(fig)
