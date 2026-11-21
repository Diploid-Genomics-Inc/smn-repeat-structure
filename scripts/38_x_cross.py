"""Pipeline step 38 (supplement): X seen as an 'X' in a self dot plot of the cenSMN arm (researcher 2026-11-18:
'it is cool X include the center for the inverted repeat that is the center of the X shape in the dot plot').
Promoted from explore/x_shape/x01_x_cross.py (same method). For each of the 270 palindrome haplotypes
(results/sv/arm_element_haps.tsv), the cenSMN arm 300-900 kb from the axis (axis = midpoint of the R_A+ R_A- junction,
results/decomp/smn.repeat_units.tsv) is dotted against itself (pangeatk.get_shmmr_dots w24 k48 r2 min_span 12, as
Fig. S10). An inverted repeat centred at c puts its dots on x+y = 2c, so the centre (x+y)/2 of off-diagonal dots
(|x-y| > 10 kb) peaks at c. Per arm: the 5-kb bin of centres with most dots; for C1 arms its distance from S
(bundle 20 beyond 300 kb on the cen arm, results/sv/arm_element_segments.tsv).
Out: results/sv/x_cross.tsv, results/sv/x_cross.json, figures/x_cross.png/.pdf."""
import json, pandas as pd, numpy as np, pangeatk as ptk, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
O = 'results/sv/'
H = pd.read_csv(O + 'arm_element_haps.tsv', sep='\t')
SEG = pd.read_csv(O + 'arm_element_segments.tsv', sep='\t')
RU = pd.read_csv('results/decomp/smn.repeat_units.tsv', sep='\t', comment='#')
rd = lambda s: ''.join(l.strip() for l in open(f'data/derived/windows/{s}.fa') if not l.startswith('>'))
P = (24, 48, 2, 12); LO, HI = 300_000, 900_000; EX = ('HG01099_mat_hprc_r2_v1.0.1', 'HG01099_pat_hprc_r2_v1.0.1')
def axis_of(src):
    u = RU[RU.row == src].sort_values('bgn').reset_index(drop=True); lab = list(u.family.astype(str) + u.strand)
    k = [j for j in range(len(u) - 1) if lab[j] == 'A+' and lab[j + 1] == 'A-']
    assert len(k) >= 1, src
    return (u.end[k[0]] + u.bgn[k[0] + 1]) / 2
rows, keep = [], {}
for src, ct in zip(H.src, H.cen_type):
    s = rd(src); ax = axis_of(src)
    a, b = max(0, int(ax - HI)), int(ax - LO)
    if b - a < 100_000: continue
    x, y = ptk.get_shmmr_dots(s[a:b].encode(), s[a:b].encode(), *P)
    x = np.asarray(x); y = np.asarray(y); off = np.abs(x - y) > 10_000
    d = (ax - (a + (x[off] + y[off]) / 2)) / 1e3
    h, e = np.histogram(d, bins=np.arange(LO / 1e3, HI / 1e3 + 5, 5)); i = int(h.argmax())
    sg = SEG[(SEG.src == src) & (SEG.side == 'cen') & (SEG.el == 'S') & (SEG.a0 > LO / 1e3)]
    sd = float(((sg.a0 + sg.a1) / 2).mean()) if len(sg) else np.nan
    rows.append(dict(src=src, cen_type=ct, arm_kb=round((ax - a) / 1e3, 1), n_off=int(off.sum()), peak_kb=e[i] + 2.5,
                     peak_share=round(h[i] / max(1, h.sum()), 4), S_kb=round(sd, 2) if sd == sd else np.nan,
                     peak_minus_S=round(e[i] + 2.5 - sd, 2) if sd == sd else np.nan))
    if src in EX: keep[src] = (x, y, a, ax, sd)
D = pd.DataFrame(rows); D.to_csv(O + 'x_cross.tsv', sep='\t', index=False)
c1 = D[D.cen_type == 'C1']; ws = c1[c1.S_kb.notna()]; Di = D.set_index('src')
S = dict(shimmer=dict(w=24, k=48, r=2, min_span=12), lo_kb=LO // 1000, hi_kb=HI // 1000, arms=len(D),
         c1_arms=len(c1), c1_with_S=len(ws), c1_peak_le20=int((ws.peak_minus_S.abs() <= 20).sum()),
         c1_peak_le10=int((ws.peak_minus_S.abs() <= 10).sum()), c1_peak_le5=int((ws.peak_minus_S.abs() <= 5).sum()),
         c1_without_S=int(c1.S_kb.isna().sum()),
         examples={k: dict(type=Di.cen_type[k], S_kb=None if keep[k][4] != keep[k][4] else round(keep[k][4], 1), arm_kb=float(Di.arm_kb[k])) for k in EX})
json.dump(S, open(O + 'x_cross.json', 'w'), indent=1); print(json.dumps(S, indent=1))
fig, axs = plt.subplots(1, 3, figsize=(18, 6.4))
for k, (src, lab) in enumerate([(EX[0], 'C1, carries X'), (EX[1], 'C0, no X')]):
    x, y, a, ax_, sd = keep[src]; A = axs[k]
    gx = (ax_ - (a + x)) / 1e3; gy = (ax_ - (a + y)) / 1e3
    A.scatter(gx, gy, s=.05, c='k', lw=0, rasterized=True)
    if sd == sd:
        A.axvline(sd, color='#d62728', lw=.7, ls='--'); A.axhline(sd, color='#d62728', lw=.7, ls='--')
    A.set_xlim(HI / 1e3, LO / 1e3); A.set_ylim(HI / 1e3, LO / 1e3); A.set_aspect('equal')
    A.set_xlabel('cenSMN arm, kb from axis'); A.set_ylabel('cenSMN arm, kb from axis')
    A.set_title(f"({'ab'[k]}) {' '.join(src.split('_')[:2])}, {lab}" + ('; dashed: S' if sd == sd else ''), loc='left', fontsize=11)
A = axs[2]; A.hist(ws.peak_minus_S, bins=40, color='#333333')
A.set_title(f'(c) {len(ws)} C1 arms: center of the cross minus S', loc='left', fontsize=11)
A.set_xlabel('strongest center of off-diagonal dots - position of S (kb)'); A.set_ylabel('arms')
fig.tight_layout(); fig.savefig('figures/x_cross.png', dpi=150); fig.savefig('figures/x_cross.pdf', metadata={'CreationDate': None})
