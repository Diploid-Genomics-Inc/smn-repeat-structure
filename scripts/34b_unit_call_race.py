"""NOTE (2026-11-15): this step reads the CLINICAL decomposition (results/bundles/), on purpose. Under the study's
decomposition (step 00) the two windows get the same unit string (C+ B- A+ A- D- B+), so the race it replays
belongs to the earlier decomposition only.
Pipeline step 34b (supplement, Fig. S1c): why PANGEA's repeat-units calls HG03816 pat 'B- + short A+' and HG01099 mat
one long A+ although both carry the same arm elements (step 28b). Replays the chaining rule (scripts/chainrule.py) on
results/bundles/smn.bed. Both windows reach a junction match J in the A core (~343 kb from the axis) from two
predecessors: the inner pairing (the A core near the axis) and the outer route (the mirror chain over Xi and the arrays,
which must skip E plus a copy of the core end on the telSMN arm to reach J). The DP keeps the higher of
best(predecessor) - step penalty; if the outer route loses, it becomes a chain of its own and the unit is cut.
Writes results/sv/unit_call_race.json and figures/unit_call_race.png/.pdf."""
import sys, json, collections
import pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
sys.path.insert(0, 'scripts'); from smnlib import FAMCOL, ELCOL, famlab
from chainrule import dp, trace, DIV
O = 'results/sv'
H = [('HG01099_mat_hprc_r2_v1.0.1', 'HG01099 mat'), ('HG03816_pat_hprc_r2_v1.0.1', 'HG03816 pat')]
rows = collections.defaultdict(list)
for l in open('results/bundles/smn.bed'):
    if l.startswith('#') or not l.strip(): continue
    f = l.rstrip('\n').split('\t'); b = f[3].split(':')
    rows[f[0]].append((int(f[1]), int(f[2]), int(b[0]), int(b[2])))
RU = pd.read_csv('results/bundles/smn.repeat_units.tsv', sep='\t', comment='#')
EB = pd.read_csv(f'{O}/arm_element_blocks.tsv', sep='\t')
out, fig = {}, plt.figure(figsize=(12, 6.9))
for k, (src, name) in enumerate(H):
    seg = sorted(rows[src], key=lambda s: (s[0], s[1]))
    A = RU[(RU.row == src) & (RU.family == 'A')].sort_values('bgn')
    pr = [(a, b) for a, b in zip(A.itertuples(), list(A.itertuples())[1:]) if a.strand == '+' and b.strand == '-']
    ax_ = (pr[0][0].end + pr[0][1].bgn) // 2
    L = lambda i: (ax_-seg[i][0])/1e3; R = lambda j: (seg[j][1]-ax_)/1e3
    best, cands = dp(seg)
    # J: the first match of the chain to the axis that has both an inner (right arm < 420 kb) and an outer (> 450 kb) predecessor
    J = None
    for e in sorted(best, key=lambda e: -e[0]):
        if L(e[0]) < 330 or L(e[0]) > 360: continue
        c = cands[e]
        if any(R(x[0][1]) > 450 for x in c) and any(R(x[0][1]) < 420 and x[2] < 1000 for x in c): J = e; break
    c = cands[J]
    inner = max((x for x in c if R(x[0][1]) < 420), key=lambda x: x[3]); outer = max((x for x in c if R(x[0][1]) > 450), key=lambda x: x[3])
    po = trace(best, outer[0]); pi = trace(best, inner[0])
    # the chain from J inward to the axis (same in both windows): best end near the axis whose traceback passes J
    pa = next(trace(best, e) for e in sorted(best, key=lambda e: -best[e][0]) if L(e[0]) < 60 and 0 < R(e[1]) < 60 and J in trace(best, e))
    pa = pa[pa.index(J):]
    steps = []
    for a, b in zip(po, po[1:]):
        gi = max(seg[b[0]][0]-seg[a[0]][1], 0); gj = max(seg[a[1]][0]-seg[b[1]][1], 0)
        steps.append(dict(l0=L(a[0]), r0=R(a[1]), l1=L(b[0]), r1=R(b[1]), gap_cen=gi/1e3, gap_tel=gj/1e3, pen=abs(gi-gj)//DIV/1e3))
    tand = max(steps, key=lambda s: s['pen'])
    won = 'outer' if best[J][1] == outer[0] else 'inner'
    out[name] = dict(src=src, axis=ax_, J=[round(L(J[0]), 1), round(R(J[1]), 1)],
        inner=dict(pred=[round(L(inner[0][0]), 1), round(R(inner[0][1]), 1)], best=inner[1]/1e3, pen=inner[2]/1e3),
        outer=dict(pred=[round(L(outer[0][0]), 1), round(R(outer[0][1]), 1)], best=outer[1]/1e3, pen=outer[2]/1e3,
                   gap_cen=outer[4]/1e3, gap_tel=outer[5]/1e3, matches=len(po), start=[round(L(po[0][0]), 1), round(R(po[0][1]), 1)]),
        tandem_step=tand, winner=won,
        unit_strings={s: ' '.join(f'{r.family}{r.strand}' for r in RU[(RU.row == src)].sort_values('bgn').itertuples()
                                  if (ax_-r.end if s == 'cen' else r.bgn-ax_) > -50000 and (ax_-r.bgn if s == 'cen' else r.end-ax_) < 900000) for s in ('cen', 'tel')})
    # drawing: x = cenSMN arm, y = telSMN arm, kb from the axis
    g = fig.add_axes([0.06+k*0.495, 0.17, 0.43, 0.72])
    mx = [(L(i), R(j)) for (i, j) in best if seg[i][1] <= ax_+5000 and seg[j][0] >= ax_-5000 and 0 < L(i) < 900 and 0 < R(j) < 900]
    g.scatter([x for x, _ in mx], [y for _, y in mx], s=4, c='#bbbbbb', lw=0, label='mirror match (same bundle, opposite strand)')
    for p, col, lab in ((po, '#e67e22', 'outer route'), (pi, '#2b6cb0', 'inner pairing')):
        g.plot([L(i) for i, _ in p], [R(j) for _, j in p], '-o', c=col, ms=2.5, lw=1.6, label=lab)
    g.plot([L(outer[0][0]), L(J[0])], [R(outer[0][1]), R(J[1])], ':', c='#e67e22', lw=1.6)
    g.plot([L(i) for i, _ in pa], [R(j) for _, j in pa], '-', c='#333333', lw=1.2, label='chain from J to the axis (both windows)')
    g.plot([L(J[0])], [R(J[1])], 'k*', ms=11, label='junction J')
    g.plot([tand['l0'], tand['l1']], [tand['r0'], tand['r1']], c='#c0392b', lw=3)
    g.annotate(f"array step: gaps {tand['gap_cen']:.0f} vs {tand['gap_tel']:.0f} kb\npenalty {tand['pen']:.1f} kb", ((tand['l0']+tand['l1'])/2, (tand['r0']+tand['r1'])/2),
               xytext=(-25, 0), textcoords='offset points', ha='right', va='center', fontsize=10, color='#c0392b', arrowprops=dict(arrowstyle='-', color='#c0392b', lw=0.6))
    g.annotate(f"E skip: penalty {outer[2]/1e3:.1f} kb", ((L(outer[0][0])+L(J[0]))/2, (R(outer[0][1])+R(J[1]))/2), xytext=(25, 5), textcoords='offset points', fontsize=10, color='#e67e22')
    vo, vi = outer[1]/1e3-outer[2]/1e3, inner[1]/1e3-inner[2]/1e3
    g.set_title(f"{name}: outer {outer[1]/1e3:.1f} − {outer[2]/1e3:.1f} = {vo:.1f} {'>' if vo > vi else '<'} inner {vi:.1f} kb\n"
                f"→ {'one chain: one long R$_A$+ unit' if won == 'outer' else 'outer route cut off: R$_B$− / R$_B$+ units beyond a short R$_A$ pair'}", fontsize=11, loc='left')
    g.set_xlim(-65, 900); g.set_ylim(-65, 900); g.set_aspect('equal', adjustable='box')
    g.set_xlabel('cenSMN arm, kb from the axis', fontsize=10); g.set_ylabel('telSMN arm, kb from the axis', fontsize=10)
    for side, (u0, w) in (('cen', (-30, 14)), ('tel', (-30, 14))):
        for r in EB[(EB.src == src) & (EB.side == side)].itertuples():
            rect = Rectangle((r.a0, u0), r.a1-r.a0, w) if side == 'cen' else Rectangle((u0, r.a0), w, r.a1-r.a0)
            rect.set_color(ELCOL.get(r.el, '#e0e0e0')); g.add_patch(rect)
        for r in RU[RU.row == src].itertuples():
            a0, a1 = (ax_-r.end, ax_-r.bgn) if side == 'cen' else (r.bgn-ax_, r.end-ax_)
            if a1 < 0 or a0 > 900000: continue
            a0, a1 = max(a0, 0)/1e3, min(a1, 900000)/1e3
            rect = Rectangle((a0, -55), a1-a0, 14) if side == 'cen' else Rectangle((-55, a0), 14, a1-a0)
            rect.set_color(FAMCOL.get(r.family, '#cccccc')); rect.set_alpha(1 if r.strand == '+' else 0.55); g.add_patch(rect)
            t = famlab(r.family, r.strand)
            if side == 'cen': g.text((a0+a1)/2, -48, t, ha='center', va='center', fontsize=7)
            else: g.text(-48, (a0+a1)/2, t, ha='center', va='center', fontsize=7, rotation=90)
    if k == 0: hl = g.get_legend_handles_labels()
fig.legend(*hl, loc='lower center', bbox_to_anchor=(0.5, 0.0), ncol=3, fontsize=9.5, frameon=False)
fig.savefig('figures/unit_call_race.png', dpi=200); fig.savefig('figures/unit_call_race.pdf', metadata={'CreationDate': None})
json.dump(out, open(f'{O}/unit_call_race.json', 'w'), indent=1, default=float)
for n, d in out.items(): print(n, d['J'], d['inner'], {k_: v for k_, v in d['outer'].items()}, d['tandem_step'], d['winner'], d['unit_strings'])
