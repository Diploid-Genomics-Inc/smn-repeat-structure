"""Pipeline step 13c: each great ape's base at every S and A/B site (step 13b), read from the ape's raw primary-assembly
sequence rather than from the multiple alignment (step 09's mafft --add --keeplength), which can misplace an ape base
next to a homopolymer or where the ape differs. Researcher 2026-11-15: 'which alignment is right at the 5 ape sites'.
For each site the 101-bp GRCh38 context (the site +-50 bp, GRCh38 chr5 read from the HPRC archive) is aligned locally (Smith-Waterman,
match 2, mismatch -3, gap 5) to the ape's 5-Mb primary region (data/outgroup/regions/, the 400 kb around its SMN copy),
both strands, from every 14-mer seed; the best-scoring alignment is kept. The read is accepted when it covers the
site, has identity >= 0.85 over >= 60 aligned columns, and has no gap within 3 columns of the site; otherwise the ape's
base there is unresolved. Writes results/sv/ape_site_reads.tsv."""
import glob, re
import pandas as pd
O = 'results/sv'; F = 50; MIN_ID = 0.85; MIN_LEN = 60; GAPWIN = 3
def rd(f): return ''.join(l.strip() for l in open(f) if not l.startswith('>')).upper()
import os, subprocess
_out = subprocess.run(['agc-rs', 'get', os.path.join(os.environ['STRATA_DATA_HPRC'], 'hprc_r2_full.agcrs'),
                       'GCA_000001405.15_GRCh38_no_alt_analysis_set.PanSN/GRCh38#0#chr5:70900000-70980000'], capture_output=True, text=True, check=True).stdout
G = ''.join(l.strip() for l in _out.splitlines() if not l.startswith('>')).upper(); G0 = 70900001   # agc-rs ranges are 0-based
assert len(G) == 80000
RC = str.maketrans('ACGTN', 'TGCAN')
APEF = {'og_chimp': 'chimp', 'og_bonobo': 'bonobo', 'og_gorilla': 'gorilla', 'og_orangutan': 'orangutan'}
regs = {}
for f in sorted(glob.glob('data/outgroup/regions/*.fa')):
    s = rd(f); regs[f.split('/')[-1].split('_')[0]] = s[2350000:2750000]   # the SMN copy sits at the region's centre (2.5 Mb)
assert set(regs) == set(APEF.values())
def sw(a, b):
    n, m = len(a), len(b); H = [[0] * (m + 1) for _ in range(n + 1)]; T = [[0] * (m + 1) for _ in range(n + 1)]; best = (0, 0, 0)
    for i in range(1, n + 1):
        Hi, Hp, Ti, ai = H[i], H[i - 1], T[i], a[i - 1]
        for j in range(1, m + 1):
            d = Hp[j - 1] + (2 if ai == b[j - 1] else -3); u = Hp[j] - 5; l = Hi[j - 1] - 5
            v = max(0, d, u, l); Hi[j] = v; Ti[j] = (0 if v == 0 else 1 if v == d else 2 if v == u else 3)
            if v > best[0]: best = (v, i, j)
    _, i, j = best; A, B = [], []
    while i > 0 and j > 0 and T[i][j]:
        t = T[i][j]
        if t == 1: A.append(a[i - 1]); B.append(b[j - 1]); i -= 1; j -= 1
        elif t == 2: A.append(a[i - 1]); B.append('-'); i -= 1
        else: A.append('-'); B.append(b[j - 1]); j -= 1
    return best[0], ''.join(A[::-1]), ''.join(B[::-1]), i
E_ = pd.read_csv(f'{O}/evolution_sites.tsv', sep='\t')
assert all(G[int(r.grch38) - G0] in 'ACGT' for r in E_.itertuples())
rows = []
for r in E_.itertuples():
    o = int(r.grch38) - G0; q = G[o - F:o + F + 1]
    for col_, ape in APEF.items():
        s = regs[ape]; hits = []
        for ss in (s, s.translate(RC)[::-1]):
            seeds = set()
            for k in range(0, len(q) - 14, 3):
                for m in re.finditer(q[k:k + 14], ss): seeds.add(m.start() - k)
            reps = []
            for c in sorted(seeds):
                if not reps or c - reps[-1] > 40: reps.append(c)
            for c in reps[:30]:
                lo = max(0, c - 40); hits.append(sw(q, ss[lo:c + len(q) + 40]))
        row = dict(kind=r.kind, grch38=int(r.grch38), ape=col_, msa=getattr(r, col_), read='', score=0, second=0, ident=0.0, alen=0, gap_near=False, accepted=False)
        if hits:
            hits.sort(key=lambda h: -h[0]); sc, A, B, i0 = hits[0]
            qi, col = i0, None
            for k, ch in enumerate(A):
                if ch != '-':
                    if qi == F: col = k; break
                    qi += 1
            ident = sum(x == y for x, y in zip(A, B)) / len(A)
            row.update(score=sc, second=hits[1][0] if len(hits) > 1 else 0, ident=round(ident, 3), alen=len(A))
            if col is not None:
                gap = '-' in A[max(0, col - GAPWIN):col + GAPWIN + 1] or '-' in B[max(0, col - GAPWIN):col + GAPWIN + 1]
                row.update(read=B[col], gap_near=gap, accepted=(B[col] in 'ACGT' and not gap and ident >= MIN_ID and len(A) >= MIN_LEN),
                           context=A[max(0, col - 8):col] + '[' + A[col] + ']' + A[col + 1:col + 9] + ' / ' + B[max(0, col - 8):col] + '[' + B[col] + ']' + B[col + 1:col + 9])
        row['checked'] = row['read'] if row['accepted'] else '?'
        rows.append(row)
D = pd.DataFrame(rows); D.to_csv(f'{O}/ape_site_reads.tsv', sep='\t', index=False)
print(len(D), 'ape-site reads;', int(D.accepted.sum()), 'accepted;', int((D.accepted & (D.read != D.msa)).sum()), 'accepted reads differ from the multiple alignment')
print(D[(~D.accepted) | (D.read != D.msa)][['kind', 'grch38', 'ape', 'msa', 'read', 'score', 'second', 'ident', 'alen', 'gap_near', 'accepted']].to_string())
