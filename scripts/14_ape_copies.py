# Pipeline step 14 (ape_copies); ported from explore/structure_vs_sequence/12_ape_copies.py (sha256 5617e6b3378d, commit 5e2b3bb),
# paths rewritten to results/sv/ and figures/; run from the project root after the steps before it.
"""Exploration step 12: every SMN copy in 5 Mb around each great ape's SMN gene, typed at c.840, c.*239 and the
13 fixed SMN1/SMN2 and 10 fixed A/B sites of step 11.
Inputs: data/outgroup/regions/*.fa (NCBI efetch, 5 Mb of the T2T primary assemblies NC_072402.2, NC_073253.2,
NC_073243.2, NC_071989.2 centred on the NCBI Gene SMN locus); out/copies.fa (GRCh38 SMN1 copy c0001, gene body
+ 2 kb flanks, GRCh38 chr5:70,923,087 onwards); results/sv/evolution_sites.tsv (step 13b).
Copies are found with 25-mers of the GRCh38 SMN1 gene body that occur once in it (either strand); hits are
scanned for dense runs on one strand (1-kb bins with >= 200 hits, joined across <= 5 kb, >= 5,000 hits). Each candidate (+2 kb) is
aligned to c0001 with mafft and read at the sites. Controls: GRCh38 SMN1 (c0001) and SMN2 (c0000) typed the same way.
Writes out/ape_copies.tsv and out/ape_copy_sites.tsv."""
import glob, os, subprocess, tempfile, collections
import pandas as pd
O = 'results/sv'
K = 25; G0 = 70923087
def readfa(p):
    out = {}; n = None
    for l in open(p):
        l = l.strip()
        if not l: continue
        if l[0] == '>': n = l[1:].split()[0]; out[n] = []
        else: out[n].append(l.upper())
    return {k: ''.join(v) for k, v in out.items()}
RC = str.maketrans('ACGTN', 'TGCAN')
rc = lambda s: s.translate(RC)[::-1]
cp = readfa(f'{O}/copies.fa')
H1 = cp['c0001']; H2 = cp['c0000']
sites = pd.read_csv(f'{O}/evolution_sites.tsv', sep='\t')
assert H1[70951946 - G0] == 'C' and H1[70952674 - G0] == 'G', 'c0001 coordinate frame wrong'
for r in sites.itertuples(): assert H1[r.grch38 - G0] == r.SMN1 if r.kind == 'S' else H1[r.grch38 - G0] in (r.A, r.B), r   # frame check; the GRCh38 copy (lineage A) need not carry A's majority base
body = H1[2000:len(H1) - 2000]
cnt = collections.Counter()
for i in range(len(body) - K + 1):
    k = body[i:i + K]; cnt[k] += 1; cnt[rc(k)] += 1
km = {}
for i in range(len(body) - K + 1):
    k = body[i:i + K]
    if cnt[k] == 1 and 'N' not in k: km[k] = (i, '+'); km[rc(k)] = (i, '-')
print('unique body 25-mers', len(km) // 2)
def find_copies(seq, bin_=1000, dens=200, min_hits=5000):
    """Dense runs of body 25-mer hits on one strand: 1-kb bins with >= `dens` hits, joined across <= 5 kb,
    kept with >= `min_hits` hits. Scattered hits (shared repeats such as Alus) fall below the density."""
    import numpy as np
    hits = []
    for j in range(len(seq) - K + 1):
        h = km.get(seq[j:j + K])
        if h: hits.append((j, h[0], h[1]))
    out = []
    for st in '+-':
        hs = [h for h in hits if h[2] == st]
        if not hs: continue
        b = np.bincount([h[0] // bin_ for h in hs], minlength=len(seq) // bin_ + 1)
        dense = np.where(b >= dens)[0]
        runs = []
        for x in dense:
            if runs and x - runs[-1][1] <= 5: runs[-1][1] = x
            else: runs.append([x, x])
        for r0, r1 in runs:
            lo, hi = r0 * bin_, (r1 + 1) * bin_
            c = [h for h in hs if lo <= h[0] < hi]
            if len(c) < min_hits: continue
            offs = sorted(x[1] for x in c)
            out.append(dict(strand=st, start=c[0][0], end=c[-1][0] + K, hits=len(c),
                            body_from=offs[0], body_to=offs[-1] + K, body_cov=round((offs[-1] + K - offs[0]) / len(body), 3),
                            kmer_share=round(len(set(offs)) / len(km) * 2, 3)))
    return sorted(out, key=lambda d: d['start'])
def align_read(q, label):
    with tempfile.TemporaryDirectory() as td:
        fa = os.path.join(td, 'in.fa'); open(fa, 'w').write(f'>ref\n{H1}\n>q\n{q}\n')
        r = subprocess.run(['mafft', '--quiet', '--thread', '4', fa], capture_output=True, text=True, check=True)
    al = {}; n = None
    for l in r.stdout.splitlines():
        if l.startswith('>'): n = l[1:]; al[n] = []
        else: al[n].append(l.upper())
    a, b = ''.join(al['ref']), ''.join(al['q'])
    pos = {}; k = -1
    for i, ch in enumerate(a):
        if ch != '-': k += 1; pos[k] = i
    row = dict(copy=label)
    for s in sites.itertuples():
        row[f'{s.kind}_{s.grch38}'] = b[pos[s.grch38 - G0]]
    return row
rows, srows = [], []
for lab, q in [('GRCh38_SMN1_c0001', H1), ('GRCh38_SMN2_c0000', H2)]:
    srows.append(align_read(q, lab)); rows.append(dict(ape='human GRCh38', copy=lab))
for f in sorted(glob.glob('data/outgroup/regions/*.fa')):
    stem, rng = os.path.basename(f)[:-3].rsplit('_', 1); ape, acc = stem.split('_', 1)
    off = int(rng.split('-')[0])
    seq = list(readfa(f).values())[0]
    cps = find_copies(seq)
    print(ape, acc, len(seq), 'candidates', len(cps))
    for n, c in enumerate(cps):
        lab = f'{ape}_{n + 1}'
        lo = max(0, c['start'] - 2500); hi = min(len(seq), c['end'] + 2500)
        q = seq[lo:hi]; q = q if c['strand'] == '+' else rc(q)
        rows.append(dict(ape=ape, copy=lab, chrom=acc, chr_start=off + c['start'], chr_end=off + c['end'] - 1, **c))
        if c['body_cov'] >= 0.5: srows.append(align_read(q, lab))
C = pd.DataFrame(rows); Sx = pd.DataFrame(srows)
sc = [c for c in Sx.columns if c.startswith('S_')]; ac = [c for c in Sx.columns if c.startswith('AB_')]
sm = sites.set_index(sites.kind + '_' + sites.grch38.astype(str))
Sx['c840'] = Sx['S_70951946']; Sx['c239'] = Sx['S_70952674']
Sx['S_smn1'] = [sum(r[c] == sm.loc[c, 'SMN1'] for c in sc) for _, r in Sx.iterrows()]
Sx['S_smn2'] = [sum(r[c] == sm.loc[c, 'SMN2'] for c in sc) for _, r in Sx.iterrows()]
Sx['AB_A'] = [sum(r[c] == sm.loc[c, 'A'] for c in ac) for _, r in Sx.iterrows()]
Sx['AB_B'] = [sum(r[c] == sm.loc[c, 'B'] for c in ac) for _, r in Sx.iterrows()]
C = C.merge(Sx, on='copy', how='left')
C.to_csv(f'{O}/ape_copies.tsv', sep='\t', index=False)
print(C.drop(columns=sc + ac).to_string())
c0 = C.set_index('copy')
assert c0.loc['GRCh38_SMN1_c0001', 'c840'] == 'C' and c0.loc['GRCh38_SMN2_c0000', 'c840'] == 'T', 'control typing failed'
assert c0.loc['GRCh38_SMN2_c0000', 'c239'] == 'A' and c0.loc['GRCh38_SMN2_c0000', 'S_smn2'] == len(sc), 'SMN2 control not SMN2 at all S sites'
# checks for the reported result: exactly one SMN copy in each 5-Mb ape region, it is the copy step 09 used
# (data/outgroup/<ape>_SMN_<acc>_<start>-<end>*.fa), and it carries SMN1's base at c.840 and c.*239
for ape in ['chimp', 'bonobo', 'gorilla', 'orangutan']:
    r = C[C.ape == ape]; assert len(r) == 1, (ape, len(r))
    f = glob.glob(f'data/outgroup/{ape}_SMN_*.fa')[0]; a, b = f.split('_')[-1].replace('_rc', '').replace('.fa', '').split('_')[-1].split('-') if False else os.path.basename(f).replace('_rc.fa', '').replace('.fa', '').split('_')[-1].split('-')
    assert int(a) <= r.chr_start.iloc[0] and r.chr_end.iloc[0] <= int(b), (ape, f)
    assert r.c840.iloc[0] == 'C' and r.c239.iloc[0] == 'G', ape
print('checks passed: one SMN copy per ape region, the step-09 copy, SMN1 bases at c.840 and c.*239')
