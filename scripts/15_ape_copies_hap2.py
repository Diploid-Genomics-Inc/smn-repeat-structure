# Pipeline step 15 (ape_copies_hap2); ported from explore/structure_vs_sequence/15_ape_copies_hap2.py (sha256 1e7a2c0ce789, commit 5e2b3bb),
# paths rewritten to results/sv/ and figures/; run from the project root after the steps before it.
"""Exploration step 15: SMN copies in each great ape's OTHER haplotype (the second T2T haplotype assembly), typed as
in step 12. Step 12 searched the primary assemblies only and found one SMN copy per ape.
Inputs: data/outgroup/regions_hap2/*.fa, NCBI efetch of the SMN-bearing chromosome of the second haplotype:
  chimp CM054415.2 (NHGRI_mPanTro3-v2.0_alt, chr 4), orangutan CM054660.2 (NHGRI_mPonAbe1-v2.0_alt, chr 4),
  bonobo CM054464.2 (mPanPan1 mat/pat haplotype whose chr 4 length differs from the primary's; the other, CM054488.2,
  has the primary NC_073253.2's length), gorilla CM054601.2 (mGorGor1 haplotype whose chr 19 differs from the primary's;
  CM054577.2 has the primary NC_073243.2's length). Windows are 8 Mb (18 Mb for gorilla, whose chr 19 is 7.7 Mb
  shorter) around the primary's SMN coordinate. Accessions and lengths: data/outgroup/assembly/*.json.
Candidates are found twice: step 12's density rule (1-kb bins >= 200 hits, >= 5,000 hits) and a relaxed rule
(>= 50 hits per bin, >= 1,000 hits) to catch partial copies. Writes out/ape_copies_hap2.tsv."""

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
for f in sorted(glob.glob('data/outgroup/regions_hap2/*.fa')):
    b = os.path.basename(f)[:-3]; ape, hapl, acc, rng = b.split('_')
    off = int(rng.split('-')[0])
    seq = list(readfa(f).values())[0]
    strict = find_copies(seq)
    relaxed = find_copies(seq, dens=50, min_hits=1000)
    print(ape, hapl, acc, len(seq), 'strict', len(strict), 'relaxed', len(relaxed), flush=True)
    for rule, cps in [('strict', strict), ('relaxed', relaxed)]:
        for n, c in enumerate(cps):
            lab = f'{ape}_{hapl}_{rule}_{n + 1}'
            rows.append(dict(ape=ape, haplotype=hapl, rule=rule, copy=lab, chrom=acc, chr_start=off + c['start'],
                             chr_end=off + c['end'] - 1, **c))
            if rule == 'strict' and c['body_cov'] >= 0.5:
                lo = max(0, c['start'] - 2500); hi = min(len(seq), c['end'] + 2500)
                q = seq[lo:hi]; q = q if c['strand'] == '+' else rc(q)
                srows.append(align_read(q, lab))
C = pd.DataFrame(rows); Sx = pd.DataFrame(srows)
sc = [c for c in Sx.columns if c.startswith('S_')]; ac = [c for c in Sx.columns if c.startswith('AB_')]
sm = sites.set_index(sites.kind + '_' + sites.grch38.astype(str))
Sx['c840'] = Sx['S_70951946']; Sx['c239'] = Sx['S_70952674']
Sx['S_smn1'] = [sum(r[c] == sm.loc[c, 'SMN1'] for c in sc) for _, r in Sx.iterrows()]
Sx['S_smn2'] = [sum(r[c] == sm.loc[c, 'SMN2'] for c in sc) for _, r in Sx.iterrows()]
Sx['AB_A'] = [sum(r[c] == sm.loc[c, 'A'] for c in ac) for _, r in Sx.iterrows()]
Sx['AB_B'] = [sum(r[c] == sm.loc[c, 'B'] for c in ac) for _, r in Sx.iterrows()]
C = C.merge(Sx, on='copy', how='left')
C.to_csv(f'{O}/ape_copies_hap2.tsv', sep='\t', index=False)
print(C.drop(columns=sc + ac).to_string())
# checks for the reported result: one SMN copy (strict rule) in each ape's other haplotype, SMN1 bases at c.840 and
# c.*239; relaxed-rule extras must overlap the strict copy (orangutan: a short inverted stretch inside its own copy)
st = C[C.rule == 'strict']
for ape in ['chimp', 'bonobo', 'gorilla', 'orangutan']:
    r = st[st.ape == ape]; assert len(r) == 1, (ape, len(r))
    assert r.c840.iloc[0] == 'C' and r.c239.iloc[0] == 'G', ape
    for x in C[(C.ape == ape) & (C.rule == 'relaxed')].itertuples():
        assert x.chr_start <= r.chr_end.iloc[0] and x.chr_end >= r.chr_start.iloc[0], ('relaxed copy outside', ape, x.copy)
print('checks passed: one SMN copy per other haplotype, SMN1 bases at c.840 and c.*239')
# whole-genome check through NCBI Gene annotation (primary assemblies): every live SMN-named gene record of the four
# apes (data/outgroup/assembly/esummary_gene_smn_all.json; query in esearch_gene_smn_all.json)
import json
g = json.load(open('data/outgroup/assembly/esummary_gene_smn_all.json'))['result']
recs = []
for u in g['uids']:
    r = g[u]; gi = (r.get('genomicinfo') or [{}])[0]
    s0, s1 = int(gi.get('chrstart', 0)), int(gi.get('chrstop', 0))
    recs.append(dict(gene_id=u, organism=r['organism']['scientificname'], symbol=r['name'], description=r['description'],
                     chromosome=r.get('chromosome'), accession=gi.get('chraccver'), start=min(s0, s1), end=max(s0, s1),
                     length=abs(s1 - s0) + 1, live=str(r.get('status', '')) in ('', '0')))
R = pd.DataFrame(recs); R.to_csv(f'{O}/ape_smn_gene_records.tsv', sep='\t', index=False)
L = R[R.live & R.description.str.startswith('survival')]
print(L.drop(columns=['live']).to_string())
full = L[L.length > 20000]
assert sorted(full.organism) == sorted(['Pan troglodytes', 'Pan paniscus', 'Gorilla gorilla gorilla', 'Pongo abelii']), full
assert set(full.accession) == {'NC_072402.2', 'NC_073253.2', 'NC_073243.2', 'NC_071989.2'}, full
other = L[(L.length <= 20000) & ~L.symbol.str.startswith('SMNDC')]
assert (other.length < 2000).all(), other
print(f'checks passed: one full-length SMN gene record per ape (the step-09 copy); {len(other)} other live SMN-like '
      f'records, all under 2 kb ({other.length.min()}-{other.length.max()} bp)')
