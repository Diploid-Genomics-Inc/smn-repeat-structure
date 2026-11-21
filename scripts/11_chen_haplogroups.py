# Pipeline step 11 (chen_haplogroups); ported from explore/structure_vs_sequence/08_chen_haplogroups.py (sha256 ba54d5964448, commit 5e2b3bb),
# paths rewritten to results/sv/ and figures/; run from the project root after the steps before it.
"""Exploration step 08: match the copy-tree lineages A and B to Chen et al. 2023's SMN1 haplogroups.
Two independent routes. (1) Copy level: every copy is genotyped, in out/copies.aln.fa, at the haplogroup-defining
variants of Chen's Table S6 (papers/supp/chen2023/mmc2.xlsx; GRCh38 position_ref_alt), placed through the GRCh38 SMN1
and SMN2 rows of the alignment; each copy gets the haplogroup whose variant set it matches best (Jaccard over the
Table S6 variants of that gene's haplogroups). (2) Person level: for people in Chen's Table S3, our A and B copy counts
against the haplogroup names Paraphase gave them (via tables/published_vs_assembly.tsv).
Writes out/chen_hg_copies.tsv, out/chen_hg_crosstab.tsv, out/chen_hg_people.tsv. Nothing here feeds the note."""
import re, numpy as np, pandas as pd
O = 'results/sv'
X = 'papers/supp/chen2023/mmc2.xlsx'
names, seqs = [], []
for l in open(f'{O}/copies.aln.fa'):
    l = l.strip()
    if l[0] == '>': names.append(l[1:]); seqs.append([])
    else: seqs[-1].append(l.upper())
A = np.array([np.frombuffer(''.join(s).encode(), dtype=np.uint8) for s in seqs])
t = pd.read_csv(f'{O}/copy_tips.tsv', sep='\t').set_index('tip').loc[names].reset_index()
cc = pd.read_csv('results/copies/copy_calls.tsv', sep='\t').set_index('copy_id')
def refrow(cls):
    i = np.where((t.copy_class == cls) & t.src.str.startswith('GCA_000001405'))[0]; assert len(i) == 1; i = i[0]
    r = cc.loc[t.copy_id[i]]; assert r.strand == '+'
    start = int(r.ctg_start) - 2000                       # sequences carry 2 kb flanks (01_gene_tree.sh)
    cols = np.where(A[i] != ord('-'))[0]
    return i, start, cols
R = {'SMN1': refrow('SMN1'), 'SMN2': refrow('SMN2')}
s6 = pd.read_excel(X, sheet_name='TableS6', header=1)
var = s6.iloc[:, 0].astype(str); hgs = [c for c in s6.columns[1:]]
rows = []
NS6 = len(var); NSNV = 0
for v, (_, r) in zip(var, s6.iterrows()):
    p, ref, alt = v.split('_'); p = int(p)
    if len(ref) != 1 or len(alt) != 1: continue          # SNVs only
    NSNV += 1
    for loc, (i, start, cols) in R.items():
        k = p - 1 - start                                 # Table S6 positions are 1-based; ctg_start is 0-based
        if 0 <= k < len(cols):
            col = cols[k]; ok = chr(A[i, col]) == ref
            rows.append(dict(var=v, locus=loc, col=col, ref=ref, alt=alt, ref_ok=ok, hgs=[h for h in hgs if r[h] == 'x']))
V = pd.DataFrame(rows)
print('SNVs placed', len(V), 'of', len(var), '; ref base agrees', int(V.ref_ok.sum()))
assert V.ref_ok.all(), 'every placed Table S6 reference base must match GRCh38'
V = V.reset_index(drop=True)
pd.DataFrame([dict(s6_variants=NS6, s6_snvs=NSNV, placed=len(V), placed_smn1=int((V.locus == 'SMN1').sum()), placed_smn2=int((V.locus == 'SMN2').sum()))]).to_csv(f'{O}/chen_hg_variants.tsv', sep='\t', index=False)
ml = pd.read_csv(f'{O}/marker_lineage.tsv', sep='\t').set_index('tip')
cls = t.copy_class.replace({'hyb_2_1': 'h21', 'hyb_1_2': 'h12'})
lin = [('B' if ml.loc[x, 'in_marker_subclade'] else 'A') if (c == 'SMN1' and x in ml.index) else c for x, c in zip(t.tip, cls)]
S1 = [h for h in hgs if h.startswith('S1-')]
GT = np.array([[A[j, c] == ord(a) for c, a in zip(V.col, V.alt)] for j in range(len(A))])   # copy x variant: alt present
CALL = np.array([[A[j, c] in b'ACGT' for c in V.col] for j in range(len(A))])
hgset = {h: np.array([h in x for x in V.hgs]) for h in hgs}
univ1 = np.any([hgset[h] for h in S1], axis=0)
out = []
for j in range(len(A)):
    if lin[j] not in ('A', 'B'): continue
    g = GT[j] & univ1; sc = {}
    for h in S1:
        m = hgset[h] & CALL[j]
        inter = (g & m).sum(); union = (g | m).sum()
        sc[h] = inter / union if union else 0.0
    bs = max(sc.values()); tied = sorted(h for h, s in sc.items() if s == bs)
    second = max([s for s in sc.values() if s < bs], default=0.0)
    out.append(dict(tip=t.tip[j], copy_id=t.copy_id[j], lineage=lin[j], hg=tied[0] if len(tied) == 1 else '|'.join(tied),
                    n_tied=len(tied), jaccard=round(bs, 3), margin=round(bs - second, 3), n_alt=int(g.sum())))
C = pd.DataFrame(out)
fam = lambda h: re.sub(r'(?<=\d)[cd]$', '', h)
C['family'] = [next(iter(s)) if len(s) == 1 else 'tie' for s in ({fam(x) for x in h.split('|')} for h in C.hg)]
C.to_csv(f'{O}/chen_hg_copies.tsv', sep='\t', index=False)
X2 = pd.crosstab(C.family, C.lineage); X2.to_csv(f'{O}/chen_hg_crosstab.tsv', sep='\t')
print(pd.crosstab(C.family, C.lineage).to_string()); print(C[C.n_tied > 1].groupby(['lineage', 'hg']).size().to_string()); print(C.groupby('lineage').jaccard.describe().round(2))
# person level
pv = pd.read_csv('tables/published_vs_assembly.tsv', sep='\t')
lh = pd.read_csv(f'{O}/lineage_haplotypes.tsv', sep='\t'); per = lh.groupby('sample')[['nA', 'nB']].sum()
P = []
for r in pv.itertuples():
    if not isinstance(r.allele1, str) or r.sample not in per.index: continue
    hg = [x for a in (r.allele1, r.allele2) for x in re.split('[,+]', str(a)) if x.startswith('S1-')]
    n = lambda f: sum(1 for x in hg if re.fullmatch(f, x))
    P.append(dict(sample=r.sample, population=r.population, nA=int(per.loc[r.sample, 'nA']), nB=int(per.loc[r.sample, 'nB']),
                  n_S1_9=n(r'S1-9[cd]?'), n_S1_6=n(r'S1-6c?'), n_S1_8=n('S1-8'), n_S1_10=n('S1-10'), n_S1=len(hg), alleles=f'{r.allele1} / {r.allele2}'))
P = pd.DataFrame(P); P.to_csv(f'{O}/chen_hg_people.tsv', sep='\t', index=False)
print('people', len(P), 'nB==S1-9*', int((P.nB == P.n_S1_9).sum()), 'nB==S1-9*+S1-6*', int((P.nB == P.n_S1_9 + P.n_S1_6).sum()))
