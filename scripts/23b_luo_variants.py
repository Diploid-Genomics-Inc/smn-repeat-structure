"""Pipeline step 23b: Luo et al. 2014's six SMN1 variants of the duplication haplotype (their Table 4,
data/papers_tables/luo2014_table4.tsv, GRCh37) placed on our copy alignment and counted by lineage (A, B1, B2, SMN2).
GRCh37 -> GRCh38 by one constant shift, checked two ways: the shifted g.27706_27707delAT must start at c.*211 of
NM_000344.4 (from the RefSeq exons, located in the GRCh38 SMN1 copy), and every shifted reference base must equal
Luo's reference base in GRCh38. Writes results/sv/luo_variants.tsv."""
import sys, collections, numpy as np, pandas as pd
sys.path.insert(0, 'scripts'); from smnlib import read_genbank
O = 'results/sv'
names, seqs = [], []
for l in open(f'{O}/copies.aln.fa'):
    l = l.strip()
    if l.startswith('>'): names.append(l[1:].split()[0]); seqs.append([])
    else: seqs[-1].append(l.upper())
M = np.array([np.frombuffer(''.join(s).encode(), dtype=np.uint8) for s in seqs]); idx = {n: i for i, n in enumerate(names)}
t = pd.read_csv(f'{O}/copy_tips.tsv', sep='\t'); ml = pd.read_csv(f'{O}/marker_lineage.tsv', sep='\t').set_index('tip')
ref = idx[t[(t.copy_class == 'SMN1') & t.src.str.startswith('GCA_000001405')].tip.iloc[0]]
rseq = bytes(M[ref][M[ref] != ord('-')]).decode(); nongap = M[ref] != ord('-')
col_of = {p: c for c, p in zip(np.where(nongap)[0], np.where(nongap)[0].cumsum() * 0 + (np.cumsum(nongap)[nongap] - 1))}
OFF = 70923086                      # GRCh38 position of the reference copy's first base minus 1 (as in steps 09 and 23)
g = read_genbank('data/refs_seq/smn1_smn2_refseq.gb')['NM_000344.4']
ex8 = g['exons'][8]; p8 = rseq.find(g['seq'][ex8[0] - 1:ex8[1]]); assert p8 >= 0
cds_end = g['cds'][1]               # last base of the stop codon in transcript coordinates
c211 = OFF + p8 + 1 + (cds_end + 211 - ex8[0])     # GRCh38 position of c.*211
L = pd.read_csv('data/papers_tables/luo2014_table4.tsv', sep='\t', comment='#')
SHIFT = c211 - int(L[L.luo_g == 'c27706_27707'].grch37.iloc[0])
G = collections.defaultdict(list)
for r in t.itertuples():
    if r.tip not in idx: continue
    if r.copy_class == 'SMN1' and r.tip in ml.index and isinstance(ml.loc[r.tip, 'sublineage'], str): G[ml.loc[r.tip, 'sublineage']].append(idx[r.tip])
    elif r.copy_class == 'SMN2': G['SMN2'].append(idx[r.tip])
W = 15                              # window either side; alleles compared as gap-free window sequences, so an indel the
                                    # aligner shifted (e.g. delAT placed at TA) still counts as the same allele
rows = []
for v in L.itertuples():
    p38 = v.grch37 + SHIFT; o = p38 - OFF - 1
    assert rseq[o:o + len(v.ref)] == v.ref, (v.luo_g, rseq[o - 3:o + 5])
    rw = rseq[o - W:o + len(v.ref) + W]; aw = rseq[o - W:o] + ('' if v.alt == '-' else v.alt) + rseq[o + len(v.ref):o + len(v.ref) + W]
    c0, c1 = col_of[o - W], col_of[o + len(v.ref) + W - 1]
    for grp in ['A', 'B1', 'B2', 'B-', 'SMN2']:
        ws = collections.Counter(bytes(M[i][c0:c1 + 1]).decode().replace('-', '') for i in G[grp])
        top = ws.most_common(1)[0]
        rows.append(dict(luo_g=v.luo_g, grch37=v.grch37, grch38=p38, region=v.region, ref=v.ref, alt=v.alt, group=grp,
                         n=len(G[grp]), n_alt=ws[aw], n_ref=ws[rw], n_other=len(G[grp]) - ws[aw] - ws[rw],
                         top_window=top[0], top_n=top[1], ref_window=rw, alt_window=aw))
R = pd.DataFrame(rows); R.to_csv(f'{O}/luo_variants.tsv', sep='\t', index=False)
print('shift', SHIFT, 'c.*211 at', c211)
print(R.pivot(index='luo_g', columns='group', values='n_alt').to_string()); print(R.groupby('group').n.first().to_dict())
