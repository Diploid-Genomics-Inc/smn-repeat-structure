"""Pipeline step 28f: where the NAIP pieces of the arms (step 28b, results/sv/arm_naip_pieces.tsv) start and end,
in the exons of the NAIP RefSeq transcript NM_004536.3 (data/refs_seq/naip_refseq_nm.gb, fetched from NCBI).

Checks that NM_004536.3 is the transcript behind results/copies/features.tsv: its base 1 maps (31-mers) to the
win_start of GRCh38's 5' NAIP piece in features.tsv, and every piece end in arm_naip_pieces.tsv is an exon junction
of NM_004536.3 or the transcript end used by features.tsv (7,088; the record is 7,093 bp).
Writes results/sv/naip_exons.json."""
import json, re
import pandas as pd

GB = 'data/refs_seq/naip_refseq_nm.gb'
rec = [x for x in open(GB).read().split('//\n') if 'VERSION     NM_004536.3' in x]
assert len(rec) == 1, 'NM_004536.3 not found'
rec = rec[0]
exons = [(int(a), int(b)) for a, b in re.findall(r'^     exon\s+(\d+)\.\.(\d+)', rec, re.M)]
cds = tuple(int(x) for x in re.search(r'^     CDS\s+(\d+)\.\.(\d+)', rec, re.M).groups())
doms = [(int(a), int(b), n) for a, b, n in re.findall(r'misc_feature\s+(\d+)\.\.(\d+)\s+/gene="NAIP"[^/]*/gene_synonym=[^/]*/note="propagated[^"]*Region: (BIR \d)', rec)]
seq = ''.join(re.sub(r'[^acgt]', '', l) for l in rec.split('ORIGIN')[1].splitlines()).upper()
assert len(seq) == exons[-1][1] == 7093 and len(exons) == 17 and len(doms) == 3, (len(seq), len(exons), doms)

# 1. is this the transcript features.tsv used?  base 1 -> GRCh38 window position of GRCh38's 5' NAIP piece
feat = pd.read_csv('results/copies/features.tsv', sep='\t')
g = feat[(feat.src.str.startswith('GCA_000001405')) & (feat.gene == 'NAIP') & (feat.tx_from == 0) & (feat.kind == 'gene')]
w = ''.join(l.strip() for l in open('data/derived/windows/GCA_000001405.15_GRCh38_no_alt_analysis_set.PanSN.fa') if not l.startswith('>')).upper()
k = 31
starts = [i for i in range(len(w) - k + 1) if w.startswith(seq[:k], i)]
assert int(g.win_start.iloc[0]) in starts, (starts, g.win_start.tolist())
assert feat[feat.gene == 'NAIP'].tx_to.max() == 7088

# 2. every piece end is an exon junction (or the transcript ends used by features.tsv)
junc = {0, 7088} | {b for _, b in exons[:-1]}
pc = pd.read_csv('results/sv/arm_naip_pieces.tsv', sep='\t')
ends = sorted(set(pc.tx_from) | set(pc.tx_to))
assert set(ends) <= junc, set(ends) - junc

def exon_of(base):  # 1-based transcript base -> exon number
    return next(i + 1 for i, (a, b) in enumerate(exons) if a <= base <= b)
def codon(base):
    return (base - cds[0]) // 3 + 1
def doms_in(a, b):  # domains fully / partly inside bases a..b (1-based)
    out = {}
    for d0, d1, n in doms:
        ov = max(0, min(b, d1) - max(a, d0) + 1)
        out[n] = 'whole' if ov == d1 - d0 + 1 else ('part' if ov else 'none')
    return out

pieces = {}
for (f, t), n in pc.groupby(['tx_from', 'tx_to']).size().items():
    a, b = f + 1, t
    pieces[f'{f}-{t}'] = dict(n=int(n), first_base=a, last_base=b, first_exon=exon_of(a), last_exon=exon_of(min(b, 7093)),
                              first_codon=max(1, codon(a)) if a <= cds[1] else None,
                              last_codon=codon(min(b, cds[1])) if b >= cds[0] else None, domains=doms_in(a, b))
out = dict(transcript='NM_004536.3', length=len(seq), n_exons=len(exons), exons=exons, cds=cds,
           cds_codons=(cds[1] - cds[0] + 1) // 3, domains=[dict(name=n, first=a, last=b, exons=[exon_of(a), exon_of(b)]) for a, b, n in doms],
           grch38_base1_at=int(g.win_start.iloc[0]), feature_tx_end=7088, piece_ends=[int(x) for x in ends], pieces=pieces,
           shared_3p_from_1611_and_5p=dict(bases=[1612, 1965], exons=[exon_of(1612), exon_of(1965)], n_bases=1965 - 1612 + 1))
json.dump(out, open('results/sv/naip_exons.json', 'w'), indent=1)
for kk, v in pieces.items(): print(kk, v)
print('domains', out['domains']); print('shared', out['shared_3p_from_1611_and_5p'], 'codons', out['cds_codons'])
