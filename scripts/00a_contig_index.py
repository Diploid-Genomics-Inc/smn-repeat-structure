"""Pipeline step 00a: the index of assembly contigs the study uses.
For each of the 440 windows (one per resolved haplotype plus GRCh38 and CHM13; the windows step 00 cuts), the contig
it is cut from in $STRATA_DATA_HPRC/hprc_r2_full.agcrs, with the contig's accession (the third field of its PanSN
name, e.g. HG00097#1#CM094063.1 -> CM094063.1), its length in the archive, and the window's 0-based half-open span
and orientation on it. The archive was built (its meta and pipeline_manifest tables) from the HPRC release-2 AGC
file HPRC_r2_assemblies_0.6.1.agc, one archive sample per AGC sample. Writes data/hprc_r2_contigs_used.tsv."""
import os, sqlite3, collections
import pandas as pd
ARC = os.path.join(os.environ['STRATA_DATA_HPRC'], 'hprc_r2_full.agcrs')
W = pd.read_csv('results/copies/windows.tsv', sep='\t'); W = W[W.n_copies > 0]
end = collections.defaultdict(int)
for l in open('results/bundles/smn.bed'):
    if l.startswith('#'): continue
    f = l.split('\t'); end[f[0]] = max(end[f[0]], int(f[2]))
def pick(src):                      # the same rule as step 00
    c = W[W.src == src]
    if len(c) == 1: return c.iloc[0]
    c = c[(c.length - end[src]).abs() < 50000]
    assert len(c) == 1, (src, len(c)); return c.iloc[0]
con = sqlite3.connect(f'file:{ARC}?mode=ro', uri=True)
meta = dict(con.execute('select key, value from meta'))
clen = {(s, c): n for s, c, n in con.execute('select s.name, ct.name, ct.length from contig ct join sample s on s.id = ct.sample_id')}
spec = dict(con.execute('select sample_name, source_spec from pipeline_manifest'))
agc = os.path.basename(meta['pipeline_source_agc'])
rows = []
for src in sorted(end):
    r = pick(src); key = (src, r.ctg)
    assert key in clen, key
    assert int(r.ctg_end) <= clen[key] and int(r.ctg_end) - int(r.ctg_start) == int(r.length), src
    assert spec[src].split()[-1] == src and agc in spec[src], src   # each archive sample is the AGC sample of that name
    parts = r.ctg.split('#')
    rows.append(dict(window=src, agc_sample=src, sample=parts[0], haplotype=parts[1], contig=r.ctg,
                     accession=parts[2] if len(parts) == 3 else '', contig_length=clen[key],
                     window_start=int(r.ctg_start), window_end=int(r.ctg_end), window_length=int(r.length),
                     orientation='-' if int(r.orientation) == 1 else '+'))
D = pd.DataFrame(rows)
assert len(D) == 440 and D.window.is_unique
noacc = set(D[~D.accession.str.match(r'^[A-Z]{2,6}\d+\.\d+$')].window)
assert noacc == {'GCA_000001405.15_GRCh38_no_alt_analysis_set.PanSN', 'chm13v2.0_maskedY_rCRS.fa.PanSN', 'hg002v1.1.mat_MT.PanSN', 'hg002v1.1.pat.PanSN'}, noacc   # these four are named chr5, not by accession
os.makedirs('data', exist_ok=True)
with open('data/hprc_r2_contigs_used.tsv', 'w') as o:
    o.write(f'# contigs holding the 440 windows of this study; source: {agc} (HPRC release 2), via hprc_r2_full.agcrs; '
            'window_start/end 0-based half-open on the contig; orientation - = window reverse-complemented; written by scripts/00a_contig_index.py\n')
    D.to_csv(o, sep='\t', index=False)
nacc = D.accession.str.match(r'^[A-Z]{2,6}\d+\.\d+$').sum()
print('windows', len(D), 'distinct contigs', D.contig.nunique(), 'with an accession', nacc, 'orientation', D.orientation.value_counts().to_dict())
print(D[~D.accession.str.match(r'^[A-Z]{2,6}\d+\.\d+$')][['window', 'contig']].to_string())
print(D.accession.str[:2].value_counts().head(8).to_dict())
