"""Pipeline step 00: the window FASTA and the study's bundle decomposition.
(1) The windows (moved here from step 28, unchanged): each row of results/bundles/smn.bed (the clinical study's
    decomposition, kept as the record of which window of results/copies/windows.tsv was decomposed for each haplotype,
    the one whose length matches the row's span where a haplotype has several) is cut from
    $STRATA_DATA_HPRC/hprc_r2_full.agcrs (agc-rs get), reverse-complemented where the window is (orientation 1), so
    BED coordinates index the FASTA directly. Written to data/derived/windows/<src>.fa and all_windows.fa (not in git).
(2) The decomposition (researcher, 2026-11-15: switch to the settings of the pangea-vista 0.27.0 recipe for this
    region, scratch/smn_window.pgrec, recipe.json decomp_params): pangea cohort bundle decomp -w 48 -k 55 -r 4
    --min-span 16 --min-cov 0 --min-branch-size 29 --bundle-length-cutoff 4238 --bundle-merge-distance 4238 on
    all_windows.fa. Working files in data/derived/decomp/ (not in git; the MAP-graph index is ~0.36 GB); the BED is
    copied to results/decomp/smn.bed.
(3) pangea cohort bundle repeat-units on it (default settings) -> results/decomp/smn.repeat_units.tsv.
(4) results/decomp/copy_units.tsv: each SMN copy of results/bundles/copy_units.tsv (copy_id, src, token = class and
    strand, mid_fwd = midpoint on the window) placed in the repeat unit whose span holds its midpoint (family, strand,
    span; none: empty and 0). The rule is checked first on the clinical study's own units: it reproduces
    results/bundles/copy_units.tsv exactly.
(5) results/decomp/haplotype_units.tsv: tables/haplotype_units.tsv (one row per window: sample, population,
    arrangement of the SMN copies) with its units string (family + strand, in window order) and family-A count taken
    from the new units. Checked the same way: the rule rebuilds the clinical table exactly from the old units."""

import os, subprocess, collections, shutil
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
WD = 'data/derived/windows'; os.makedirs(WD, exist_ok=True)
ARC = os.path.join(os.environ['STRATA_DATA_HPRC'], 'hprc_r2_full.agcrs')
W = pd.read_csv('results/copies/windows.tsv', sep='\t'); W = W[W.n_copies > 0]
end = collections.defaultdict(int)
for l in open('results/bundles/smn.bed'):
    if l.startswith('#'): continue
    f = l.split('\t'); end[f[0]] = max(end[f[0]], int(f[2]))
COMP = str.maketrans('ACGTN', 'TGCAN')
def pick(src):
    c = W[W.src == src]
    if len(c) == 1: return c.iloc[0]
    c = c[(c.length - end[src]).abs() < 50000]
    assert len(c) == 1, (src, len(c)); return c.iloc[0]
def seqlen(f):
    return sum(len(l.strip()) for l in open(f) if not l.startswith('>'))
def get(src):
    f = f'{WD}/{src}.fa'; r = pick(src)
    if os.path.exists(f) and seqlen(f) == int(r.length): return src, 'kept'
    out = subprocess.run(['agc-rs', 'get', ARC, f'{src}/{r.ctg}:{int(r.ctg_start)}-{int(r.ctg_end)}'], capture_output=True, text=True, check=True).stdout
    s = ''.join(l.strip() for l in out.splitlines() if not l.startswith('>')).upper()
    assert len(s) == int(r.length), (src, len(s), r.length)
    if int(r.orientation) == 1: s = s.translate(COMP)[::-1]
    open(f, 'w').write(f'>{src}\n{s}\n'); return src, 'new'
srcs = sorted(end)
with ThreadPoolExecutor(6) as ex: res = list(ex.map(get, srcs))
print(collections.Counter(x[1] for x in res), len(srcs))
assert all(end[s] <= int(pick(s).length) for s in srcs), 'a BED row runs past its window'
with open(f'{WD}/all_windows.fa', 'w') as o:
    for s in srcs: o.write(open(f'{WD}/{s}.fa').read())

DD, RD = 'data/derived/decomp', 'results/decomp'; os.makedirs(DD, exist_ok=True); os.makedirs(RD, exist_ok=True)
SET = ['-w', '48', '-k', '55', '-r', '4', '--min-span', '16', '--min-cov', '0', '--min-branch-size', '29',
       '--bundle-length-cutoff', '4238', '--bundle-merge-distance', '4238']
print(subprocess.run(['pangea', '--version'], capture_output=True, text=True).stdout.strip())
subprocess.run(['pangea', 'cohort', 'bundle', 'decomp', '-f', f'{WD}/all_windows.fa', '-o', f'{DD}/smn'] + SET, check=True)
shutil.copyfile(f'{DD}/smn.bed', f'{RD}/smn.bed')
subprocess.run(['pangea', 'cohort', 'bundle', 'repeat-units', '-b', f'{RD}/smn.bed', '-o', f'{RD}/smn'], check=True)
for ext in ['repeat_units.bed']:
    if os.path.exists(f'{RD}/smn.{ext}'): os.remove(f'{RD}/smn.{ext}')   # the BED copy of the units is not used
rows = sorted({l.split('\t')[0] for l in open(f'{RD}/smn.bed') if not l.startswith('#')})
assert rows == srcs, 'one decomposed row per window'
def place(CU, RU):
    out = []
    for r in CU.itertuples():
        u = RU[(RU.row == r.src) & (RU.bgn <= r.mid_fwd) & (RU.end >= r.mid_fwd)]
        assert len(u) <= 1, r.copy_id
        out.append(dict(copy_id=r.copy_id, src=r.src, token=r.token, mid_fwd=r.mid_fwd,
                        unit_family=u.family.iloc[0] if len(u) else '', unit_strand=u.strand.iloc[0] if len(u) else '',
                        unit_span=int(u.span_bp.iloc[0]) if len(u) else 0))
    return pd.DataFrame(out)
CU0 = pd.read_csv('results/bundles/copy_units.tsv', sep='\t', keep_default_na=False)
old = place(CU0, pd.read_csv('results/bundles/smn.repeat_units.tsv', sep='\t', comment='#'))
assert (old.astype(str).values == CU0.astype(str).values).all(), 'the placement rule reproduces the clinical copy_units.tsv'
RU = pd.read_csv(f'{RD}/smn.repeat_units.tsv', sep='\t', comment='#')
new = place(CU0, RU); new.to_csv(f'{RD}/copy_units.tsv', sep='\t', index=False)
print('units', len(RU), 'families', RU.family.nunique(), 'copies in a unit', int((new.unit_family != '').sum()), 'of', len(new))
def unit_strings(HU, RU):
    H = HU.copy(); us, na = [], []
    for src in H.src:
        g = RU[RU.row == src].sort_values(['bgn', 'end'])
        us.append(' '.join(g.family + g.strand)); na.append(int((g.family == 'A').sum()))
    H['units'] = us; H['n_unit_A'] = na; return H
HU0 = pd.read_csv('tables/haplotype_units.tsv', sep='\t', keep_default_na=False)
chk = unit_strings(HU0, pd.read_csv('results/bundles/smn.repeat_units.tsv', sep='\t', comment='#'))
assert (chk.astype(str).values == HU0.astype(str).values).all(), 'the units rule rebuilds tables/haplotype_units.tsv'
HN = unit_strings(HU0, RU); HN.to_csv(f'{RD}/haplotype_units.tsv', sep='\t', index=False)
print('haplotype rows', len(HN), 'without units', int((HN.units == '').sum()), 'changed unit strings', int((HN.units != HU0.units).sum()))
