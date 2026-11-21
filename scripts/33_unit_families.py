"""Pipeline step 33: what defines each repeat-unit family, from PANGEA's own output.
Reads results/decomp/smn.repeat_units.tsv (pangea cohort bundle repeat-units, default settings; its '# cmd' line) and
results/decomp/smn.bed (the bundle decomposition it was run on). A unit occurrence is a run of at least
--min-matches bundles that recurs within one window (direct or inverted); occurrences share a family letter when the
bp-weighted Jaccard similarity of their bundle ids to the family's mean content is >= --family-jaccard.
For each family: occurrences, windows carrying it, span (median, range), bundles per occurrence, strand, relation,
and its core bundles (ids present in >= 90% of its occurrences) with the share of the family's bp they cover.
Also checks the automatic minimum unit size (median window span / 50, clamped to 5-40 kb) against the shortest unit.
Also measures overlap between units (researcher 2026-11-14: 'mention that different units may still have overlapping
with each other? (check pangea code, and explain why)'): in position (pairs of units in one window sharing a base; the
code, repeat_units.rs, gives a new copy only segments no known copy holds) and in content (bundle ids shared between
families' core sets, and bundle ids found in units of two or more families).
Writes results/sv/unit_families.tsv and results/sv/unit_families.json."""
import json, re, subprocess
import collections, numpy as np, pandas as pd
RU, BED, O = 'results/decomp/smn.repeat_units.tsv', 'results/decomp/smn.bed', 'results/sv'
cmd = open(RU).readline().strip()
assert cmd.startswith('# cmd: pangea cohort bundle repeat-units')
opts = re.findall(r'--([a-z-]+)', cmd); assert opts == [], f'non-default options: {opts}'
ru = pd.read_csv(RU, sep='\t', comment='#')
bed = pd.read_csv(BED, sep='\t', comment='#', header=None, names=['row', 'bgn', 'end', 'name'])
bed['bid'] = bed.name.str.split(':').str[0].astype(int)
bed['len'] = bed.end - bed.bgn
span = bed.groupby('row').agg(b=('bgn', 'min'), e=('end', 'max')); span = span.e - span.b
auto = int(min(max(span.median() / 50, 5000), 40000))
by_row = {r: g.sort_values('bgn') for r, g in bed.groupby('row')}
def content(r):
    g = by_row[r.row]; g = g[(g.bgn >= r.bgn - 1) & (g.end <= r.end + 1)]
    return g.groupby('bid').len.sum()
rows = []
for fam, f in ru.groupby('family'):
    cs = [content(r) for r in f.itertuples()]
    pres = pd.Series([b for c in cs for b in c.index]).value_counts() / len(cs)
    core = sorted(pres[pres >= 0.9].index)
    tot = sum(c.sum() for c in cs); core_bp = sum(c.reindex(core).fillna(0).sum() for c in cs)
    rows.append(dict(family=fam, occurrences=len(f), windows=f.row.nunique(), max_per_window=int(f.groupby('row').size().max()),
                     span_kb_median=round(f.span_bp.median() / 1e3), span_kb_min=round(f.span_bp.min() / 1e3), span_kb_max=round(f.span_bp.max() / 1e3),
                     bundles_median=int(f.n_segments.median()), forward=int((f.strand == '+').sum()), reverse=int((f.strand == '-').sum()),
                     inverted=int((f.relation == 'inverted').sum()), direct=int((f.relation == 'direct').sum()),
                     core_bundles=len(core), core_ids=' '.join(map(str, core)), core_bp_share=round(core_bp / tot, 2) if tot else np.nan))
T = pd.DataFrame(rows).sort_values(['occurrences', 'family'], ascending=[False, True])
assert T.occurrences.sum() == len(ru)
T.to_csv(f'{O}/unit_families.tsv', sep='\t', index=False)
hlp = subprocess.run(['pangea', 'cohort', 'bundle', 'repeat-units', '--help'], capture_output=True, text=True).stdout
ver = subprocess.run(['pangea', '--version'], capture_output=True, text=True).stdout.strip()
defaults = {k.replace('-', '_'): float(v) for k, v in re.findall(r'--([a-z-]+) <[A-Z_]+>\n[ \t]+[^\n]*\[default: ([0-9.]+)\]', hlp)}
assert set(defaults) == {'max_gap', 'max_indel', 'min_unit_bp', 'min_matches', 'family_jaccard'}, defaults
assert 'median row span / 50, 5-40 kb' in hlp and defaults['min_unit_bp'] == 0, 'automatic minimum unit size'
below_min = int((ru.n_matched < defaults['min_matches']).sum())   # PANGEA writes some units with fewer matched bundles than --min-matches (reported, not asserted)
print('units with fewer than min_matches matched bundles:', below_min, 'of', len(ru))
# overlap in position: units of one window sharing a base
KMER = 55                            # the decomposition's k: neighboring bundle segments can share up to k bases
pos_ov, kmer_ov = 0, 0
for r_, g in ru.groupby('row'):
    g = g.sort_values('bgn'); o = np.maximum.accumulate(g.end.values)[:-1] - g.bgn.values[1:]
    pos_ov += int((o > KMER).sum()); kmer_ov += int(((o > 0) & (o <= KMER)).sum())
print('unit pairs overlapping by more than k:', pos_ov, '; by 1-k bases:', kmer_ov)
# overlap in content: core ids shared between the families listed in Table S1 (>= 10 occurrences)
big = T[T.occurrences >= 10]; core = {r.family: set(map(int, r.core_ids.split())) for r in big.itertuples()}
fams = list(big.family); pairs = [(a, b) for i, a in enumerate(fams) for b in fams[i + 1:]]
shared = {f'{a}/{b}': len(core[a] & core[b]) for a, b in pairs if core[a] & core[b]}
within_A = [f for f in fams if f != 'A' and core[f] and core[f] <= core['A']]
fam_of_bid = {}
for r in ru.itertuples():
    for b_ in content(r).index: fam_of_bid.setdefault(int(b_), set()).add(r.family)
multi = sum(len(v) >= 2 for v in fam_of_bid.values())
overlap = dict(position_overlapping_pairs=pos_ov, kmer_edge_overlaps=kmer_ov, families_listed=len(fams), family_pairs=len(pairs), family_pairs_sharing_core=len(shared),
               shared_core=shared, core_within_A=within_A, bundle_ids_in_units=len(fam_of_bid), bundle_ids_in_2plus_families=multi)
assert pos_ov == 0, 'units overlap in position'
print(overlap)
# families A and B share their core bundle ids (Table S1): what separates them. Each occurrence read in its own orientation
# (strand '-' reversed and flipped), consecutive repeats of a token collapsed.
def otoks(r):
    g = by_row[r.row]; g = g[(g.bgn >= r.bgn - 1) & (g.end <= r.end + 1)]
    d_ = g.name.str.split(':').str[2].astype(int).tolist(); b_ = g.bid.tolist(); t = list(zip(b_, d_))
    if r.strand == '-': t = [(b, 1 - d) for b, d in t[::-1]]
    out = []
    for x in t:
        if not out or out[-1] != x: out.append(x)
    return out
def subseq(small, big):
    it = iter(big); return all(any(x == y for y in it) for x in small)
fa, fb = ru[ru.family == 'A'], ru[ru.family == 'B']
ta, tb = [otoks(r) for r in fa.itertuples()], [otoks(r) for r in fb.itertuples()]
ca = collections.Counter(tuple(t) for t in ta).most_common(1)[0][0]
shared_core = sorted(core['A'] & core['B'])
order = tuple(x for x in ca if x[0] in shared_core)
assert len(order) == len(shared_core) and all(x[1] == 0 for x in order), order
bpA = pd.concat([content(r) for r in fa.itertuples()], axis=1).T; bpB = pd.concat([content(r) for r in fb.itertuples()], axis=1).T
first2 = lambda ts, pat: sum(t[:len(pat)] == list(pat) for t in ts)
ab = dict(same_core=core['A'] == core['B'], core_ids=shared_core, core_order=' '.join(f'{b}+' for b, _ in order),
          A_n=len(fa), B_n=len(fb), A_in_order=sum(subseq(order, t) for t in ta), B_in_order=sum(subseq(order, t) for t in tb),
          bp_median_kb={str(b): [round(float(bpA[b].median()) / 1e3, 1), round(float(bpB[b].median()) / 1e3, 1)] for b in shared_core},
          span_median_kb=[round(fa.span_bp.median() / 1e3), round(fb.span_bp.median() / 1e3)],
          B_starts_9p8p=first2(tb, [(9, 0), (8, 0)]), A_starts_9p8p=first2(ta, [(9, 0), (8, 0)]),
          A_has_10m_before_core=sum((10, 1) in t[:t.index(order[0])] for t in ta if order[0] in t), B_has_10m_before_core=sum((10, 1) in t[:t.index(order[0])] for t in tb if order[0] in t))
print('A vs B:', ab)

# Why R_A and R_B are two families although they hold the same core bundles: replay PANGEA's family rule
# (copies largest span first, ties row then start; each joins the family whose mean content, bp per bundle id
# with direction ignored, is most similar by bp-weighted Jaccard if that is >= family_jaccard, else founds one).
def bwj(a, m):
    ks = set(a) | set(m); lo = sum(min(a.get(k, 0), m.get(k, 0)) for k in ks); hi = sum(max(a.get(k, 0), m.get(k, 0)) for k in ks)
    return lo / hi if hi else 0.0
JT = defaults['family_jaccard']
rep = sorted(((int(r.span_bp), r.row, int(r.bgn), r.family, content(r).to_dict()) for r in ru.itertuples()), key=lambda x: (-x[0], x[1], x[2]))
fams, labs, founders = [], [], {}
for sp, row, bgn, lab, c in rep:
    best = max(((bwj(c, {k: v / n for k, v in s_.items()}), -k) for k, (s_, n) in enumerate(fams)), default=None)
    if best is not None and best[0] + 1e-9 >= JT:
        k = -best[1]; s_, n = fams[k]
        for b, v in c.items(): s_[b] = s_.get(b, 0) + v
        fams[k][1] = n + 1
    else:
        fams.append([dict(c), 1]); k = len(fams) - 1; founders[k] = (row, bgn, sp, None if best is None else best[0], None if best is None else -best[1])
    labs.append((lab, k))
fam_lab = {}
for lab, k in labs: fam_lab.setdefault(k, collections.Counter())[lab] += 1
assert all(len(v) == 1 for v in fam_lab.values()), fam_lab  # the replay reproduces every unit's family
k_of = {v.most_common(1)[0][0]: k for k, v in fam_lab.items()}
kA, kB = k_of['A'], k_of['B']; fB = founders[kB]; assert fB[4] == kA and fB[3] < JT
mA = {b: v / fams[kA][1] for b, v in fams[kA][0].items()}; mB = {b: v / fams[kB][1] for b, v in fams[kB][0].items()}
jA = [bwj(content(r).to_dict(), mA) for r in fa.itertuples()]; jB = [bwj(content(r).to_dict(), mB) for r in fb.itertuples()]
jBA = [bwj(content(r).to_dict(), mA) for r in fb.itertuples()]; jAB = [bwj(content(r).to_dict(), mB) for r in fa.itertuples()]
ab_family = dict(replay_matches_all=len(rep), family_jaccard=JT, A_founder=list(founders[kA][:3]),
                 B_founder=[fB[0], fB[1], fB[2]], B_founder_sim_to_A=round(fB[3], 3),
                 means_sim=round(bwj(mA, mB), 3), A_own_median=round(float(np.median(jA)), 3), B_own_median=round(float(np.median(jB)), 3),
                 A_closer_own=int(sum(a > b for a, b in zip(jA, jAB))), B_closer_own=int(sum(a > b for a, b in zip(jB, jBA))),
                 B_vs_A_median=round(float(np.median(jBA)), 3), B_ge_thr_vs_A=int(sum(x >= JT for x in jBA)),
                 mean_kb={str(b): [round(mA.get(b, 0) / 1e3, 1), round(mB.get(b, 0) / 1e3, 1)] for b in (1, 15)})
print('A/B families:', ab_family)
json.dump({'ab': ab, 'ab_family': ab_family, 'overlap': overlap, 'cmd': cmd, 'pangea_version': ver, 'defaults': defaults,
           'auto_min_unit_bp': auto, 'units_below_min_matches': below_min, 'median_window_span_bp': int(span.median()), 'shortest_unit_bp': int(ru.span_bp.min()),
           'min_matched_bundles': int(ru.n_matched.min()), 'occurrences': len(ru), 'windows': int(ru.row.nunique()), 'families': len(T)},
          open(f'{O}/unit_families.json', 'w'), indent=1)
print(T.to_string(index=False)); print(open(f'{O}/unit_families.json').read())
