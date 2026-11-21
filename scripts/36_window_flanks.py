"""Pipeline step 36: how each haplotype's SMN window was taken, and whether it reaches single-copy flank on both sides.
Inputs (clinical study, read-only): results/copies/windows.tsv (PANGEA query hits of a GRCh38 stretch against the HPRC
release-2 index; q_from/q_to are the hit's span on that query), results/decomp/smn.bed (the decomposition, one window
per haplotype, strand-normalised to GRCh38). A bundle is single-copy in a window when its id occurs once there.
Output: results/sv/window_flanks.tsv (one row per decomposed window) and window_flanks.json."""
import collections, json, statistics
import pandas as pd
W = pd.read_csv('results/copies/windows.tsv', sep='\t'); W = W[W.n_copies > 0]
rows = collections.defaultdict(list)
for l in open('results/decomp/smn.bed'):
    if l.startswith('#'): continue
    s, b, e, n = l.rstrip('\n').split('\t')[:4]
    rows[s].append((int(b), int(e), n.split(':')[0]))
def pick(src, end):
    c = W[W.src == src]
    if len(c) == 1: return c.iloc[0]
    c = c[(c.length - end).abs() < 50000]; assert len(c) == 1, src; return c.iloc[0]
def unique_run(v, cnt, from_start):
    """bp from the window end to the far end of the last single-copy bundle before the first bundle that recurs."""
    seq = v if from_start else v[::-1]
    last = None
    for b, e, i in seq:
        if cnt[i] > 1: break
        last = (b, e)
    if last is None: return 0
    return last[1] - v[0][0] if from_start else v[-1][1] - last[0]
out = []
for s, v in sorted(rows.items()):
    v.sort(); cnt = collections.Counter(x[2] for x in v); end = max(x[1] for x in v)
    r = pick(s, end)
    out.append(dict(src=s, ctg=r.ctg, ctg_start=int(r.ctg_start), ctg_end=int(r.ctg_end), length=int(r.length),
                    orientation=int(r.orientation), q_from=int(r.q_from), q_to=int(r.q_to),
                    first_single=cnt[v[0][2]] == 1, last_single=cnt[v[-1][2]] == 1,
                    cen_unique_bp=unique_run(v, cnt, True), tel_unique_bp=unique_run(v, cnt, False)))
D = pd.DataFrame(out); D.to_csv('results/sv/window_flanks.tsv', sep='\t', index=False)
g = D[D.src.str.startswith('GCA_000001405')].iloc[0]
qlen = int(W.q_to.max())
J = dict(n_windows=len(D), query_len=qlen,
         grch38_start1=int(g.ctg_start) + 1, grch38_end=int(g.ctg_end), grch38_q_from=int(g.q_from), grch38_q_to=int(g.q_to),
         n_rev=int((D.orientation == 1).sum()),
         n_both_single=int((D.first_single & D.last_single).sum()),
         n_reach_query_end=int((D.q_to == qlen).sum()), n_start_le_2kb=int((D.q_from <= 2000).sum()),
         cen_unique_kb_min=round(D.cen_unique_bp.min() / 1e3, 1), cen_unique_kb_med=round(statistics.median(D.cen_unique_bp) / 1e3, 1),
         tel_unique_kb_min=round(D.tel_unique_bp.min() / 1e3, 1), tel_unique_kb_med=round(statistics.median(D.tel_unique_bp) / 1e3, 1),
         len_kb_min=round(D.length.min() / 1e3), len_kb_med=round(statistics.median(D.length) / 1e3), len_kb_max=round(D.length.max() / 1e3))
assert J['grch38_q_to'] - J['grch38_q_from'] == int(g.length)
json.dump(J, open('results/sv/window_flanks.json', 'w'), indent=1, sort_keys=True); print(json.dumps(J, indent=1))
