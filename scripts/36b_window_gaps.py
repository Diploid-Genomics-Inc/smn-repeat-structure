"""Step 36b: assembly gaps (runs of N) inside each SMN window.

Reads data/derived/windows/all_windows.fa (step 00) and tables/haplotype_copy_counts.tsv;
writes results/sv/window_gaps.tsv (one row per window: N bases, number of N runs, longest run,
window length, the clinical status and resolved flag).
"""
import re, pandas as pd
FA = 'data/derived/windows/all_windows.fa'; HC = 'tables/haplotype_copy_counts.tsv'; OUT = 'results/sv/window_gaps.tsv'
seqs = {}; name = None; buf = []
for line in open(FA):
    if line.startswith('>'):
        if name: seqs[name] = ''.join(buf)
        name = line[1:].split()[0]; buf = []
    else: buf.append(line.strip())
seqs[name] = ''.join(buf)
h = pd.read_csv(HC, sep='\t'); st = dict(zip(h.src, h.status)); rs = dict(zip(h.src, h.resolved))
rows = []
for n, s in seqs.items():
    runs = [m.end() - m.start() for m in re.finditer(r'[Nn]+', s)]
    rows.append(dict(src=n, length=len(s), n_bases=sum(runs), n_runs=len(runs), longest_run=max(runs) if runs else 0,
                     status=st.get(n, 'reference'), resolved=rs.get(n, False)))
d = pd.DataFrame(rows).sort_values('src')
# every window with N is a haplotype the clinical table calls 'gapped', and the reverse
assert set(d[d.n_bases > 0].src) == set(h[h.status == 'gapped'].src)
d.to_csv(OUT, sep='\t', index=False)
g = d[d.n_bases > 0]
print(len(d), 'windows;', len(g), 'with N; resolved among them', int(g.resolved.sum()))
print(g.sort_values('n_bases')[['src','n_bases','n_runs','longest_run']].to_string())
