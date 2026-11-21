"""Pipeline step 28: PANGEA's tandem-repeat finder over every window (the window FASTA is written by step 00).
pangea cohort bundle tandem-repeats -f data/derived/windows/all_windows.fa --min-unit-bp 1000 --max-unit-bp 30000
-> results/sv/smn.tandem_repeats.tsv / .bed (one line per array per window: start, end, unit size, copies)."""
import subprocess
import pandas as pd
WD = 'data/derived/windows'
srcs = sorted({l.split('\t')[0] for l in open('results/decomp/smn.bed') if not l.startswith('#')})
subprocess.run(['pangea', 'cohort', 'bundle', 'tandem-repeats', '-f', f'{WD}/all_windows.fa', '-o', 'results/sv/smn',
                '--min-unit-bp', '1000', '--max-unit-bp', '30000'], check=True)
T = pd.read_csv('results/sv/smn.tandem_repeats.tsv', sep='\t', comment='#')
assert set(T.row) <= set(srcs)
print('arrays', len(T), 'windows with an array', T.row.nunique(), 'of', len(srcs))
