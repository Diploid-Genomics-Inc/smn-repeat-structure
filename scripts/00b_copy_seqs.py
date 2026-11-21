"""Pipeline step 00b: the SMN copy sequences, cut from each haplotype's window (step 00).

For every copy of results/copies/copy_calls.tsv the sequence is the window from win_start - 2000 to
win_end + 2000 (clamped to the window), reverse-complemented when the copy's strand is '-', so every copy reads on
the SMN1-gene strand: the gene body with 2 kb of flank on each side. One record per copy, named by its copy_id, in
data/derived/copy_seqs/<src>.fa; step 02 aligns them.

The source study wrote these files. This step states the rule and checks it: by default it rebuilds each file from
the window and the copy table and compares (sha256) with the file already there, writing nothing. Two kinds of file
are not rebuilt, and the check allows for them:
  - a haplotype with no window of its own (step 00 writes a window only where the query hit one contig): 24 of the
    haplotypes step 02 aligns, holding 51 of its 872 copies. Their copy sequences cannot be rebuilt here;
  - HG00597 mat, whose SMN2 copy is joined from two pieces (status composite), so it is not one stretch of the
    window. Step 02 leaves joined copies out of the alignment.
With --write it writes the files it can rebuild instead of checking.
"""
import argparse, hashlib, os, sys
import pandas as pd

WD, OD, CC, FLANK = 'data/derived/windows', 'data/derived/copy_seqs', 'results/copies/copy_calls.tsv', 2000
RC = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
EXPECTED = {'HG00597_mat_hprc_r2_v1.0.1': 'SMN2 copy joined from two pieces (composite)'}

def build(src, rows):
    w = ''.join(l.strip() for l in open(f'{WD}/{src}.fa') if not l.startswith('>'))
    out = []
    for r in rows.itertuples():
        a, b = max(0, int(r.win_start) - FLANK), min(len(w), int(r.win_end) + FLANK)
        s = w[a:b]
        out.append('>%s\n%s\n' % (r.copy_id, s[::-1].translate(RC) if r.strand == '-' else s))
    return ''.join(out)

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--write', action='store_true'); a = ap.parse_args()
    cc = pd.read_csv(CC, sep='\t'); os.makedirs(OD, exist_ok=True)
    same = written = nowin = 0; unexpected = []; expected = []
    for src, rows in cc.groupby('src', sort=True):
        if not os.path.exists(f'{WD}/{src}.fa'):
            nowin += 1; continue
        text = build(src, rows.sort_values('copy_id')); p = f'{OD}/{src}.fa'
        if a.write:
            open(p, 'w').write(text); written += 1; continue
        if not os.path.exists(p):
            unexpected.append((src, 'not there')); continue
        if hashlib.sha256(text.encode()).hexdigest() == hashlib.sha256(open(p, 'rb').read()).hexdigest():
            same += 1
        elif src in EXPECTED:
            expected.append((src, EXPECTED[src]))
        else:
            unexpected.append((src, 'differs'))
    if a.write:
        print(f'{written} files written to {OD}; {nowin} haplotype(s) with no window, left alone'); return
    print(f'copy sequences rebuilt from the windows: {same} the same, {len(expected)} different as expected, '
          f'{len(unexpected)} unexplained; {nowin} haplotype(s) with no window')
    for s, why in expected: print(f'  expected: {s} ({why})')
    for s, why in unexpected: print(f'  UNEXPECTED: {s} ({why})')
    assert not unexpected, 'copy sequences differ from the rule in this step'

if __name__ == '__main__':
    main()
