"""Pipeline step 13b: the sites the history (step 16) and the ape-copy typing (steps 14, 15) read, moved out of step 16
so that 14 and 15 can run after it. S sites: gene-body columns where SMN1 (A+B) and SMN2 are each >=90% fixed for
different bases; AB sites: where lineages A and B (step 05) are. The ancestral base at a site is the base of at least
3 of the 4 apes (else unresolved). Reads results/sv/copies_with_outgroups.aln.fa, copy_tips.tsv, marker_lineage.tsv,
copy_labels.tsv. Writes results/sv/evolution_sites.tsv."""
import json, numpy as np, pandas as pd
O = 'results/sv'
THR = 0.9
names, seqs = [], []
for l in open(f'{O}/copies_with_outgroups.aln.fa'):
    l = l.strip()
    if l[0] == '>': names.append(l[1:].split()[0]); seqs.append([])
    else: seqs[-1].append(l.upper())
M = np.array([np.frombuffer(''.join(s).encode(), dtype=np.uint8) for s in seqs]); del seqs
idx = {n: i for i, n in enumerate(names)}
t = pd.read_csv(f'{O}/copy_tips.tsv', sep='\t')
ml = pd.read_csv(f'{O}/marker_lineage.tsv', sep='\t').set_index('tip')
lab = pd.read_csv(f'{O}/copy_labels.tsv', sep='\t').set_index('tip')
ref = idx[t[(t.copy_class == 'SMN1') & t.src.str.startswith('GCA_000001405')].tip.iloc[0]]
rpos = np.cumsum(M[ref] != ord('-')) - 1; L = rpos[-1] + 1
gene = (rpos >= 2000) & (rpos < L - 2000) & (M[ref] != ord('-'))
grp = {'A': [], 'B': [], 'SMN2': []}
for r in t.itertuples():
    if r.copy_class == 'SMN1' and r.tip in ml.index: grp['B' if ml.loc[r.tip, 'in_marker_subclade'] else 'A'].append(idx[r.tip])
    elif r.copy_class == 'SMN2': grp['SMN2'].append(idx[r.tip])
grp['SMN1'] = grp['A'] + grp['B']
APES = ['og_chimp', 'og_bonobo', 'og_gorilla', 'og_orangutan']
APE_LABEL = {'og_chimp': 'chimpanzee', 'og_bonobo': 'bonobo', 'og_gorilla': 'gorilla', 'og_orangutan': 'orangutan'}
assert all(a in idx for a in APES)
def maj(rows, c):
    col = M[rows, c]; col = col[np.isin(col, list(b'ACGT'))]
    if len(col) < 0.8 * len(rows): return None, 0.0
    v, n = np.unique(col, return_counts=True); k = n.argmax(); return chr(v[k]), n[k] / len(col)
def ancestral(c):
    b = [chr(M[idx[a], c]) for a in APES]; b = [x for x in b if x in 'ACGT']
    if not b: return None, 0
    v, n = np.unique(b, return_counts=True); k = n.argmax()
    return (str(v[k]), int(n[k])) if n[k] >= 3 else (None, int(n[k]))
rows = []
for c in np.where(gene)[0]:
    s1, f1 = maj(grp['SMN1'], c); s2, f2 = maj(grp['SMN2'], c)
    a, fa = maj(grp['A'], c); b, fb = maj(grp['B'], c)
    kinds = []
    if s1 and s2 and s1 != s2 and f1 >= THR and f2 >= THR: kinds.append('S')
    if a and b and a != b and fa >= THR and fb >= THR: kinds.append('AB')
    for k in kinds:
        anc, nag = ancestral(c)
        st1, st2 = (s1, s2) if k == 'S' else (a, b)
        r = dict(kind=k, col=int(c), grch38=int(70923086 + rpos[c] + 1), state1=st1, state2=st2,
                 A=a, B=b, SMN2=s2, SMN1=s1, ancestral=anc, apes_agreeing=nag)
        for ap in APES: r[ap] = chr(M[idx[ap], c])
        rows.append(r)
S = pd.DataFrame(rows).sort_values(['kind', 'grch38'], ascending=[False, True]).reset_index(drop=True)
S.to_csv(f'{O}/evolution_sites.tsv', sep='\t', index=False)
