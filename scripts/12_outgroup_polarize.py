# Pipeline step 12 (outgroup_polarize); ported from explore/structure_vs_sequence/09_outgroup_polarize.py (sha256 7fdf9cee9e1e, commit 5e2b3bb),
# paths rewritten to results/sv/ and figures/; run from the project root after the steps before it.
"""Exploration step 09: polarize the sites that separate SMN1 lineage A from B with four great-ape SMN gene bodies
(data/outgroup/*.fa, fetched from NCBI nuccore; NCBI Gene 461829, 100984482, 101148291, 100171813), added to
out/copies.aln.fa with: mafft --thread 4 --add out/outgroups.fa --keeplength out/copies.aln.fa > out/copies_with_outgroups.aln.fa
Argument: the fixation threshold within each lineage (default 0.9). Writes out/outgroup_polarized_<thr>.tsv."""
import numpy as np, pandas as pd, sys
O='results/sv'
names,seqs=[],[]
for l in open('results/sv/copies_with_outgroups.aln.fa'):
    l=l.strip()
    if l[0]=='>': names.append(l[1:].split()[0]); seqs.append([])
    else: seqs[-1].append(l.upper())
A=np.array([np.frombuffer(''.join(s).encode(),dtype=np.uint8) for s in seqs])
idx={n:i for i,n in enumerate(names)}
t=pd.read_csv(f'{O}/copy_tips.tsv',sep='\t')
ml=pd.read_csv('results/sv/marker_lineage.tsv',sep='\t').set_index('tip')
ref=idx[t[(t.copy_class=='SMN1')&t.src.str.startswith('GCA_000001405')].tip.iloc[0]]
rpos=np.cumsum(A[ref]!=ord('-'))-1; L=rpos[-1]+1
gene=((rpos>=2000)&(rpos<L-2000))
grp={'A':[],'B':[],'SMN2':[]}
for r in t.itertuples():
    if r.copy_class=='SMN1' and r.tip in ml.index: grp['B' if ml.loc[r.tip,'in_marker_subclade'] else 'A'].append(idx[r.tip])
    elif r.copy_class=='SMN2': grp['SMN2'].append(idx[r.tip])
OG=[n for n in names if n.startswith('og_')]
def maj(rows,c):
    col=A[rows,c]; col=col[np.isin(col,list(b'ACGT'))]
    if len(col)<0.8*len(rows): return None,0
    v,n=np.unique(col,return_counts=True); k=n.argmax(); return chr(v[k]), n[k]/len(col)
thr=float(sys.argv[1]) if len(sys.argv)>1 else 0.9
rows=[]
for c in np.where(gene)[0]:
    a,fa=maj(grp['A'],c); b,fb=maj(grp['B'],c)
    if a is None or b is None or a==b or fa<thr or fb<thr: continue
    s2,_=maj(grp['SMN2'],c)
    r=dict(col=int(c),grch38_smn1_pos=int(70923086+rpos[c]+1),A=a,fA=round(fa,3),B=b,fB=round(fb,3),SMN2=s2)
    for o in OG:
        x=chr(A[idx[o],c]); r[o]='A' if x==a else 'B' if x==b else ('gap' if x=='-' else 'other')
    rows.append(r)
d=pd.DataFrame(rows); print('diagnostic A/B sites (both >=%.2f fixed):'%thr,len(d))
d['SMN2_state']=np.where(d.SMN2==d.A,'A',np.where(d.SMN2==d.B,'B','other'))
for o in OG+['SMN2_state']: print(o, d[o].value_counts().to_dict())
d.to_csv('results/sv/outgroup_polarized_%s.tsv'%thr,sep='\t',index=False)
print(d.to_string())
