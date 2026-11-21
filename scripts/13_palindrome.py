# Pipeline step 13 (palindrome); ported from explore/structure_vs_sequence/10_palindrome.py (sha256 741e00808eb1, commit 5e2b3bb),
# paths rewritten to results/sv/ and figures/; run from the project root after the steps before it.
"""Exploration step 10: where each SMN copy sits inside its family-A repeat unit, and how far the unit order mirrors
around an adjacent inverted A pair. Reads results/decomp/smn.repeat_units.tsv, copy_units.tsv, results/decomp/haplotype_units.tsv."""
import pandas as pd, numpy as np, collections
RU=pd.read_csv('results/decomp/smn.repeat_units.tsv',sep='\t',comment='#')
CU=pd.read_csv('results/decomp/copy_units.tsv',sep='\t')
HU=pd.read_csv('results/decomp/haplotype_units.tsv',sep='\t')
hc=pd.read_csv('tables/haplotype_copy_counts.tsv',sep='\t')
print(RU.row.nunique(), 'rows;', RU.columns.tolist()); print(hc.columns.tolist()[:12])
res=set(hc[hc.resolved].src) if 'src' in hc.columns else None
# (1) relative position of each copy midpoint inside the family-A unit containing it
rel=[]
for r in CU.itertuples():
    u=RU[(RU.row==r.src)&(RU.family=='A')&(RU.bgn<=r.mid_fwd)&(RU.end>=r.mid_fwd)]
    if len(u)!=1: continue
    u=u.iloc[0]; x=(r.mid_fwd-u.bgn)/(u.end-u.bgn)
    if u.strand=='-': x=1-x
    rel.append(dict(src=r.src,copy_id=r.copy_id,token=r.token,rel=x,unit_span=u.span_bp,dist_to_centre_kb=(r.mid_fwd-(u.bgn+u.end)/2)/1e3))
R=pd.DataFrame(rel); print('copies inside one A unit:',len(R),'of',len(CU))
print(R.rel.describe().round(3)); print(np.histogram(R.rel,bins=10,range=(0,1))[0])
print('|dist to unit centre| kb quantiles:',R.dist_to_centre_kb.abs().quantile([.1,.5,.9]).round(1).to_dict(), ' unit span kb median',round(R.unit_span.median()/1e3))
# (2) palindrome depth around the axis of an adjacent inverted A pair, by unit order
def toks(s): return [(x[:-1],x[-1]) for x in s.split()]
rows=[]
for h in HU.itertuples():
    if not isinstance(h.units,str): continue
    T=toks(h.units); best=None
    for i in range(len(T)-1):
        (f1,s1),(f2,s2)=T[i],T[i+1]
        if f1=='A' and f2=='A' and s1!=s2:
            d=0
            while i-d-1>=0 and i+2+d<len(T):
                (fl,sl),(fr,sr)=T[i-d-1],T[i+2+d]
                if fl==fr and sl!=sr: d+=1
                else: break
            side=min(i, len(T)-i-2)
            kind='head' if (s1,s2)==('+','-') else 'tail'
            # deeper mirror wins; on equal depth a head-to-head junction (A+ A-, the axis part 3 measures from) wins
            if best is None or d>best[0] or (d==best[0] and kind=='head' and best[2]=='tail'): best=(d,side,kind)
    rows.append(dict(src=h.src,arr=h.arrangement,units=h.units,n=len(T),has_pair=best is not None,depth=best[0] if best else None,room=best[1] if best else None,kind=best[2] if best else None))
P=pd.DataFrame(rows)
if res: P=P[P.src.isin(res)|P.src.str.startswith('GCA_')|P.src.str.contains('CHM13')]
print('haplotypes',len(P),' with adjacent inverted A pair',int(P.has_pair.sum()))
Q=P[P.has_pair]
print('pair kind',Q.kind.value_counts().to_dict())
print('mirror depth (unit pairs beyond the A pair) x room on the shorter side'); print(pd.crosstab(Q.depth,Q.room))
print('fully mirrored to the shorter end:',int((Q.depth==Q.room).sum()),'of',len(Q))
print(Q.groupby('arr').agg(n=('src','size'),full=('depth',lambda s: 0)).sort_values('n',ascending=False).head(0))
for a in Q.arr.value_counts().index[:6]:
    q=Q[Q.arr==a]; print(a,len(q),'median depth',q.depth.median(),'full',int((q.depth==q.room).sum()),'eg',q.units.iloc[0])
print('no pair examples:',P[~P.has_pair].arr.value_counts().head(8).to_dict())
R.to_csv('results/sv/copy_in_unit.tsv',sep='\t',index=False); P.to_csv('results/sv/palindrome_units.tsv',sep='\t',index=False)
# (3) distance of each copy from the axis of an adjacent head-to-head A pair (A+ then A-), resolved haplotypes
out = []; arm = []
for src, g in RU.groupby('row'):
    if res and src not in res: continue
    g = g.sort_values('bgn').reset_index(drop=True)
    for i in range(len(g) - 1):
        a, b = g.iloc[i], g.iloc[i + 1]
        if a.family == 'A' and b.family == 'A' and a.strand == '+' and b.strand == '-':
            axis = (a.end + b.bgn) / 2
            for c in CU[CU.src == src].itertuples():
                if a.bgn <= c.mid_fwd <= b.end: out.append(dict(src=src, token=c.token, d_kb=(c.mid_fwd - axis) / 1e3))
            arm.append(dict(src=src, left_kb=a.span_bp / 1e3, right_kb=b.span_bp / 1e3, gap_kb=(b.bgn - a.end) / 1e3))
OA = pd.DataFrame(out); AR = pd.DataFrame(arm)
OA.to_csv('results/sv/palindrome_axis_copies.tsv', sep='\t', index=False); AR.to_csv('results/sv/palindrome_arms.tsv', sep='\t', index=False)
print('A+ A- haplotypes', AR.src.nunique(), 'copies on the pair', len(OA), '|d| kb', OA.d_kb.abs().quantile([0, .5, 1]).round(1).to_dict(),
      'arms median', AR.left_kb.median(), AR.right_kb.median(), '|diff| median', (AR.left_kb - AR.right_kb).abs().median(), 'gap median', AR.gap_kb.median())
