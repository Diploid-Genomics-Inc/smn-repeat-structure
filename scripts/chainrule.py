"""Python port of the chaining rule of PANGEA's `repeat-units` (scratch/Pangea, pangea-bin/src/bundle/repeat_units.rs
@1d4b511c; ported in explore/a_break/a10_chain_rule.py, which reproduces the long/short call on 261 of 264 pair
haplotypes). Inverted chains only. A step from match (i', j') to (i, j) adds the bp of segment i and subtracts
|gap on arm 1 - gap on arm 2| / 4 (GAP_PENALTY_DIV); a step may skip at most max_indel+1 = 13 segments; an indel step
is refused when the skipped segments carry the bundle of either end (tandem guard); two-sided skips must form one
inverted block. dp() returns, for every match, its best score and predecessor, and every legal candidate predecessor."""
import collections
MAX_GAP, MAX_INDEL, DIV = 3, 12, 4
def inverted_block(a, b):
    n, m = len(a), len(b)
    if min(n, m) < 3: return False
    L = [[0]*(m+1) for _ in range(n+1)]
    for x in range(n):
        for y in range(m):
            L[x+1][y+1] = L[x][y]+1 if a[x][2] == b[y][2] and a[x][3] != b[y][3] else max(L[x][y+1], L[x+1][y])
    k = L[n][m]; s, l = min(n, m), max(n, m)
    return k >= 3 and k*4 >= s*3 and s-k <= MAX_GAP and l-k <= MAX_INDEL
def dp(seg):
    ln = lambda i: max(seg[i][1]-seg[i][0], 0)
    pos = collections.defaultdict(list)
    for i, s in enumerate(seg): pos[s[2]].append(i)
    step = MAX_GAP+1; reach = max(step, MAX_INDEL+1)
    m = sorted((i, j) for ix in pos.values() for i in ix for j in ix if i < j and seg[i][3] != seg[j][3])
    lead = collections.defaultdict(list)
    for i, j in m: lead[i].append(j)
    best, cands = {}, {}
    for i, j in m:
        sc, pv, cl = ln(i), None, []
        for di in range(1, reach+1):
            if di > i: break
            for jj in lead.get(i-di, ()):
                dj = jj-j; lo, hi = min(di, dj), max(di, dj)
                if lo < 1 or hi > reach: continue
                if lo > step:
                    if not inverted_block(seg[i-di+1:i], seg[j+1:jj]): continue
                elif hi > step:
                    a_, b_ = (i-di, i) if di > dj else (j, jj)
                    if any(s[2] in (seg[i][2], seg[i-di][2]) for s in seg[a_+1:b_]): continue
                gi = max(seg[i][0]-seg[i-di][1], 0); gj = max(seg[jj][0]-seg[j][1], 0)
                pen = abs(gi-gj)//DIV; s = best[(i-di, jj)][0]+ln(i)-pen
                cl.append(((i-di, jj), best[(i-di, jj)][0], pen, s, gi, gj))
                if s > sc: sc, pv = s, (i-di, jj)
        best[(i, j)] = (sc, pv); cands[(i, j)] = cl
    return best, cands
def trace(best, e):
    p = []
    while e is not None: p.append(e); e = best[e][1]
    return p[::-1]
