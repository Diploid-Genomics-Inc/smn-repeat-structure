"""Pipeline step 41: how the data flow through the scripts to the figures and the numbers, for the supplement.
Reads the files each Python step opened (results/sv/provenance_io.tsv, traced by scripts/trace_io.py) and the
declared inputs and outputs of the four shell steps; reads note/numbers.json (written by step 40: each number, the script line and
the files it reads) and results/ledger.tsv (where each result is used).
Writes figures/provenance.png/.pdf (the data-flow diagram), results/sv/provenance_edges.tsv (step, path, read/write)
and note/tab_numsrc.tex (one row per file the numbers read: the steps that write it, how many numbers read it, the
sections that use it)."""
import re, json, glob, os, collections
import pandas as pd
import matplotlib; matplotlib.use('Agg'); matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt

# Python steps: the files each one actually opened, traced by scripts/trace_io.py (results/sv/provenance_io.tsv).
# Shell steps run external programs (PANGEA, MAFFT, IQ-TREE), which the trace cannot see: their files are declared here
# from the scripts and checked to exist.
SHELL = {'01': (['results/decomp/smn.bed'], ['results/sv/smn_struct.dist', 'results/sv/smn_struct.nwk']),
         '01b': (['results/bundles/smn.bed'], ['results/sv/fine/smn_struct.dist', 'results/sv/fine/smn_struct.nwk']),
         '02': (['results/copies/copy_calls.tsv', 'data/derived/copy_seqs'], ['results/sv/copies.fa', 'results/sv/copy_tips.tsv', 'results/sv/copies.aln.fa']),
         '03': (['data/outgroup', 'results/sv/copies.aln.fa'], ['results/sv/outgroups.fa', 'results/sv/copies_with_outgroups.aln.fa']),
         '04a': (['results/sv/copies_with_outgroups.aln.fa', 'results/sv/copy_tips.tsv'], ['results/sv/copies_body_og.fa', 'results/sv/mltree.treefile', 'results/sv/mltree.log'])}
for st_, (ins, outs) in SHELL.items():
    for f in ins + outs: assert os.path.exists(f), f
T = pd.read_csv('results/sv/provenance_io.tsv', sep='\t')
T['step'] = T.step.str.split('_').str[0]
rows = [dict(step=st_, path=f, kind=k) for st_, (ins, outs) in SHELL.items() for k, fs in (('read', ins), ('write', outs)) for f in fs]
rows += [dict(step=r.step, path=r.path, kind=r.mode) for r in T.itertuples() if r.step not in ('40', '41')]
E = pd.DataFrame(rows).drop_duplicates().sort_values(['step', 'kind', 'path'])
E = E[~E.path.str.startswith('note/')]
E.to_csv('results/sv/provenance_edges.tsv', sep='\t', index=False)
steps = sorted({p for p in glob.glob('scripts/[0-9][0-9]*_*') if p.endswith(('.py', '.sh'))} - {'scripts/40_numbers.py', 'scripts/41_provenance.py', 'scripts/42_revision.py'})   # 42 writes only the title-page revision label, from git
assert {os.path.basename(p).split('_')[0] for p in steps} == set(E.step), set(os.path.basename(p).split('_')[0] for p in steps) ^ set(E.step)
writers = E[E.kind == 'write'].groupby('path').step.apply(lambda s: sorted(set(s)))
assert writers['results/sv/marker_lineage.tsv'] == ['05'] and writers['results/sv/evolution_sites.tsv'] == ['13b'], 'each result has one writer'
assert (writers.map(len) == 1).all(), writers[writers.map(len) > 1]

# ---- numbers by source file (note/numbers.json, step 40)
R = json.load(open('note/numbers.json'))
src = collections.defaultdict(set)
for n in R['numbers']:
    for f in n['reads'] or ['(program versions and method settings)']:
        src[os.path.normpath(f)].add(n['name'])
led = pd.read_csv('results/ledger.tsv', sep='\t').fillna('')
def where(f):   # the most specific ledger row (longest pattern) that matches
    hit = [r for r in led.itertuples() if re.match('^' + re.escape(r.result).replace('\\*', '.*') + '$', f)]
    return max(hit, key=lambda r: len(r.result.replace('*', ''))).where if hit else ''
def wrote(f):
    for k, v in writers.items():
        if re.fullmatch(re.escape(k).replace('\\*', '.*'), f) or k == f: return ', '.join(v)
    if f.startswith('scripts/'): return 'setting'
    return 'input' if f.startswith(('tables/', 'results/copies', 'results/bundles', 'data/', 'papers/')) else ''
tab = sorted(src.items(), key=lambda kv: (-len(kv[1]), kv[0]))
tex = lambda s: s.replace('_', r'\_').replace('%', r'\%').replace('&', r'\&')
lines = []
for f, ms in tab:
    ex = ', '.join(r'\texttt{\textbackslash ' + m + '}' for m in sorted(ms)[:2])
    lines.append(f'\\texttt{{{tex(f.replace("results/sv/", ""))}}} & {tex(wrote(f))} & {len(ms)} & {ex} & {tex(where(f))} \\\\')
n_all = len(R['numbers']); n_multi = sum(1 for n in R['numbers'] if len(n['reads']) > 1)
open('note/tab_numsrc.tex', 'w').write('% written by scripts/41_provenance.py from note/numbers.json, results/ledger.tsv and the traced file accesses (results/sv/provenance_io.tsv)\n'
    '\\newcommand{\\tabNumSrc}{%\n' + '\n'.join(lines) + '}\n'
    f'\\newcommand{{\\pvNumbers}}{{{n_all}}}\n\\newcommand{{\\pvMulti}}{{{n_multi}}}\n\\newcommand{{\\pvFiles}}{{{len(tab)}}}\n'
    f'\\newcommand{{\\pvSteps}}{{{len(steps)}}}\n\\newcommand{{\\pvFigures}}{{{len(R["figures"])}}}\n')

# ---- the diagram, one row per step in run order: which inputs it reads, which earlier steps' files it reads,
# what it writes (figures), and how many of the document's numbers read its files
st_of = {os.path.basename(p).split('_')[0]: os.path.basename(p).rsplit('.', 1)[0] for p in steps}
order = sorted(st_of)
def group(p):
    if p.startswith(('results/bundles', 'results/copies', 'tables/')): return 0
    if p.startswith('data/'): return 1
    if p.startswith('papers/'): return 2
    return None
GN = ['clinical-study\ntables', 'data\n(NCBI, refs)', 'papers\n(Chen 2023)']
prod = {}
for s_, p_ in zip(E[E.kind == 'write'].step, E[E.kind == 'write'].path): prod.setdefault(p_, s_)
from_steps = collections.defaultdict(set)
for r in E[E.kind == 'read'].itertuples():
    for p_, s_ in prod.items():
        if s_ != r.step and s_ < r.step and (p_ == r.path or re.fullmatch(re.escape(p_).replace('\\*', '.*'), r.path)): from_steps[r.step].add(s_)
num_by_step = collections.Counter()
for n in R['numbers']:
    for s_ in {prod.get(os.path.normpath(f)) for f in n['reads']} - {None}: num_by_step[s_] += 1
figs = E[(E.kind == 'write') & E.path.str.startswith('figures/') & E.path.str.endswith('.png')].groupby('step').path.apply(lambda x: ', '.join(sorted(os.path.basename(v)[:-4] for v in x)))
fig, ax = plt.subplots(figsize=(13, 0.33 * len(order) + 1.6)); ax.axis('off')
X = dict(g=0.0, step=1.25, frm=4.0, fig=7.6, num=11.2); top = len(order)
for k, t in [('g', 'inputs read'), ('step', 'step (scripts/)'), ('frm', 'reads files written by steps'), ('fig', 'figures drawn'), ('num', 'numbers in the note\nthat read its files')]:
    ax.text(X[k], top + 1.9, t, fontsize=8.5, weight='bold', va='bottom')
for gi, g in enumerate(GN): ax.text(X['g'] + gi * 0.38, top + 0.15, g, fontsize=5.5, rotation=60, va='bottom', ha='left', color='#ca6f1e')
for i, s_ in enumerate(order):
    y = top - 1 - i
    if i % 2 == 0: ax.add_patch(plt.Rectangle((-0.1, y - 0.45), 12.6, 0.9, color='#f4f6f6', lw=0))
    gs = set(filter(lambda v: v is not None, E[(E.step == s_) & (E.kind == 'read')].path.map(group)))
    for gi in range(3):
        ax.plot(X['g'] + gi * 0.38 + 0.1, y, 'o', ms=5, color='#ca6f1e' if gi in gs else 'white', mec='#ca6f1e', mew=0.6)
    ax.text(X['step'], y, st_of[s_], fontsize=7.2, va='center', family='monospace', color='#1a5276' if s_ in figs.index else 'k')
    ax.text(X['frm'], y, ', '.join(sorted(from_steps[s_])) or '-', fontsize=7, va='center')
    ax.text(X['fig'], y, figs.get(s_, '-'), fontsize=6.6, va='center')
    n_ = num_by_step.get(s_, 0)
    ax.barh(y, n_ / max(num_by_step.values()) * 1.0, left=X['num'], height=0.5, color='#1e8449')
    ax.text(X['num'] + n_ / max(num_by_step.values()) * 1.0 + 0.05, y, str(n_), fontsize=6.6, va='center')
ax.text(0, -0.9, f'Step 40 (scripts/40_numbers.py) reads these files and writes all {n_all} numbers as macros with checks (note/numbers.tex) '
        'and their record (note/numbers.json); step 41 draws this figure; note.tex and supplement.tex are compiled to note.pdf.', fontsize=7.5)
ax.set_xlim(-0.2, 12.6); ax.set_ylim(-1.3, top + 3.0)
for ext in ('png', 'pdf'): fig.savefig(f'figures/provenance.{ext}', dpi=200, bbox_inches='tight', **({'metadata': {'CreationDate': None}} if ext == 'pdf' else {}))
print(len(steps), 'steps;', len(E), 'edges;', sum(len(v) for v in from_steps.values()), 'step-to-step;', len(tab), 'source files;', n_all, 'numbers;', n_multi, 'read several files')
