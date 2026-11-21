"""Run one pipeline step (a .py file) and record every project file it opens, with read or write, using Python's audit
hook on 'open'. Usage: python3 scripts/trace_io.py scripts/NN_name.py. Appends rows (step, path, mode) to
results/sv/provenance_io.tsv, after removing that step's earlier rows. Used by step 41; not a numbered step itself."""
import sys, os, runpy
step = os.path.basename(sys.argv[1]).rsplit('.', 1)[0]
ROOT = os.getcwd(); seen = {}
def hook(ev, args):
    if ev != 'open' or not args or not isinstance(args[0], (str, bytes, os.PathLike)): return
    p = os.fspath(args[0]); p = p.decode() if isinstance(p, bytes) else p
    ap = os.path.abspath(p)
    if not ap.startswith(ROOT + os.sep): return
    rel = os.path.relpath(ap, ROOT)
    if rel.startswith(('scratch', '.', 'scripts')) or '__pycache__' in rel or os.path.isdir(ap): return
    mode = args[1] if len(args) > 1 and isinstance(args[1], str) else 'r'
    flags = args[2] if len(args) > 2 and isinstance(args[2], int) else 0
    w = any(c in mode for c in 'wax+') or bool(flags & (os.O_WRONLY | os.O_RDWR | os.O_CREAT))
    k = seen.get(rel)
    seen[rel] = 'write' if (w or k == 'write') else 'read'
sys.addaudithook(hook)
sys.argv = sys.argv[1:]
try:
    runpy.run_path(sys.argv[0], run_name='__main__')
finally:
    out = 'results/sv/provenance_io.tsv'
    rows = [l for l in open(out).read().splitlines()[1:] if not l.startswith(step + '\t')] if os.path.exists(out) else []
    rows += [f'{step}\t{p}\t{m}' for p, m in sorted(seen.items()) if p != out]
    with open(out, 'w') as fh: fh.write('step\tpath\tmode\n' + '\n'.join(sorted(rows)) + '\n')
