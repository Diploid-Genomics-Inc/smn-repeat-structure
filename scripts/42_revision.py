#!/usr/bin/env python3
"""Write note/revision.tex: the draft's revision label.

revision = number of earlier commits of note/note.pdf + 1 (this build).
hash     = short git hash of HEAD, the commit the sources are built from.
If note/ or scripts/ differ from HEAD (other than revision.tex and note.pdf),
it stops: the hash would not describe the sources.
Run after committing the sources and before building the PDF.
"""
import subprocess, sys, pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent
def git(*a):
    return subprocess.run(["git", *a], cwd=ROOT, check=True,
                          capture_output=True, text=True).stdout.strip()

dirty = [l for l in git("status", "--porcelain", "--", "note", "scripts").splitlines()
         if not l.endswith(("note/revision.tex", "note/note.pdf")) and not l.startswith("??")]
if dirty:
    sys.exit("uncommitted changes to note/ or scripts/; commit them first:\n" + "\n".join(dirty))
n_builds = len(git("log", "--format=%H", "--", "note/note.pdf").splitlines())
rev = n_builds + 1
h = git("rev-parse", "--short=7", "HEAD")
out = ROOT / "note" / "revision.tex"
out.write_text(
    "% written by scripts/42_revision.py; do not edit\n"
    f"\\newcommand{{\\draftrev}}{{{rev}}}\n"
    f"\\newcommand{{\\draftgit}}{{{h}}}\n")
print(f"revision {rev}, sources at {h} -> {out.relative_to(ROOT)}")
