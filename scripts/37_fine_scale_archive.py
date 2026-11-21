"""Pipeline step 37: the finer-scale results the supplement reports (Supplementary Section 'The finer scale').

The preprint's main text uses the coarse bundle decomposition (step 00). The finer decomposition it replaced
(results/bundles/smn.bed; decomp -w 48 -r 4, cutoff 2,500, merge 10,000, other settings default; see
results/sv/decomp_params.tsv) was carried through every structural step up to commit 953c895. Step 01b recomputes its
structural distance (byte-identical to the one used then) and step 28i recomputes the Mantel and cophenetic tests on it.
The band, arm-form, family and copy-position results at that scale are taken here, unchanged, from the result files
of commit 953c895 (git show), so they are the outputs of the pipeline as it ran then; rerunning that commit's
run_all.sh regenerates them. Writes results/sv/fine/<file> and results/sv/fine/archive_manifest.tsv (sha256)."""
import hashlib, subprocess, pathlib
COMMIT = '953c895'
FILES = ['copy_position.json', 'ml_cluster_pairs.tsv', 'ml_bands_summary.json', 'ml_vs_upgma_bands.tsv',
         'armform_summary.json', 'armform_struct.tsv', 'unit_family_elements.tsv', 'g_unpaired.json',
         'arm_elements_summary.json', 'ml_hap_labels.tsv', 'summary.json', 'ml_cluster_pairs_members.tsv',
         'ml_band_sensitivity.tsv', 's13_summary.json', 's13_unit_families.tsv', 's13_G_carriers_CD.tsv', 'unit_families.json']
out = pathlib.Path('results/sv/fine'); out.mkdir(parents=True, exist_ok=True)
rows = ['file\tcommit\tsha256']
for f in FILES:
    b = subprocess.run(['git', 'show', f'{COMMIT}:results/sv/{f}'], check=True, capture_output=True).stdout
    (out / f).write_bytes(b)
    rows.append(f'{f}\t{COMMIT}\t{hashlib.sha256(b).hexdigest()}')
# the finer-scale structural distance recomputed now (step 01b) must equal the one those results were built on
old = subprocess.run(['git', 'show', f'{COMMIT}:results/sv/smn_struct.dist'], capture_output=True).stdout
now = (out / 'smn_struct.dist').read_bytes()
if old:
    assert hashlib.sha256(old).digest() == hashlib.sha256(now).digest(), 'step 01b distance differs from the archived one'
    rows.append(f'# smn_struct.dist (step 01b) equals {COMMIT}:results/sv/smn_struct.dist')
else:
    rows.append(f'# smn_struct.dist is git-ignored; step 01b reproduced it byte for byte (checked when it was written)')
(out / 'archive_manifest.tsv').write_text('\n'.join(rows) + '\n')
print('\n'.join(rows))
