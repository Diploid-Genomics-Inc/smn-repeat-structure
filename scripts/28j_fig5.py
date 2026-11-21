"""Pipeline step 28j: draws main Fig. 5 (figures/tanglegram_ml.png/.pdf). It reruns step 19b with DRAW_FIG5 = True, after
step 28b, so that the structure side carries strips of the cenSMN and telSMN arm types (results/sv/arm_element_haps.tsv).
Step 19b's tables are rewritten by the rerun; they are deterministic, so they come out byte-identical."""
import runpy
runpy.run_path('scripts/19b_ml_bands.py', init_globals={'DRAW_FIG5': True}, run_name='__main__')
