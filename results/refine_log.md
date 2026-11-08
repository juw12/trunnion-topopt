# Refine log - 5-axis trunnion upright topology optimisation

Rubric: rubric_v1.md (frozen before results). Gates are code-checked; every number below is from
linear FE on the exported geometry at 1.25 mm elements (about 1.16 M DOF), not from the optimiser.

| round | change | vol cm3 | worst C N.mm | defl vert / horiz / axial um | belt tilt arcmin | peak vM env / belt MPa | decision |
|---|---|---|---|---|---|---|---|
| R0 | conventional 5 mm web + edge flanges (hand design) | 99.7 | 6.28 | 13.7 / 66.4 / 209.4 | 70.2 | 2.23 / 11.84 | baseline; fails stress (11.8 MPa > 10) and overhang (3.6 %) |
| v0 | single-load SIMP, vertical 80 N, no overhang control | 103.5 | 6.83 | 5.5 / 50.5 / 227.6 | 28.3 | 2.51 / 6.36 | baseline; 3.0 % unsupported overhang |
| v1 | eigen p-norm over 80 N any-tilt + 30 N axial, 45 deg AM filter, keep-outs | 103.2 | 2.59 | 7.7 / 31.9 / 86.5 | 12.3 | 1.46 / 3.50 | kept |
| v2 | v1 + overhung belt case (100 N, 15 mm outboard) in the envelope | 103.2 | 2.70 | 8.3 / 33.8 / 69.0 | 8.9 | 0.97 / 3.31 | kept; belt tilt 8.9' > 5' |
| v3 | belt case weighted 1.7x, 130 cm3 budget, warm start from v2 | 128.5 | 2.20 | 8.2 / 27.5 / 57.8 | 5.6 | 0.82 / 2.40 | kept; belt tilt 5.6' > 5' |
| v4 | belt case weighted 2.2x, 145 cm3 budget, warm start from v3 | 143.5 | 1.96 | 8.2 / 24.6 / 52.9 | 4.6 | 0.70 / 2.04 | RECOMMENDED; all 6 gates pass |

Mesh convergence: v2 2.5 mm -> 1.25 mm within about 3 %; v3 belt tilt 5.43 -> 5.61 arcmin (+3.3 %).
Overhang audit (45 deg, build +Z): v1, v2, v3 0.0 % unsupported, one connected body each.
v4: 2.5 mm -> 1.25 mm belt tilt 4.42 -> 4.56 arcmin (+3.2 %); overhang 0.0 %, one body; STL watertight, all 8 interface probes pass; STEP valid single solid, 141.7 cm3, bbox 120.1 x 40.1 x 142.6 mm.

## Bugs found and fixed
- Evaluation bug (v1 report): bearing loads were spread over seat nodes that touched only void elements,
  adding a constant offset that hid the differences between designs. Loads now go on solid-contact nodes;
  all tables above use the corrected method. Corrected v1 vs v0: worst-case compliance 6.83 -> 2.59 N.mm,
  axial deflection 228 -> 87 um, belt tilt 28 -> 12 arcmin.
- STL export filled the bearing seats (band intersection instead of union). Fixed; interface probes now
  check seats, bore, bolt holes and key access on every export.
- OCC STL reader failed on sliver triangles; fixed by clamping |field| >= 0.02 voxel before meshing.
- v3 CAD step crashed silently inside the pipeline (no output, no STEP) while the parent still held the
  1.16 M DOF FE model; a standalone rerun produced a valid single solid. Cause: memory pressure.

## Results viewer (/design)
- Four directions (A workbench, B datasheet, C report, D spec sheet), all from the same data.
- slop_check: 0 blockers, 0 warnings on all four.
- Bracket (3 council judges, both orders): A vs B 1/1/1 (A on tie-break), A vs D 2/0/1 (A on votes).
  C disqualified: 390 px viewport scrolled to 396 px.
- Screenshots checked by eye at 390 and 1440 px. Fixed: A phone overflow (grid min-content), blank
  viewer after WebGL context loss, wrong '105 cm3' label.

## Stop reason
All six gates met by v4 (179 g). Mass was the price: v2 103 cm3 -> v4 144 cm3 to bring belt tilt from 8.9 to 4.6 arcmin. v2 stays available as the lighter option if the belt pull is kept below 56 N.
