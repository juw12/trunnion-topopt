# Rubric v1 — 5-axis trunnion upright topology optimisation (frozen before results)

Part: A-axis trunnion upright (one per side) for an AC-trunnion 5-axis FDM printer.
Domain 120 (X) x 40 (Y, axle) x 150 (Z, build) mm. 2x 6001-2RS bearings (28 mm seats, 8 mm deep,
one per face), axle centre at X 60 / Z 120. 4x M5 to the frame extrusion through two 20 mm base pads.
Material assumption: PETG-CF / PA-CF, E = 3000 MPa effective, nu = 0.35. Mass budget 105 cm3 for every design.

Compared designs (same 105 cm3, same mesh, same FE):
- R0 conventional: hand-style trapezoid web, centred
- v0 naive TO: single vertical bearing load, no overhang control (typical "first pass" TO)
- v1 improved TO: worst-case over the full bearing-load envelope (80 N any direction in the tilt
  plane + 30 N axial) + Langelaar 45 deg overhang filter + bolt/tool keep-outs

## Hard gates (pass/fail, tool-checked)
G1 Mass: binary design volume within +-5 % of 105 cm3.
G2 Connectivity: one connected solid, load path from bearing seats to base pads.
G3 Interfaces intact: bearing seats, shoulder bore, bolt holes and hex-key access columns are void;
   bosses and pads solid (checked on the exported field).
G4 Printability: unsupported down-facing elements (45 deg rule, excluding the bearing bore arc) <= 2 %.
G5 Mesh: exported STL watertight and consistently oriented; bbox within the 120 x 40 x 150 envelope.
G6 Strength: worst-case von Mises (p99) <= 10 MPa (>= 3x margin on ~30 MPa FDM inter-layer strength).
G7 Robustness: v1 worst-case envelope compliance lower than v0 and R0.

## Soft criteria (score 0-10, weights)
S1 Worst-case bearing deflection vs R0 (30 %): 4 = equal, 7 = 30 % lower, 10 = 50 % lower.
S2 Tilt uniformity, max/min seat deflection over A = -90..+90 deg (20 %): 4 = 3x, 7 = 2x, 10 = 1.3x.
S3 Printability margin, unsupported share (20 %): 4 = 5 %, 7 = 1 %, 10 = 0 %.
S4 Interface/assembly quality: real bearing seats, shoulders, bolt + key access (15 %):
   4 = loads only, 7 = seats + bolts, 10 = seats, shoulders, bolts and tool access verified.
S5 Evidence quality (15 %): 4 = picture only, 7 = FE numbers, 10 = FE + gradient checks + baselines.

Done = all gates pass and weighted score >= 8.5. Scores are self-judged unless stated otherwise;
gate results are objective (code-checked).
