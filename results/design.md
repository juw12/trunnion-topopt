# design.md - Trunnion Upright results viewer

Audience: a mechanical engineer deciding whether to print and fit this part.
One job: show which design to print, why, and what was checked - then hand over the files.
Voice: plain, technical, specific numbers with units. No hype words ("revolutionary", "seamless",
"supercharge", "unlock", "next-gen"), no exclamation marks, no emoji anywhere.

## References (what is taken from each)
- McMaster-Carr product pages: dense neutral tables, small type, data before decoration.
- Onshape / Altair Inspire result panels: model viewer as the main object, inspector beside it,
  load cases listed as a tree, one colour scale for results.
- SKF bearing catalogue pages: specification tables with units in their own column, footnotes
  for assumptions, sober two-weight typography.

## Palette (one accent, AA contrast checked on the listed backgrounds)
| token | hex | use | contrast |
|---|---|---|---|
| ink | #15181C | headings | 17.9:1 on #FFFFFF |
| text | #2A2F35 | body | 13.6:1 on #FFFFFF |
| muted | #5B6570 | secondary text | 5.9:1 on #FFFFFF, 5.4:1 on #F3F4F5 |
| line | #D7DBDF | rules, table borders | decorative |
| canvas | #F3F4F5 | page background | - |
| surface | #FFFFFF | panels | - |
| accent | #B3460E | the recommended design, primary action, part colour | 5.2:1 on #FFFFFF |
| pass | #2E6B39 | gate passed | 6.3:1 on #FFFFFF |
| fail | #A3261B | gate failed | 6.9:1 on #FFFFFF |
| baseline | #6E7C8A | comparison designs in charts | 4.5:1 on #FFFFFF |

No purple, violet or indigo anywhere. No gradients. No glow, no glass blur.

## Type
- Display / UI: "IBM Plex Sans Condensed" 600 for headings, 500 for labels.
- Body: "IBM Plex Sans" 400/500.
- Numbers and units: "IBM Plex Mono" 400, tabular figures.
- Scale (px): 12 / 13 / 14 / 16 / 20 / 28. Line height 1.45 body, 1.15 headings.

## Space, shape, motion
- 4 px grid: 4, 8, 12, 16, 24, 32, 48.
- Radius 4 px on panels and buttons, 2 px on tags. Borders 1 px #D7DBDF instead of shadows.
- Motion: none except the viewer rotation and 120 ms colour transitions on hover/focus.
- Focus ring: 2 px solid #B3460E with 2 px offset on every interactive element.

## Signature element
The tilt dial: a polar plot of bearing-seat deflection against A-axis tilt (-90 to +90 degrees),
one trace per design, the recommended design in the accent colour. It is the one picture that
says why this part is built for a 5-axis machine.

## Banned defaults
Purple/indigo gradients, Inter as the only face, emoji or emoji-like icons, sparkles/rocket/zap
icons, glow blobs, glassmorphism, gradient text, rows of identical rounded icon cards, centred
everything, fake testimonials or logos, decorative status dots, "Get started" as the only action,
stock 3D blobs, em-dash heavy copy.

## Layout (direction A, "workbench", chosen by the bracket)
- Desktop (> 1100 px): three columns, 280 px design list | viewer | 360 px inspector (spec table + checks),
  then a two-column lower band: comparison table | tilt dial, files | assumptions.
- Tablet and phone (<= 1100 px): one column with `minmax(0, 1fr)` tracks so wide tables scroll inside
  their panel instead of widening the page. The design list becomes a 2-up grid.
- Header: part name, one-line context, secondary "Download STEP" and the single primary
  "Download STL to print" in the accent colour.
- The recommended design carries a "Recommended" tag (2 px radius, accent fill) in the list and on its
  own line under its name in the comparison table.

## Decision record
- Four directions built from the same data: A workbench (IBM Plex, #B3460E), B datasheet
  (Source Serif/Sans, #1D4E89), C report (Archivo, #0E6B5C), D spec sheet (Roboto Condensed, #8A5A00).
- slop_check: 0 blockers, 0 warnings on all four.
- Council bracket (3 judges, both orders): R1 A vs B 1/1/1, A on tie-break (fewer warnings);
  final A vs D 2/0/1, A on votes. C was disqualified for a 6 px phone overflow.
- Fixes after review: A phone overflow (grid track min-content), blank viewer after a lost WebGL
  context (now restored automatically and redrawn by a ResizeObserver), "105 cm3" label that
  disagreed with the 103 cm3 measured volume.
