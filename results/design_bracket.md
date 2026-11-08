# Design bracket

| direction | blockers | warnings | findings | status |
|---|---|---|---|---|
| A-dir-A | 0 | 0 | - | in |
| B-dir-B | 0 | 0 | - | in |
| C-dir-C | 1 | 0 | OVERFLOW | out (blocker) |
| D-dir-D | 0 | 0 | - | in |

| round | match | winner | how |
|---|---|---|---|
| R1 | A-dir-A vs B-dir-B | A-dir-A | votes A/B/tie 1/1/1 (council-judged (3 judges)); tie-break on slop warnings |
| R2 | A-dir-A vs D-dir-D | A-dir-A | votes A/B/tie 2/0/1 (council-judged (3 judges)) |

**Winner:** A-dir-A
