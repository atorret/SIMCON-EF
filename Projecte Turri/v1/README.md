# v1 — archived (do not use for the report)

This folder keeps the first round of work so it can still be looked at.
It is **not current**: the system was wrong.

## Why it is archived

`charmm-gui-DMPC-NaCl-membrane/` was built with **KCl**, not NaCl:
the PSF has 15 `POT` (K⁺) + 15 `CLA` (Cl⁻) and no `SOD` (Na⁺).
That is why every `*_Na` g(r) in `Plots/Step3/G(r)/` is all zeros.

The current work uses the new CHARMM-GUI builds (with `SOD`) in
`../charmm-ORG`, `../charmm-N-64` and `../charmm-N-192`.

## What is here

| Folder | Content |
|---|---|
| `charmm-gui-DMPC-NaCl-membrane/` | old KCl box and all its NAMD run folders |
| `Plots/Step3/` | data exported from VMD/NAMDPlot by hand (303.15 K run) |
| `Plots/_output/` | figures and `Plots_report.pdf` made from that data |
| `Plots/make_plots.py` | copy of the plotting script as it was, to regenerate the old report |
| `gR/` | an extra g(r) export (`POT-NTL.dat`) |
| `__MACOSX/` | macOS leftovers from unzipping the original `.tgz` |

### Old NAMD run folders

Only three of them are really different simulations:

| Folder | What it is |
|---|---|
| `namd-ORG` | 128 DMPC, 303.15 K (source of `Plots/Step3`) |
| `namd-TEMP-283` | same box, 283 K |
| `namd-TEMP-323` | same box, 323 K |

`namd`, `namd-FINAL`, `namd-N-64` and `namd-N-192` have exactly the same
`output.log` and `step7_production.dcd` as `namd-ORG` (copies, not separate runs).

## Regenerating the old report

```bash
python v1/Plots/make_plots.py
```
