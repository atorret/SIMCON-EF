---
name: simcom-lab-report
description: >-
  Scaffolds an English LaTeX lab report for SIMCOM practice sessions (P1, P2, ...).
  Copies Fortran sources, figures, and short outputs into P*/report/material/,
  emits a cover page and per-exercise Objective/Code/Results/Discussion sections
  with visible prose placeholders, and lints asset paths. Use when generating,
  writing, or scaffolding a SIMCOM practice report, lab report, P*_report.tex,
  or a practice write-up.
---

# SIMCOM Lab Report

Scaffold an Overleaf-uploadable English report for a practice directory (`P1`, `P2`, ...). Wire up structure, listings, and figures mechanically. Leave every word of judgement as a visible placeholder.

This machine has no LaTeX compiler. Do not attempt to build a PDF. `gfortran` and `python` are available. Every Fortran program `INCLUDE`s `chdir_to_code.inc`, which `CHDIR`s to the executable's directory, so binaries must be built and run **inside** the exercise directory.

## Output layout

```
P<N>/report/P<N>_report.tex
P<N>/report/material/code/     # Fortran sources only
P<N>/report/material/images/   # figures (and optional upc-logo.png)
P<N>/report/material/data/     # short console/text output
```

Copy the skeleton from [assets/report-template.tex](assets/report-template.tex). After writing, run [scripts/lint_report.py](scripts/lint_report.py).

## Workflow

Copy this checklist and track it:

```
- [ ] 1. Discover exercises
- [ ] 2. English pass on plot scripts
- [ ] 3. Rebuild stale Fortran outputs
- [ ] 4. Regenerate figures
- [ ] 5. Copy assets into material/
- [ ] 6. Write (or refuse to overwrite) the report
- [ ] 7. Lint
- [ ] 8. Git hygiene
- [ ] 9. Summarize
```

### 1. Discover exercises

Natural-sort `P<N>/ex*` subdirectories (`ex1`, `ex2`, … `ex10`). Skip `report/`.

If none exist, treat each `.f90` under `P<N>` (not under `report/`) as one exercise.

### 2. English pass on plot scripts

In each `plot_ex*.py`, translate **user-facing strings only**: titles, axis labels, legends, `print`. Minimal diff. Do not restyle figsize, dpi, grid, or fonts.

Never edit `.f90` files. If comments or `PRINT`/`WRITE` strings are not English, flag them in the run summary.

### 3. Rebuild stale Fortran outputs

A result file (`*_results.dat`, `*_results.txt`, or a file the source `OPEN`s for write) is stale if missing or older than its `.f90`.

From **inside** the exercise directory:

```bash
gfortran -O2 -o ex<N>.exe ex<N>.f90
./ex<N>.exe
```

On Windows the binary is `ex<N>.exe`; run it from that same directory so `chdir_to_code.inc` resolves data paths.

This overwrites `.exe` files that may be tracked in git. Note every overwrite in the summary.

**Do not run** a program that `READ`s from stdin (it will hang). Compile if useful, skip execution, and flag it.

### 4. Regenerate figures

Never emit tables for bulk numeric data. Always a figure.

Run each plot script with a non-interactive backend so trailing `plt.show()` is a no-op:

```bash
MPLBACKEND=Agg python plot_ex<N>.py
```

On Windows PowerShell: `$env:MPLBACKEND='Agg'; python plot_ex<N>.py`.

If a plottable numeric data file exists and there is no `plot_ex<N>.py`, **create one in the exercise directory** following existing conventions, in English, then run it:

```python
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

DATA_FILE = Path(__file__).with_name("exN_results.dat")
PLOT_FILE = Path(__file__).with_name("exN_plot.png")

cols = np.loadtxt(DATA_FILE, unpack=True)
# plot columns with English labels
fig.savefig(PLOT_FILE, dpi=150)
plt.show()
```

### 5. Copy assets

Refresh derived copies freely:

| Source | Destination |
|---|---|
| `ex<N>.f90` | `material/code/` |
| `*_plot.png` | `material/images/` |
| Short console/text output (a few lines, not a numeric matrix) | `material/data/` |

Do not copy plot scripts into `material/code/` (Code sections are Fortran only).

Leave orphaned files in `material/` whose source vanished; warn in the summary.

### 6. Write the report

If `P<N>/report/P<N>_report.tex` **already exists**, do not touch it. Report what looks stale (source newer than `material/` copy, missing exercises, missing figures) and **ask**.

Otherwise:

1. Copy [assets/report-template.tex](assets/report-template.tex) to `P<N>/report/P<N>_report.tex`.
2. Set `\practiceid` to `N` and `\practicetitle` to `Practice N`. Leave `\reportauthors` as `[Author Name]` and `\reportdate` as `[Date]`.
3. Replace the `BEGIN EXERCISES` / `END EXERCISES` region with one section per exercise (pattern below).
4. No table of contents. No introduction. Cover, then Exercise 1.

**Code:** `\lstinputlisting` of `material/code/ex<N>.f90` only.

**Results:**

- Figure: `\includegraphics{material/images/ex<N>_plot.png}` with caption `Output of exercise N`.
- Short text: `\lstinputlisting` of `material/data/...` (empty language, not Fortran).
- Both, if both exist (figure first).

**Objective / Discussion:** visible placeholders via `\ph{...} % TODO`. Do not draft prose.

Exercise section pattern:

```latex
\section{Exercise 1}
\subsection{Objective}
\ph{State the problem and the intended goal} % TODO
\subsection{Code}
\lstinputlisting[language={[free]Fortran},label={lst:ex1}]{material/code/ex1.f90}
\subsection{Results}
\begin{figure}[H]
  \centering
  \includegraphics[width=0.9\linewidth]{material/images/ex1_plot.png}
  \caption{Output of exercise 1}
  \label{fig:ex1}
\end{figure}
\subsection{Discussion \& Notes}
\ph{Comment on critical or non-trivial components} % TODO
```

For text-only results, replace the figure with:

```latex
\lstinputlisting[language={},label={lst:ex1-out}]{material/data/ex1_results.txt}
```

### 7. Lint

```bash
python .cursor/skills/simcom-lab-report/scripts/lint_report.py P<N>/report/P<N>_report.tex
```

Fix missing-asset errors before finishing. Remaining `\ph` / `% TODO` markers are expected on a fresh scaffold.

### 8. Git hygiene

Ensure the repo `.gitignore` contains these patterns (add any that are missing). Keep the PDF tracked:

```
*.aux
*.log
*.out
*.toc
*.synctex.gz
*.fls
*.fdb_latexmk
```

### 9. Summarize

Report: overwritten binaries, translated plot scripts, generated plot scripts, skipped stdin programs, orphaned `material/` files, non-English Fortran, remaining placeholders, lint result.

## Cover rules

Replicate the ETSETB title page in the template. Keep institutional names in Catalan as proper nouns: Universitat Politècnica de Catalunya, Escola Tècnica Superior d'Enginyeria de Telecomunicacions de Barcelona (ETSETB), Grau en Enginyeria Física. All other cover labels in English (Authors, Subject: SIMCOM, Date).

Logo is optional: drop `material/images/upc-logo.png` in later; the template already guards with `\IfFileExists`.
