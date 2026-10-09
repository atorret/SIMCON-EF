---
name: simcom-lab-report
description: >-
  Writes an English LaTeX lab report for SIMCOM practice sessions (P1, P2, ...)
  in the style of the finished P1 report: ETSETB cover, then per exercise a
  verbatim Question transcribed from the practice PDF, the Fortran listing,
  figures with descriptive captions, and a quantitative Discussion & Notes drafted
  from the real program output. Listings and figures point at the original
  exercise files (no copies); a lint script checks paths and style. Use when
  generating, writing, drafting, or scaffolding a SIMCOM practice report, lab
  report, P*_report.tex, or a practice write-up.
---

# SIMCOM Lab Report

Produce an English report for a practice directory (`P1`, `P2`, ...) that reads like a finished submission, not a skeleton. **The gold-standard reference is [P1/report/P1_report.tex](../../../P1/report/P1_report.tex)**. Read it before writing and match its structure, tone, and LaTeX idioms.

Environment: Windows, PowerShell. `gfortran` and `python` are available; there is **no LaTeX compiler** on PATH (the user compiles the PDF themselves). Every Fortran program `INCLUDE`s `chdir_to_code.inc`, which `CHDIR`s to the executable's directory, so build and run binaries **inside** the exercise directory.

## Output layout

```
P<N>/report/P<N>_report.tex     # the only file this skill creates under report/
upc-logo.png                    # repo root, included as ../../upc-logo.png
```

Do not create `material/` or copy any Fortran, PNG, or data file into `report/`. Every `\lstinputlisting` / `\includegraphics` points at the exercise directory (`../ex<N>/...`).

## Workflow

Copy this checklist and track it:

```
- [ ] 1. Discover exercises and read the practice statement PDF
- [ ] 2. English pass on plot scripts
- [ ] 3. Rebuild stale Fortran outputs
- [ ] 4. Regenerate figures
- [ ] 5. Extract the numbers the discussion will cite
- [ ] 6. Write (or refuse to overwrite) the report
- [ ] 7. Lint
- [ ] 8. Git hygiene
- [ ] 9. Summarize
```

### 1. Discover exercises and read the statement

Exercise units, in natural sort order (`ex2` before `ex10`):

1. `P<N>/ex*` subdirectories, if any exist.
2. Otherwise, every subdirectory of `P<N>` containing a `.f90`, skipping `report/`, `data/`, and hidden folders (e.g. P2's `leapfrog_step1`, `leapfrog_step2`, ...). A pristine `*_org` folder handed out by the teacher is not an exercise unless the statement asks for its output.
3. Otherwise, each top-level `.f90` under `P<N>`.

Then read the practice statement PDF in `P<N>/` (e.g. `SIMCON1_..._Prac1.pdf`, `Practice2_MD_Lennard-Jones.pdf`) with the Read tool. Map each exercise unit to its question in the PDF. If a reference-results PDF exists (e.g. `Practice2_results.pdf`), read it too: the discussion should compare against it.

### 2. English pass on plot scripts

In each `plot_*.py`, translate **user-facing strings only** (titles, axis labels, legends, `print`). Minimal diff; do not restyle figsize, dpi, grid, or fonts. Axis labels must carry units, e.g. `r"$t$ (s)"`.

Never edit `.f90` files. Flag non-English comments or `PRINT`/`WRITE` strings in the summary.

### 3. Rebuild stale Fortran outputs

A result file (`*_results.dat`, `*_results.txt`, `*.dat`, or anything the source `OPEN`s for write) is stale if missing or older than its `.f90`. From inside the exercise directory:

```powershell
gfortran -O2 -o ex<N>.exe ex<N>.f90; .\ex<N>.exe
```

This overwrites `.exe` files that may be tracked in git; note every overwrite.

**Do not run** a program that `READ`s from stdin (`READ(*,...)` / `READ *`); it will hang. Compile only and flag it. Its result is documented with a terminal screenshot `ex<N>_photo.png` that the user provides (see step 6).

### 4. Regenerate figures

Never put bulk numeric data in tables; always use a figure. Run each plot script with a non-interactive backend:

```powershell
$env:MPLBACKEND='Agg'; python plot_ex<N>.py
```

PNGs are written only next to the plot script. If a plottable data file exists without a plot script, create `plot_ex<N>.py` in the exercise directory (English labels, units, `dpi=150`), following existing scripts:

```python
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

DATA_FILE = Path(__file__).with_name("exN_results.dat")
PLOT_FILE = Path(__file__).with_name("exN_plot.png")

cols = np.loadtxt(DATA_FILE, unpack=True)
fig, ax = plt.subplots(figsize=(8, 5))
# plot columns with English labels and units
fig.savefig(PLOT_FILE, dpi=150)
plt.show()
```

When an exercise compares several methods or parameters, prefer one figure per physical quantity (trajectory, error, energy, ...) with stacked panels (for example, the top panel showing the parameter sweep and the bottom panel the method comparison) over a single crowded plot. P1 ex5 is the model.

### 5. Extract the numbers

The discussion is quantitative. Before writing, pull the concrete values it will cite from the actual output files with a short throwaway Python one-liner. Do not save it in the repo. Typical values:

- Converged value, iteration count, or `N` at which a tolerance is met.
- Errors at representative points (e.g. relative error at \(N = 10^2, 10^4, 10^6, 10^7\)).
- Peak positions, maxima, drift magnitudes, final values at `t_max`.
- Grid facts that justify the method (number of points, step size, parity of intervals).

Every number written in the report must come from this step, from the statement, or from a derivation shown in the text. Never invent or eyeball values.

### 6. Write the report

If `P<N>/report/P<N>_report.tex` **already exists**, do not overwrite it unless the user asked. Report what looks stale (missing exercises, figures not included, placeholders left) and **ask**.

Otherwise copy [assets/report-template.tex](assets/report-template.tex) to `P<N>/report/P<N>_report.tex`, set `\practiceid` / `\practicetitle` (`Practice N`), keep authors (`Arnau Torret \& Arnau Montull`) and `\today`. Replace the `BEGIN EXERCISES` / `END EXERCISES` region. No table of contents, no introduction: the cover is followed by Exercise 1. Number sections `Exercise 1, 2, ...` in order even if folders are named differently.

Each exercise has exactly four subsections: **Question**, **Code**, **Results**, **Discussion \& Notes**. Put `\newpage` before every `\section` except the first.

#### Question

Transcribe the exercise statement from the PDF **verbatim** (fix only obvious typos), converted to proper LaTeX:

- Formulas become `equation` environments (with `\label` if the discussion cites them).
- Sub-parts become `\paragraph{Procedure:}` + `enumerate` (`\item[A.]`, ...) or `itemize`, mirroring the PDF.
- File names in `\texttt{}` with escaped underscores; units as `$0.0005\text{ ps}$`.
- Wrap long statements containing unbreakable numbers in `{\sloppy ... \par}`.

#### Code

```latex
\lstinputlisting[style=fortranstyle,label={lst:ex3}]{../ex3/ex3.f90}
```

Add a `caption={...}` naming the algorithm when the code is non-trivial (e.g. "Fortran implementation of composite Simpson's quadrature for the Fourier cosine transform."). Skip the caption for trivial programs.

#### Results

- Simple exercise with one figure: just the figure (`width=0.75\linewidth`).
- Several figures or important parameters: open with a sentence stating the setup (physical parameters, ranges, step sizes, methods), then introduce **each** figure with a sentence citing it via `Figure~\ref{fig:...}` before the float. Use `width=0.78\linewidth` for multi-figure exercises and `0.9\linewidth` for tall stacked panels.
- Captions are **descriptive**: what is plotted, with ranges and parameters, panel by panel ("Top: ... Bottom: ..."). Never write "Output of exercise N".
- Stdin/interactive programs: include the user's terminal screenshot `../ex<N>/ex<N>_photo.png` captioned "Terminal session of exercise N: ...". If it does not exist yet, write `\ph{Screenshot ex<N>\_photo.png pending} % TODO` instead of the `\includegraphics` and ask the user to provide it.
- Short text output (a few lines): `\lstinputlisting[style=datastyle]{../ex<N>/ex<N>_results.txt}`.
- Include every meaningful PNG in the exercise folder. Skip PNGs superseded by a split set of plots (e.g. a combined `ex5_plot.png` when `ex5_pos_plot.png`, `ex5_error_plot.png`, etc. exist), and mention the skip in the summary.
- Labels: `fig:ex<N>` for a single figure, `fig:ex<N>_<quantity>` for several.

Use `\clearpage` before Discussion when an exercise has three or more figures, so floats don't drift past it.

#### Discussion \& Notes

Draft the prose yourself, grounded in step 5. Tone: first-person plural ("we"), concise, technical. Scale with the exercise:

- **Trivial exercise** (I/O, basic syntax): one short paragraph describing what the program does and where the output goes.
- **Non-trivial exercise**: an `itemize` whose items open with a bold lead, `\item \textbf{Topic:} ...`. Cover, in this order where applicable:
  1. **Method and setup**: algorithm, why its preconditions hold (e.g. even number of intervals for Simpson), sample sizes, repetitions, the estimator formula.
  2. **Code details**: non-obvious implementation choices (tolerances and *why* that tolerance, division-by-zero guards, subroutine decomposition).
  3. **Results vs expectations**: quantitative statements with numbers from step 5, compared against the exact or reference value and the expected theoretical behaviour (e.g. MC error \(\propto 1/\sqrt{N}\), Euler energy drift, Verlet symplecticity). Use nested `itemize` with `\emph{Method:}` leads when comparing methods or sweeping a parameter.
  4. **Trade-off or conclusion**: the recommended parameter or method and why (accuracy vs cost).
- Show short derivations inline with `equation` when they back a number (e.g. \(E_0 = \tfrac12 k x_0^2\)).
- If a claim cannot be verified from the data, the statement, or a derivation, leave `\ph{...} % TODO` rather than guessing.

#### LaTeX conventions

- Units: `$0.02\text{ s}$`, `$\text{ ps}^{-1}$`. Absolute values: `\vert{}x\vert{}`. Subscripts that are words: `E_{\text{tot}}`.
- Cite with `Figure~\ref{...}`. Don't stack `\newpage\newpage`.
- Escape `_`, `%`, `&`, `#` in text.

### 7. Lint

```powershell
python .cursor/skills/simcom-lab-report/scripts/lint_report.py P<N>/report/P<N>_report.tex
```

Fix every missing asset (exit code 1). Then address the style warnings: generic captions, figures never cited with `\ref` (acceptable only for single-figure trivial exercises), and plots in exercise folders that the report does not include (justify any that remain). Remaining `\ph` / `% TODO` markers must each correspond to something genuinely unknown.

### 8. Git hygiene

Ensure `.gitignore` contains these patterns (add missing ones). Keep the PDF tracked:

```
*.aux
*.log
*.out
*.toc
*.synctex.gz
*.fls
*.fdb_latexmk
```

Bulk outputs larger than about 10 MB (like `P1/ex2/ex2_results.dat`) get their own `.gitignore` entry with a comment saying how to regenerate them.

### 9. Summarize

Report: overwritten binaries, translated or created plot scripts, skipped stdin programs and pending screenshots, PNGs deliberately left out, non-English Fortran, remaining placeholders (and what info each needs), and the lint result. Remind the user to compile the PDF.

## Cover rules

The template already contains the ETSETB title page; don't modify it. Institutional names stay in Catalan as proper nouns (Universitat Politècnica de Catalunya, Escola Tècnica Superior d'Enginyeria de Telecomunicacions de Barcelona (ETSETB), Grau en Enginyeria Física); all other labels are in English (Authors, Subject: SIMCOM, Date). The logo is `../../upc-logo.png` at the repo root, guarded by `\IfFileExists`.
