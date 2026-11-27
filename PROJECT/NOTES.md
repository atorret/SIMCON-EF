# Project notes: Monte Carlo simulation of solid Ar

Running log of the project: plan, decisions, every code change and why it was made,
validations and pending work. Newest steps go at the end of section 4.

---

## 1. Assignment and scope

Assignment: [report/project_MC-solid-Ar.pdf](report/project_MC-solid-Ar.pdf)

- MC (N,V,T) code for Lennard-Jones (LJ) atoms, **N = 500**, U(r) = 4ε[(σ/r)¹² − (σ/r)⁶]
- Initial configuration: **fcc lattice** (4 superposed simple-cubic lattices)
- **ρ\* = 1.1**, T\* = 1.0, then T\* = 2, 3, 4. Compute U\*, P\*, g(r), S(k)
- Ar: m = 39.948 u, ε/k_B = 119.8 K, σ = 3.405 Å
- Questions:
  1. Metropolis algorithm and the hypotheses behind it
  2. Reduced units, and T\*, ρ\*, U\*, P\* converted to real units
  3. How Δ (trial-move size) was fixed
  4. Qualitative and quantitative changes in g(r) and S(k) when heating

**Extras chosen** (beyond the assignment):

- **Extra 1, translational order parameter** Φ = ⟨|ρ(G₁₁₁)|⟩/N. Equals 1 in a perfect crystal and
  about 0.04 in a liquid (see 5.3). Used to detect melting.
- **Extra 3, Debye–Waller check**: Bragg peaks S(G)/N ≈ exp(−G²⟨u²⟩/3). Fit ⟨u²⟩ from
  ln S(G)/N vs G².
- **Extra 5, fine temperature scan, heating and cooling**: T\* from about 0.5 to 4 in steps of about 0.25,
  heating from the fcc crystal and cooling from the melted liquid. U(T), P(T) and Φ(T) should show the
  jump at melting and a hysteresis loop (superheating vs supercooling in a small NVT box).

---

## 2. Reuse from P3 (`P3/report/code/mcnvt.f90`)

Kept unchanged because it was already validated in P3: `lj`, `eneri`, `energy`, `mcmove`
(Metropolis, single-atom moves, minimum image), `sample` (U, virial P, g(r) histogram,
energy-drift check), `stats` (block-averaged error of the mean), Δ tuning towards 50 %
acceptance during equilibration (frozen in production), tail corrections, `mc-summary.dat`.

Also reused: the workflow (one folder per run, notebook → figures → LaTeX), the report template
and logos, and text for Questions 1 and 2 from the P3 report.

---

## 3. Folder layout and workflow

```
PROJECT/
  NOTES.md                this file
  simulation/
    mcsolid.f90           MC code (started as an exact copy of P3 mcnvt.f90)
    mc-solid.data         input template
    mcsolid.exe           compiled with: gfortran -O2 -Wall -fcheck=all -o mcsolid.exe mcsolid.f90
    checks/               independent Python validation scripts (see section 5)
    runs/                 one subfolder per run (input copied in, exe run from there)
  report/
    project_MC-solid-Ar.pdf
    images/               UPC/EF logos (from P3)
    figures/              PDFs exported by the notebook
    code/                 final code copy for the appendix (filled at the end)
```

Workflow rule: plots go into the notebook (`report/plots_project.ipynb`) first, then into the LaTeX report.

### Input file `mc-solid.data`

```
ncycl nsamp nequil      MC moves in production, U/P/g(r) every nsamp moves, equilibration moves
natoms                  must be 4 n^3 (108, 256, 500, 864, ...)
tref                    T*
sigma, epsil            3.405 (A), 119.8 (K)  -- only used to write newconf.data in A
deltax                  initial Delta* (tuned during equilibration)
rho                     rho*
iseed                   random seed (reproducibility)
nsampk kmax             S(k)/Phi every nsampk moves, |k| <= kmax (1/sigma)   [new line in the project]
```

Note: one "cycle" = one trial move of one random atom; 1 sweep = N = 500 moves.

### Output files (per run folder)

| file | content |
|---|---|
| `energy-mc_TX.XX.dat` | move index, U\*/N (with tail) |
| `pressure-mc_TX.XX.dat` | move index, P\* (with tail) |
| `rdf-mc_TX.XX.dat` | r\*, g(r) (200 bins up to L/2) |
| `sk-mc_TX.XX.dat` | k\*, S(k), number of k vectors in the shell |
| `sk-lattice.dat` | the same for the perfect lattice (T\* = 0 reference) |
| `bragg-mc_TX.XX.dat` | h²+k²+l², G\*, S(G)/N, σ, error of the mean, number of vectors |
| `order-mc_TX.XX.dat` | move index, Φ |
| `mc-summary.dat` | one line per run: parameters, U, P, Φ with σ and errors |
| `newconf.data` | last configuration (Å) + box length |

---

## 4. Code log (what changed vs P3 and why)

### Step 0: setup

- Created the folder layout above and copied P3 `mcnvt.f90` to `simulation/mcsolid.f90` **unchanged first**,
  so every later change is a clean diff against validated code:
  `git diff --no-index P3/report/code/mcnvt.f90 PROJECT/simulation/mcsolid.f90`

### Step 1: perfect fcc lattice as the initial configuration

| Change | Why |
|---|---|
| New subroutine `fcc_lattice(natoms, rho, r, boxlength, ncell)`: 4 simple-cubic lattices of side a shifted by (0,0,0), (½,½,0), (½,0,½), (0,½,½)·a. 4 atoms per cube, so ρ\* = 4/a³ and a = (4/ρ\*)^(1/3); L = ncell·a | The project needs a solid starting from fcc sites (PDF hint). P3 read a liquid configuration from MD, which would never crystallise by itself. |
| Stops if N ≠ 4n³ | A crystal that does not fill the box is not periodic: it gives defects and silently wrong results. A crash is better. |
| Removed reading `leap-conf.data`, the `reduced` subroutine, the density rescaling and the "ρ ≤ 0 → use the file's density" fallback | The lattice is built directly in σ units at the requested ρ\*. No hidden default for ρ\*. |
| `r(3,1000)` → allocatable `r(3,natoms)`; subroutines use `r(3,natoms)` | The fixed size would silently overflow for N > 1000 (finite-size runs). |
| Output files labelled by T (`*_TX.XX.dat`); T\* first column in `mc-summary.dat` | The scan is in T at fixed ρ (P3 scanned ρ at fixed T). |
| Stops if T\* ≤ 0 or ρ\* ≤ 0 | Fail loudly on unphysical input. |
| Prints U\*/N of the perfect lattice before any move | Validation hook (see 5.1). |

Bug found while compiling: the program was renamed but not its `end program` line. This was fixed.

### Step 2: structure factor S(k), Bragg peaks S(G), order parameter Φ

Definitions: ρ(**k**) = Σⱼ exp(i**k**·**r**ⱼ),  S(k) = ⟨|ρ(k)|²⟩/N,
S(G) = S at fcc reciprocal-lattice vectors, Φ = ⟨|ρ(G₁₁₁)|⟩/N.

| Change | Why |
|---|---|
| New input line `nsampk kmax` | S(k) is expensive, so it is sampled less often than U and P. Both values are recorded in the input file. Stops if `ncycl` is not a multiple of `nsampk` (array sizing). |
| New subroutine `kvectors`: **k** = (2π/L)(nₓ,n_y,n_z), integer n, \|n\|² ≤ nmax², k = 0 excluded, **half space only** | Only these k are compatible with periodic boundaries. S(−k) = S(k) for real positions, so the half space halves the cost. |
| Shells labelled by the integer `kn2` = n² (no histogram bins) | All vectors with the same n² have exactly the same \|k\|, so the averaging is exact and needs no bin width. |
| Bragg vectors tagged with `km2` = h²+k²+l² when **n** = ncell·(h,k,l), h,k,l all odd or all even | Since L = 5a, the fcc reciprocal-lattice vectors G = (2π/a)(h,k,l) lie exactly on the box k grid. `fcc_lattice` now returns `ncell` for this. |
| New subroutine `sample_sk` returns the S(k) shell averages, S(G) per {hkl} family, and Φ (average over the 4 half-space (111) vectors) | One sum ρ(k) gives everything needed for Question 4 and extras 1 and 3. |
| exp(i**k**·**r**) = eₓ(nₓ)·e_y(n_y)·e_z(n_z), 1D factors computed once per atom | Only 2 complex multiplications per atom per k vector, no exponential inside the k loop. |
| Φ uses \|ρ(G)\|, not its phase | Single-atom MC moves let the whole crystal drift slowly. That changes only the phase of ρ(G), so \|ρ\| is unaffected. |
| New subroutine `write_sk` | Writes k\*, S(k), number of vectors. |
| Before the first move: S(k) of the perfect lattice → `sk-lattice.dat`, prints Φ, S(111)/N, S(200)/N | T\* = 0 reference for the plots, plus a built-in check (all must be 1). |
| Production: every `nsampk` moves accumulate S(k), store the S(G)/N series per family and Φ; Φ written to `order-mc_T*.dat` | Time series for equilibration checks and block errors. |
| End of run: `sk-mc_T*.dat`, `bragg-mc_T*.dat` (σ and block error from `stats`), Φ ± error to screen and to `mc-summary.dat` (3 new columns) | Same error analysis as U and P. |

Bugs found in testing:
- The Bragg output format had 5 reals for 4 values. Fixed.
- A runtime warning about an array temporary came from `sfam_series(family, sample)`: the slice for one family is not contiguous in memory (Fortran stores by column). Swapped to `sfam_series(sample, family)`.

**Choice k_max = 20 σ⁻¹** → nmax = 24, 28,888 half-space vectors. It covers the liquid S(k) peaks and
8 Bragg families, (111) to (420), for the Debye–Waller fit. Suggested `nsampk = 5000` (every 10 sweeps).

---

## 5. Validations done

### 5.1 fcc lattice and energy (`checks/check_fcc.py`, independent direct lattice sum, no PBC code)

| Quantity | Fortran | Independent |
|---|---|---|
| a\* | 1.53776 | 1.53776 |
| r_nn\* = a/√2, number of nearest neighbours | 1.08736 | 1.08736, 12 |
| L\* = 5a, r_c\* = L/2 | 7.6888, 3.8444 | |
| U\*/N truncated at r_c | −8.436868 | −8.436868 |
| U\*/N + tail | −8.599041 | −8.599042 |
| U\*/N infinite lattice, 2[A₁₂ r_nn⁻¹² − A₆ r_nn⁻⁶], A₁₂ = 12.13188, A₆ = 14.45392 | | −8.6081 |

**Report point:** truncated + tail differs from the exact infinite-crystal energy by **+0.009 (0.1 %)**.
The tail correction assumes g(r) = 1 beyond r_c, which is true for a liquid but not for a crystal (still discrete
shells at r_c). It should matter less at T\* = 3–4 once the system has melted.
(Note: the script's own A₆ sum gives 14.448 because r⁻⁶ converges slowly within the summed cells.
Use the literature value.)

### 5.2 S(k), Bragg peaks, Φ (`checks/check_sk.py <run_dir>`, numpy, full k sphere)

Test run: 5000 moves, `nsampk = 5000`, so S(k) is sampled exactly once, on the final configuration,
and can be compared with `newconf.data`.

| Check | Expected | Result |
|---|---|---|
| Perfect lattice: Φ, S(111)/N, S(200)/N | 1 | 1.000000 |
| Perfect lattice: S(k) away from the Bragg peaks | 0 | about 1e-28 (rounding) |
| Thermal configuration: S(k) shells vs numpy | identical | max diff 5e-6 on S up to 87 (7 significant digits, limited by the output) |
| Bragg families, Φ vs numpy | identical | agree to 7–8 digits; vector counts × 2 = full sphere |

### 5.3 Early physics hints (10 sweeps only, NOT equilibrated, do not quote as results)

- Bragg peaks already follow Debye–Waller: −ln(S/N)/G² = 0.00210 for both (111) and (420),
  so ⟨u²⟩ ≈ 0.0063 σ². Extra 3 looks promising.
- Φ in a liquid does not go to 0. With random phases, |ρ|/N has a Rayleigh distribution with mean
  √π/(2√N) ≈ **0.040** for N = 500. That is the "melted" baseline.

---

## 6. Performance (measured)

- About 4×10⁻⁵ s per MC move including U/P/g(r) every sweep (50,000 moves in 2.1 s), so **10⁷ moves ≈ 7 min**.
- One S(k) sample (k_max = 20) ≈ **26 ms**, about one sweep. Every 10 sweeps adds about 13 %.
- Pilot estimate (10⁶ equilibration + 10⁷ production + 2000 S(k) samples): **about 8–9 min**.
- The PC has **20 logical cores**. Runs are single-threaded and independent, so temperatures run in parallel
  (T\* = 1–4 in about 10 min, heating + cooling scan of about 30 runs in about 20–30 min).

---

## 7. Key numbers for the analysis and report

**Lattice at ρ\* = 1.1:** a\* = 1.53776, 2π/a = 4.08593 σ⁻¹, box k spacing dk = 2π/L = 0.81719 σ⁻¹.

Ideal fcc neighbour shells, r = r_nn√m (mark them on g(r)):

| m | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| r\* | 1.0874 | 1.5378 | 1.8834 | 2.1747 | 2.4314 | 2.6635 |
| neighbours | 12 | 6 | 24 | 12 | 24 | 8 |

Bragg peaks, G = (2π/a)√(h²+k²+l²):

| hkl | 111 | 200 | 220 | 311 | 222 | 400 | 331 | 420 |
|---|---|---|---|---|---|---|---|---|
| G\* | 7.077 | 8.172 | 11.557 | 13.551 | 14.154 | 16.344 | 17.810 | 18.273 |

**Real units (Question 2):**
- ε = 1.654×10⁻²¹ J = 0.996 kJ/mol
- Temperature unit ε/k_B = 119.8 K, so T\* = 1, 2, 3, 4 → **119.8, 239.6, 359.4, 479.2 K**
- Pressure unit ε/σ³ = **41.90 MPa**
- ρ\* = 1.1 → n = 2.786×10²⁸ m⁻³ → **1848 kg/m³ (1.85 g/cm³)**. That is denser than solid Ar at ambient
  pressure (about 1.6–1.7 g/cm³), so this is a compressed solid.
- U\*: multiply by ε (per atom) or 0.996 kJ/mol.

**Expected phase behaviour at ρ\* = 1.1** (my interpolation of the Hansen–Verlet 1969 freezing/melting data;
**verify against the papers before writing it in the report**):
- solid up to about T\* ≈ 1.9
- solid–liquid coexistence about 1.9–2.6
- liquid above about 2.6

So T\* = 1 is crystal, T\* = 2 is near or inside coexistence (the interesting case), and T\* = 3 and 4 are fluid
(consistent with the reference g(r) in the PDF). An NVT box started from fcc may **superheat**, staying
crystalline above melting, which is what extra 5's hysteresis loop shows.

**Harmonic-solid check (for U(T)):** at low T, U/N ≈ U_lattice + (3/2)T\* with U_lattice ≈ −8.60
(equipartition of the potential energy). Deviations measure anharmonicity.

---

## 8. Roadmap

| Step | Type | Content | Status |
|---|---|---|---|
| 0 | setup | folder layout, copy of P3 code | ✅ |
| 1 | code | fcc lattice + energy check | ✅ |
| 2 | code | S(k), Bragg peaks, Φ + numpy check | ✅ |
| 3 | pilot simulation | T\* = 1: equilibration length, production length (block-error plateau), Δ convergence (Q3), cost; start `plots_project.ipynb` | **next** |
| 4 | code | restart from a saved configuration (needed for the cooling branch of extra 5) | todo |
| 5 | production | T\* = 1, 2, 3, 4: U, P, g(r), S(k), Φ | todo |
| 6 | production | fine T scan, heating and cooling (extra 5) | todo |
| 7 | analysis | Debye–Waller fit (extra 3), Φ(T), melting and hysteresis, real units (Q2), S(k) direct vs FT of g(r) | todo |
| 8 | report | LaTeX, Questions 1–4 | todo |

### Pending before the pilot (Step 3)

- **Code change agreed in principle, not done yet:** nothing is recorded during equilibration (U, P and Φ are only
  sampled in production, and Δ tuning is not logged). Add an equilibration log written every time Δ is adjusted
  (every 100 sweeps): move index, Δ\*, acceptance, U\*/N (from the running `epot`, no extra cost), Φ. This is
  needed to see the relaxation away from the perfect crystal and to answer Question 3.
- Then run the pilot: T\* = 1, nequil = 10⁶, ncycl = 10⁷, nsamp = 500, nsampk = 5000, kmax = 20, Δ₀ = 0.1,
  in the background (about 8–9 min).

### Ideas and discussion points for the report

- Tail-correction error in a crystal (+0.009, section 5.1)
- Φ baseline 0.04 in a liquid, not 0 (section 5.3)
- S(k) from the direct k sum vs Fourier transform of g(r): truncation artefacts because g(r) only reaches L/2
- First-peak height and width of g(r), first minimum, coordination number (12 for the first fcc shell)
- Bragg peak heights vs T, main S(k) peak vs the Hansen–Verlet freezing value S(k_max) ≈ 2.85
- Superheating and hysteresis in a small periodic NVT box

### References to look up

- Hansen & Verlet, Phys. Rev. 184, 151 (1969): LJ phase transitions, freezing criterion
- Agrawal & Kofke, Mol. Phys. 85, 43 (1995): LJ solid–fluid coexistence
- van der Hoef, J. Chem. Phys. 113, 8142 (2000): free energy and EOS of the LJ solid
- Johnson, Zollweg & Gubbins, Mol. Phys. 78, 591 (1993): LJ fluid EOS (already used in P3; for T\* = 3, 4)
