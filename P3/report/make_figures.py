"""Figures and numbers for the Practice 3 report (run from P3/report).

Reads the runs in ../simulation/runs/rhoX.X/ and writes vector PDFs to figures/
plus a small text file (figures/numbers.txt) with the values quoted in the report.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RUNS = HERE.parent / "simulation" / "runs"
FIG = HERE / "figures"
FIG.mkdir(exist_ok=True)
T_STAR = 2.0

# Johnson, Zollweg & Gubbins, Mol. Phys. 78, 591 (1993), T* = 2.0 (U per atom)
REF = pd.DataFrame({
    "rho": [0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8],
    "P_ref": [0.1777, 0.3290, 0.705, 1.069, 1.756, 3.024, 5.28],
    "U_ref": [-0.667, -1.306, -2.538, -3.149, -3.746, -4.300, -4.753],
})

plt.rcParams.update({
    "font.size": 12, "axes.labelsize": 13, "axes.titlesize": 12, "legend.fontsize": 11,
    "xtick.labelsize": 11, "ytick.labelsize": 11, "axes.grid": True, "grid.alpha": 0.3,
    "lines.linewidth": 1.6, "savefig.bbox": "tight",
})
C_SIM, C_REF, C_IDEAL = "#2a6fdb", "#d9822b", "#555555"

COLS = ["rho", "T", "N", "nequil", "ncycl", "nsamp", "deltax", "acc", "rc",
        "U", "sigma_U", "err_U", "P", "sigma_P", "err_P"]
frames = [pd.read_csv(f, sep=r"\s+", comment="#", names=COLS)
          for f in sorted(RUNS.glob("*/mc-summary.dat"))]
S = pd.concat(frames).round({"rho": 4}).drop_duplicates("rho", keep="last")
S = S.sort_values("rho").reset_index(drop=True)
rhos = list(S.rho)


def load(kind, rho):
    return np.loadtxt(RUNS / f"rho{rho:g}" / f"{kind}-mc_rho{rho:4.2f}.dat")


# ---------- Fig. 1: U and P along the production run ----------
show = [0.1, 0.4, 0.8]
fig, axes = plt.subplots(len(show), 2, figsize=(11, 7.5), sharex=True)
for i, rho in enumerate(show):
    row = S[S.rho == rho].iloc[0]
    for j, (kind, m, s, lab) in enumerate([("energy", row.U, row.sigma_U, r"$U^*/N$"),
                                           ("pressure", row.P, row.sigma_P, r"$P^*$")]):
        d = load(kind, rho)
        ax = axes[i, j]
        ax.plot(d[:, 0] / 1e6, d[:, 1], color=C_SIM, lw=0.4, alpha=0.8)
        ax.axhline(m, color="k", lw=1.3)
        ax.axhspan(m - s, m + s, color=C_REF, alpha=0.2, lw=0)
        span = 0.6 if kind == "energy" else 3.0     # same y span in every panel of a column
        ax.set_ylim(m - span / 2, m + span / 2)
        ax.set_ylabel(lab)
        ax.set_title(rf"$\rho^*={rho:g}$:  {lab} $= {m:.3f} \pm {s:.3f}$ ($\sigma$)", loc="left")
for ax in axes[-1]:
    ax.set_xlabel(r"production trial moves ($\times 10^6$)")
fig.tight_layout()
fig.savefig(FIG / "traces.pdf")
plt.close(fig)

# all densities together: one panel per magnitude
cols_all = {r: plt.cm.viridis(0.9 * i / (len(rhos) - 1)) for i, r in enumerate(rhos)}
fig, (aU, aP) = plt.subplots(1, 2, figsize=(12, 5), sharex=True)
for rho in rhos:
    for ax, kind in ((aU, "energy"), (aP, "pressure")):
        d = load(kind, rho)
        ax.plot(d[:, 0] / 1e6, d[:, 1], color=cols_all[rho], lw=0.4, alpha=0.9)
    aU.plot([], [], color=cols_all[rho], lw=2, label=rf"$\rho^*={rho:g}$")
aU.set_ylabel(r"$U^*/N$")
aP.set_ylabel(r"$P^*$")
for ax in (aU, aP):
    ax.set_xlabel(r"production trial moves ($\times 10^6$)")
fig.legend(*aU.get_legend_handles_labels(), loc="lower center", ncol=8,
           bbox_to_anchor=(0.5, -0.04), frameon=False)
fig.tight_layout(rect=(0, 0.05, 1, 1))
fig.savefig(FIG / "traces_all.pdf")
plt.close(fig)

# ---------- Fig. 2: g(r) for all densities ----------
cols = {r: plt.cm.viridis(0.9 * i / (len(rhos) - 1)) for i, r in enumerate(rhos)}
fig, ax = plt.subplots(figsize=(9, 5))
for rho in rhos:
    g = load("rdf", rho)
    ax.plot(g[:, 0], g[:, 1], color=cols[rho], label=rf"$\rho^*={rho:g}$")
ax.axhline(1.0, color=C_IDEAL, lw=0.9, ls="--")
ax.set_xlim(0, 4)
ax.set_xlabel(r"$r^* = r/\sigma$")
ax.set_ylabel(r"$g(r)$")
ax.legend(ncol=2)
fig.savefig(FIG / "rdf_all.pdf")
plt.close(fig)

fig, axes = plt.subplots(2, 4, figsize=(13, 5.5), sharex=True, sharey=True)
for ax, rho in zip(axes.flat, rhos):
    g = load("rdf", rho)
    ax.plot(g[:, 0], g[:, 1], color=C_SIM)
    ax.axhline(1.0, color=C_IDEAL, lw=0.9, ls="--")
    ax.set_title(rf"$\rho^*={rho:g}$", loc="left")
    ax.set_xlim(0, 4)
for ax in axes[-1]:
    ax.set_xlabel(r"$r^*$")
for ax in axes[:, 0]:
    ax.set_ylabel(r"$g(r)$")
fig.tight_layout()
fig.savefig(FIG / "rdf_grid.pdf")
plt.close(fig)

# ---------- Fig. 3: U*, P*, Z vs rho* ----------
rr = np.linspace(0, 0.85, 100)
Z = S.P / (S.rho * T_STAR)
# one individual figure per magnitude
fig, a1 = plt.subplots(figsize=(7, 4.5))
a1.errorbar(S.rho, S.U, yerr=S.err_U, fmt="o-", color=C_SIM, ms=6, label="MC (this work)")
a1.plot(REF.rho, REF.U_ref, "s", mfc="none", color=C_REF, ms=9, mew=1.5, label="Johnson et al.")
a1.set_xlabel(r"$\rho^*$"); a1.set_ylabel(r"$U^*/N$"); a1.legend()
fig.savefig(FIG / "eos_U.pdf")
plt.close(fig)

fig, a2 = plt.subplots(figsize=(7, 4.5))
a2.errorbar(S.rho, S.P, yerr=S.err_P, fmt="o-", color=C_SIM, ms=6, label="MC (this work)")
a2.plot(REF.rho, REF.P_ref, "s", mfc="none", color=C_REF, ms=9, mew=1.5, label="Johnson et al.")
a2.plot(rr, rr * T_STAR, "--", color=C_IDEAL, label=r"ideal gas $P^*=\rho^* T^*$")
a2.set_xlabel(r"$\rho^*$"); a2.set_ylabel(r"$P^*$"); a2.legend()
fig.savefig(FIG / "eos_P.pdf")
plt.close(fig)

fig, a3 = plt.subplots(figsize=(7, 4.5))
a3.plot(S.rho, Z, "o-", color=C_SIM, ms=6, label="MC (this work)")
a3.axhline(1.0, ls="--", color=C_IDEAL, label="ideal gas $Z=1$")
a3.axhspan(0.9, 1.1, color=C_IDEAL, alpha=0.12, lw=0, label=r"$\pm 10\,\%$")
a3.set_xlabel(r"$\rho^*$"); a3.set_ylabel(r"$Z = P^*/(\rho^* T^*)$"); a3.legend()
fig.savefig(FIG / "eos_Z.pdf")
plt.close(fig)

# zoom at low density, where the ideal gas should work
fig, ax = plt.subplots(figsize=(6.5, 4.3))
ax.errorbar(S.rho[:4], S.P[:4], yerr=S.err_P[:4], fmt="o-", color=C_SIM, ms=6, label="MC (this work)")
r2 = np.linspace(0, 0.42, 50)
ax.plot(r2, r2 * T_STAR, "--", color=C_IDEAL, label=r"ideal gas $P^*=\rho^* T^*$")
ax.set_xlabel(r"$\rho^*$"); ax.set_ylabel(r"$P^*$"); ax.set_xlim(0, 0.42); ax.set_ylim(0, 0.85)
ax.legend()
fig.savefig(FIG / "eos_zoom.pdf")
plt.close(fig)


# ---------- Ex. 7: U and P from g(r) ----------
def trapz(y, x):
    return np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x))


rows = []
for _, s in S.iterrows():
    g = load("rdf", s.rho)
    r, gr = g[:, 0], g[:, 1]
    m = (r > 0) & (r <= s.rc)
    r, gr = r[m], gr[m]
    ut = (8 * np.pi / 3) * s.rho * ((1 / 3) * s.rc**-9 - s.rc**-3)
    pt = (16 * np.pi / 3) * s.rho**2 * ((2 / 3) * s.rc**-9 - s.rc**-3)
    u = 4 * (r**-12 - r**-6)
    du = -48 * r**-13 + 24 * r**-7
    U_int = 2 * np.pi * s.rho * trapz(u * gr * r**2, r) + ut
    P_int = s.rho * T_STAR - (2 * np.pi / 3) * s.rho**2 * trapz(du * gr * r**3, r) + pt
    rows.append(dict(rho=s.rho, U=s.U, U_int=U_int, dU=100 * (U_int - s.U) / abs(s.U),
                     P=s.P, P_int=P_int, dP=100 * (P_int - s.P) / abs(s.P), ut=ut, pt=pt))
C = pd.DataFrame(rows)

fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.3))
a1.plot(C.rho, C.U, "o-", color=C_SIM, ms=6, label="direct average")
a1.plot(C.rho, C.U_int, "D", mfc="none", color=C_REF, ms=9, mew=1.5, label=r"integral of $g(r)$")
a1.set_xlabel(r"$\rho^*$"); a1.set_ylabel(r"$U^*/N$"); a1.legend()
a2.plot(C.rho, C.P, "o-", color=C_SIM, ms=6, label="direct average")
a2.plot(C.rho, C.P_int, "D", mfc="none", color=C_REF, ms=9, mew=1.5, label=r"integral of $g(r)$")
a2.set_xlabel(r"$\rho^*$"); a2.set_ylabel(r"$P^*$"); a2.legend()
fig.tight_layout()
fig.savefig(FIG / "rdf_integration.pdf")
plt.close(fig)

# ---------- numbers quoted in the text ----------
T = S.merge(REF, on="rho", how="left")
T["Z"] = Z
with open(FIG / "numbers.txt", "w") as f:
    f.write(T.to_string() + "\n\n" + C.to_string() + "\n")
    for _, s in S.iterrows():
        g = load("rdf", s.rho)
        k = np.argmax(g[:, 1])
        f.write(f"rho={s.rho}: g max {g[k,1]:.3f} at r={g[k,0]:.3f}\n")
print(open(FIG / "numbers.txt").read())
