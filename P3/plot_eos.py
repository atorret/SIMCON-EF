"""Graficas de la practica 3: U*/N y P* frente a rho*, y g(r) por densidad.

Lee runs/summary_n*.csv y runs/n*/rho*/gdr-mc.dat (los crea run_series.py) y
guarda los PNG en figures/.   Uso:  python plot_eos.py
"""
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
RUNS, FIG = HERE / "runs", HERE / "figures"
FIG.mkdir(exist_ok=True)


def load(path):
    d = np.genfromtxt(path, delimiter=",", names=True)
    return {k: np.atleast_1d(d[k]) for k in d.dtype.names}


series = {}
for p in sorted(RUNS.glob("summary_n*.csv")):
    series[int(re.search(r"n(\d+)", p.name).group(1))] = load(p)
if not series:
    raise SystemExit("No hay runs/summary_n*.csv: ejecuta run_series.py primero")
tstar = next(iter(series.values()))["T"][0]

for key, ylab, name in (("U_per_atom", r"$U^*/N$", "U_vs_rho"), ("P", r"$P^*$", "P_vs_rho")):
    fig, ax = plt.subplots(figsize=(6, 4))
    for n, d in series.items():
        sig = d["sigma_U" if key == "U_per_atom" else "sigma_P"]
        ax.errorbar(d["rho"], d[key], yerr=sig, marker="o", capsize=3, label=f"$n_{{cycl}}$={n:.0e}")
    if key == "P":
        rho = np.linspace(0, max(d["rho"].max() for d in series.values()), 50)
        ax.plot(rho, rho * tstar, "k--", label=r"ideal gas $\rho^*T^*$")
    else:
        ax.axhline(0, color="k", ls="--", label="ideal gas")
    ax.set_xlabel(r"$\rho^*$")
    ax.set_ylabel(ylab)
    ax.set_title(rf"$T^*={tstar:g}$ (with tail corrections)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG / f"{name}.png", dpi=200)
    plt.close(fig)

n = max(series)
fig, ax = plt.subplots(figsize=(6, 4))
dirs = sorted((RUNS / f"n{n}").glob("rho*"))
colors = plt.cm.turbo(np.linspace(0.05, 0.95, len(dirs)))
for d, col in zip(dirs, colors):
    r, g = np.loadtxt(d / "gdr-mc.dat", unpack=True)
    ax.plot(r, g, color=col, label=rf"$\rho^*={float(d.name[3:]):.1f}$")
ax.set_xlabel(r"$r^*$")
ax.set_ylabel(r"$g(r)$")
ax.set_title(rf"$T^*={tstar:g}$, $n_{{cycl}}$={n:.0e}")
ax.legend(fontsize=7, ncol=2)
fig.tight_layout()
fig.savefig(FIG / "gdr.png", dpi=200)
plt.close(fig)

# Evolucion de U y P en una sola ejecucion (Exercise 1)
d = RUNS / f"n{n}" / "rho0.50"
fig, axs = plt.subplots(1, 2, figsize=(10, 3.6))
for ax, fname, ylab in ((axs[0], "energy-mc.dat", r"$U^*/N$"), (axs[1], "pressure-mc.dat", r"$P^*$")):
    c, v = np.loadtxt(d / fname, unpack=True)
    ax.plot(c / 500, v, lw=0.5, color="0.6", label="sampled value")
    ax.plot(c / 500, np.cumsum(v) / np.arange(1, len(v) + 1), color="C3", label="running mean")
    ax.set_xscale("log")
    ax.set_xlabel("MC moves per atom")
    ax.set_ylabel(ylab)
    ax.legend(fontsize=8)
fig.suptitle(rf"$\rho^*=0.5$, $T^*={tstar:g}$, $n_{{cycl}}$={n:.0e}")
fig.tight_layout()
fig.savefig(FIG / "traces_rho05.png", dpi=200)
plt.close(fig)

# U*/N frente a icycl para todas las densidades (serie de 1e6)
n1 = min(series)
fig, ax = plt.subplots(figsize=(7, 4.2))
dirs = sorted((RUNS / f"n{n1}").glob("rho*"))
for d, col in zip(dirs, plt.cm.turbo(np.linspace(0.05, 0.95, len(dirs)))):
    c, u = np.loadtxt(d / "energy-mc.dat", unpack=True)
    ax.plot(c, u, lw=0.6, color=col, label=rf"$\rho^*={float(d.name[3:]):.1f}$")
ax.set_xlim(0, n1)
ax.set_xlabel(r"icycl (MC trial moves)")
ax.set_ylabel(r"$U^*/N$")
ax.set_title(rf"$T^*={tstar:g}$, $n_{{cycl}}$={n1:.0e}")
ax.legend(fontsize=7, ncol=4, loc="lower right")
fig.tight_layout()
fig.savefig(FIG / "U_vs_icycl.png", dpi=200)
plt.close(fig)

# P* frente a icycl para todas las densidades (serie de 1e6)
fig, ax = plt.subplots(figsize=(7, 4.2))
for d, col in zip(dirs, plt.cm.turbo(np.linspace(0.05, 0.95, len(dirs)))):
    c, p = np.loadtxt(d / "pressure-mc.dat", unpack=True)
    ax.plot(c, p, lw=0.6, color=col, label=rf"$\rho^*={float(d.name[3:]):.1f}$")
ax.set_xlim(0, n1)
ax.set_xlabel(r"icycl (MC trial moves)")
ax.set_ylabel(r"$P^*$")
ax.set_title(rf"$T^*={tstar:g}$, $n_{{cycl}}$={n1:.0e}")
ax.legend(fontsize=7, ncol=4, loc="lower right")
fig.tight_layout()
fig.savefig(FIG / "P_vs_icycl.png", dpi=200)
plt.close(fig)

# U*/N frente a icycl a rho*=0.8, con la media directa y la de integrar g(r)
c, u = np.loadtxt(RUNS / f"n{n1}" / "rho0.80" / "energy-mc.dat", unpack=True)
s = series[n1]
k = int(np.argmin(abs(s["rho"] - 0.8)))
uint = load(RUNS / f"integral_n{n1}.csv")
ki = int(np.argmin(abs(uint["rho"] - 0.8)))
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(c, u, lw=0.6, color="0.55", label="sampled $U^*/N$")
ax.axhline(s["U_per_atom"][k], color="C3", lw=1.5, label=rf"mean value $\langle U^*\rangle/N={s['U_per_atom'][k]:.4f}$")
ax.axhline(uint["U_integral"][ki], color="C0", lw=1.5, ls="--", label=rf"integral of $g(r)$: $U^*/N={uint['U_integral'][ki]:.4f}$")
ax.set_xlim(0, n1)
ax.set_xlabel(r"icycl (MC trial moves)")
ax.set_ylabel(r"Interaction energy $U^*/N$ (reduced units)")
ax.set_title(rf"$\rho^*=0.8$, $T^*={tstar:g}$, $n_{{cycl}}$={n1:.0e}")
ax.legend(fontsize=8, loc="lower right")
fig.tight_layout()
fig.savefig(FIG / "U_icycl_rho08.png", dpi=200)
print("Figuras en", FIG)
