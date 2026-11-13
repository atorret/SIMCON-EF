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
        ax.plot(rho, rho * tstar, "k--", label=r"gas ideal $\rho^*T^*$")
    else:
        ax.axhline(0, color="k", ls="--", label="gas ideal")
    ax.set_xlabel(r"$\rho^*$")
    ax.set_ylabel(ylab)
    ax.set_title(rf"$T^*={tstar:g}$ (con correcciones de cola)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG / f"{name}.png", dpi=200)
    plt.close(fig)

n = max(series)
fig, ax = plt.subplots(figsize=(6, 4))
for d in sorted((RUNS / f"n{n}").glob("rho*")):
    r, g = np.loadtxt(d / "gdr-mc.dat", unpack=True)
    ax.plot(r, g, label=rf"$\rho^*={float(d.name[3:]):.1f}$")
ax.set_xlabel(r"$r^*$")
ax.set_ylabel(r"$g(r)$")
ax.set_title(rf"$T^*={tstar:g}$, $n_{{cycl}}$={n:.0e}")
ax.legend(fontsize=7, ncol=2)
fig.tight_layout()
fig.savefig(FIG / "gdr.png", dpi=200)
print("Figuras en", FIG)
