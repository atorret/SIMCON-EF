from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

here = Path(__file__).resolve().parent
t, energy = np.loadtxt(here / "energy-leap.dat", unpack=True)
_, temperature = np.loadtxt(here / "temp-leap.dat", unpack=True)

fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
plots = (
    (energy, "Energía total $E_{tot}$", "#d62728"),
    (temperature, "Temperatura $T$", "#1f77b4"),
)
for ax, (y, label, color) in zip(axes, plots):
    ax.plot(t, y, color=color, lw=1)
    ax.set_ylabel(label)
    ax.grid(True, alpha=0.25)

axes[0].set_title("Evolución temporal del sistema (Lennard-Jones NVE)")
axes[1].set_xlabel("Tiempo")
fig.tight_layout()
fig.savefig(here / "md_results_plot.png", dpi=200)
plt.show()
