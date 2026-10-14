from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

here = Path(__file__).resolve().parent
t, energy = np.loadtxt(here / "energy-leap.dat", unpack=True)
_, temperature = np.loadtxt(here / "temp-leap.dat", unpack=True)

fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
plots = (
    (energy, r"$E_{\mathrm{tot}}$ (reduced units)", "#d62728"),
    (temperature, r"$T$ (reduced units)", "#1f77b4"),
)
for ax, (y, label, color) in zip(axes, plots):
    ax.plot(t, y, color=color, lw=1)
    ax.set_ylabel(label)
    ax.grid(True, alpha=0.25)

axes[0].set_title("Time evolution of the Lennard-Jones NVE liquid")
axes[1].set_xlabel(r"$t$ (reduced units)")
fig.tight_layout()
fig.savefig(here / "md_results_plot.png", dpi=200)
plt.show()
