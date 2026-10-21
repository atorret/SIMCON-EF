from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

here = Path(__file__).resolve().parent
t, energy = np.loadtxt(here / "energy-leap.dat", unpack=True)
_, temperature = np.loadtxt(here / "temp-leap.dat", unpack=True)

# Wide, short canvas: T and E side by side across the report text width.
plt.rcParams.update({"font.size": 9})
fig, axes = plt.subplots(1, 2, figsize=(6.8, 2.1), sharex=True)
plots = (
    (temperature, r"$T$ (reduced units)", "#1f77b4"),
    (energy, r"$E_{\mathrm{tot}}$ (reduced units)", "#d62728"),
)
for ax, (y, label, color) in zip(axes, plots):
    ax.plot(t, y, color=color, lw=0.6)
    ax.set_xlim(0.0, t[-1])
    ax.set_xlabel(r"$t$ (reduced units)")
    ax.set_ylabel(label)
    ax.grid(True, alpha=0.25)

fig.tight_layout()
fig.savefig(here / "md_results_plot.png", dpi=250)
plt.show()
