from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

here = Path(__file__).resolve().parent
r, gr = np.loadtxt(here / "g-leap.dat", unpack=True)

# Small canvas: the three Exercise 2 plots sit side by side in the report.
plt.rcParams.update({"font.size": 9, "legend.fontsize": 8})
fig, ax = plt.subplots(figsize=(3.4, 2.7))
ax.plot(r, gr, color="#2ca02c", lw=1.2)
ax.axhline(1.0, color="0.45", ls="--", lw=0.8, label=r"$g(r) = 1$")
ax.set_xlabel(r"$r$ (reduced units)")
ax.set_ylabel(r"$g(r)$")
ax.set_xlim(0.0, float(r.max()))
ax.set_ylim(0.0, 3.5)
ax.grid(True, alpha=0.25)
ax.legend()
fig.tight_layout()
fig.savefig(here / "rdf_plot.png", dpi=250)
plt.show()
