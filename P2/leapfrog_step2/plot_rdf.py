from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

here = Path(__file__).resolve().parent
r, gr = np.loadtxt(here / "g-leap.dat", unpack=True)

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(r, gr, color="#2ca02c", lw=1.5)
ax.axhline(1.0, color="0.45", ls="--", lw=1, label=r"$g(r) = 1$")
ax.set_xlabel(r"$r$ (unidades reducidas)")
ax.set_ylabel(r"$g(r)$")
ax.set_title("Función de distribución radial")
ax.set_xlim(0.0, float(r.max()))
ax.set_ylim(0.0, 4)
ax.grid(True, alpha=0.25)
ax.legend()
fig.tight_layout()
fig.savefig(here / "rdf_plot.png", dpi=200)
plt.show()
