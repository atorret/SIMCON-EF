from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

here = Path(__file__).resolve().parent
# t, ekin, epot, etail, etot, pkin, pvir, ptail, ptot
data = np.loadtxt(here / "thermo-leap.dat")
t = data[:, 0]
ekin = data[:, 1]
epot = data[:, 2]
etail = data[:, 3]
etot = data[:, 4]

fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(t, ekin, color="red", lw=0.7, label="kinetic")
ax.plot(t, epot, color="green", lw=0.7, label="potential")
ax.plot(t, etail, color="purple", lw=0.8, ls="--", label="tail")
ax.plot(t, etot, color="blue", lw=0.7, label="total")
ax.set_xlabel(r"$t$ (reduced units)")
ax.set_ylabel(r"Energy (reduced units)")
ax.set_xlim(t[0], t[-1])
ax.legend()
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(here / "energy_plot.png", dpi=150)
plt.show()
