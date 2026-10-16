from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

here = Path(__file__).resolve().parent
# t, etot, pkin, pvir, ptail, ptot
data = np.loadtxt(here / "thermo-leap.dat")
t = data[:, 0]
pkin = data[:, 2]
# Configurational pressure: virial plus the constant tail correction.
ppot = data[:, 3] + data[:, 4]
ptot = data[:, 5]

fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(t, pkin, color="red", lw=0.7, label="kinetic")
ax.plot(t, ppot, color="green", lw=0.7, label="potential")
ax.plot(t, ptot, color="blue", lw=0.7, label="total")
ax.set_xlabel("time (reduced units)")
ax.set_ylabel("Pressure (reduced units)")
ax.set_xlim(0, 50)
ax.set_ylim(-1.5, 2)
ax.legend()
fig.tight_layout()
fig.savefig(here / "pressure_plot.png", dpi=100)
plt.show()
