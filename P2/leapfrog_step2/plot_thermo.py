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

# Small canvas: the three Exercise 2 plots sit side by side in the report.
plt.rcParams.update({"font.size": 9, "legend.fontsize": 8})
fig, ax = plt.subplots(figsize=(3.4, 2.7))
ax.plot(t, pkin, color="red", lw=0.5, label="kinetic")
ax.plot(t, ppot, color="green", lw=0.5, label="potential")
ax.plot(t, ptot, color="blue", lw=0.5, label="total")
ax.set_xlabel(r"$t$ (reduced units)")
ax.set_ylabel(r"$P$ (reduced units)")
ax.set_xlim(0, 50)
ax.set_ylim(-1.5, 2)
ax.legend(ncol=3, loc="upper center", handlelength=1.2, columnspacing=0.8)
fig.tight_layout()
fig.savefig(here / "pressure_plot.png", dpi=250)
plt.show()
