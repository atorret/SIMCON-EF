from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

here = Path(__file__).resolve().parent

# Columnas de thermo-leap.dat:
# 0: t, 1: ekin, 2: epot, 3: etail, 4: etot, 5: pkin, 6: pvir, 7: ptail, 8: ptot
data = np.loadtxt(here / "thermo-leap.dat")
t = data[:, 0]
pkin = data[:, 5]
pvir = data[:, 6]
ptail = data[:, 7]
ptot = data[:, 8]

fig, ax = plt.subplots(figsize=(8, 6))

ax.plot(t, pkin, color="red", lw=0.8, label=r"$P_{\mathrm{kin}}$")
ax.plot(t, pvir, color="green", lw=0.7, alpha=0.85, label=r"$P_{\mathrm{pot}}$ (virial)")
ax.plot(t, ptail, color="purple", lw=0.8, linestyle="--", label=r"$\Delta P$ (tail)")
ax.plot(t, ptot, color="blue", lw=0.8, label=r"$P_{\mathrm{tot}}$")

ax.axhline(0.0, color="black", lw=0.5, linestyle=":", alpha=0.6)
ax.set_xlabel("time (reduced units)")
ax.set_ylabel("Pressure (reduced units)")
ax.set_xlim(t[0], t[-1])
ax.grid(True, alpha=0.3)
ax.legend(loc="upper right")

fig.tight_layout()
fig.savefig(here / "pressure_plot.png", dpi=150)
plt.show()  