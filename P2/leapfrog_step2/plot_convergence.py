from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

here = Path(__file__).resolve().parent
# t, etot, pkin, pvir, ptail, ptot
data = np.loadtxt(here / "thermo-leap.dat")
t = data[:, 0]
etot = data[:, 1]
ptot = data[:, 5]
n = np.arange(1, len(t) + 1)
etot_run = np.cumsum(etot) / n
ptot_run = np.cumsum(ptot) / n

fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
axes[0].plot(t, etot_run, color="#d62728", lw=1.2)
axes[0].set_ylabel(r"$\langle E_{\mathrm{tot}}\rangle$ (reduced units)")
axes[0].grid(True, alpha=0.25)

axes[1].plot(t, ptot_run, color="#1f77b4", lw=1.2)
axes[1].set_xlabel(r"$t$ (reduced units)")
axes[1].set_ylabel(r"$\langle P_{\mathrm{tot}}\rangle$ (reduced units)")
axes[1].grid(True, alpha=0.25)

fig.tight_layout()
fig.savefig(here / "convergence_plot.png", dpi=150)
plt.show()
