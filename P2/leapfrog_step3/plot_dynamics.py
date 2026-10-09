from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

here = Path(__file__).resolve().parent
t, vacf = np.loadtxt(here / "c-leap.dat", unpack=True)
t_msd, msd = np.loadtxt(here / "r2-leap.dat", unpack=True)
omega, sdw = np.loadtxt(here / "c_spectrum-leap.dat", unpack=True)
t_d, dinst = np.loadtxt(here / "D-leap.dat", unpack=True)
d_einstein, slope, intercept, tmin, tmax, d_gk = np.loadtxt(
    here / "diffusion-leap.dat"
)

fig, axes = plt.subplots(2, 2, figsize=(10, 7))

ax = axes[0, 0]
ax.plot(t, vacf, color="#1f77b4", lw=1.2)
ax.axhline(0.0, color="0.45", ls="--", lw=0.8)
ax.set_xlim(0.0, 5.0)
ax.set_xlabel(r"time $t$ (reduced units)")
ax.set_ylabel(r"$C(t) = \langle \vec{v}(t)\cdot\vec{v}(0) \rangle$")
ax.set_title("Velocity autocorrelation")
ax.grid(True, alpha=0.25)

ax = axes[0, 1]
ax.plot(t_msd, msd, color="#d62728", lw=1.2, label="MSD")
t_fit = np.linspace(tmin, tmax, 200)
ax.plot(
    t_fit,
    slope * t_fit + intercept,
    color="black",
    ls="--",
    lw=1.2,
    label=rf"$6Dt+b$, $D={d_einstein:.4f}$",
)
ax.set_xlabel(r"time $t$ (reduced units)")
ax.set_ylabel(r"$r^{2}(t) = \langle [\vec{r}(t)-\vec{r}(0)]^{2} \rangle$")
ax.set_title("Mean square displacement")
ax.legend()
ax.grid(True, alpha=0.25)

ax = axes[1, 0]
ax.plot(omega, sdw, color="#2ca02c", lw=1.2)
ax.axhline(0.0, color="0.45", ls="--", lw=0.8)
ax.set_xlabel(r"frequency $\omega$ (reduced units)")
ax.set_ylabel(r"$S(\omega)$")
ax.set_title("Spectral density of the VACF")
ax.grid(True, alpha=0.25)

ax = axes[1, 1]
ax.plot(t_d, dinst, color="#9467bd", lw=1.0, label=r"$r^{2}(t)/(6t)$")
ax.axhline(d_einstein, color="black", ls="--", lw=1.2, label=rf"Einstein $D={d_einstein:.4f}$")
ax.axvspan(tmin, tmax, color="0.85", zorder=0, label="fit window")
ax.set_ylim(0.0, dinst.max() * 1.35)
ax.set_xlabel(r"time $t$ (reduced units)")
ax.set_ylabel(r"$D(t)$")
ax.set_title("Self-diffusion coefficient")
ax.legend()
ax.grid(True, alpha=0.25)

fig.tight_layout()
fig.savefig(here / "dynamics_plot.png", dpi=200)
plt.show()
