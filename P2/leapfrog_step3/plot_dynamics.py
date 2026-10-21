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


# Small canvas: the three required plots sit side by side in the report.
FIGSIZE = (3.4, 2.7)
plt.rcParams.update({"font.size": 9, "legend.fontsize": 8})


def save(fig, name):
    fig.tight_layout()
    fig.savefig(here / name, dpi=250)
    plt.close(fig)


fig, ax = plt.subplots(figsize=FIGSIZE)
ax.plot(t, vacf, color="#1f77b4", lw=1.2)
ax.axhline(0.0, color="0.45", ls="--", lw=0.8)
ax.set_xlim(0.0, 5.0)
ax.set_xlabel(r"$t$ (reduced units)")
ax.set_ylabel(r"$C(t)$")
ax.grid(True, alpha=0.25)
save(fig, "vacf_plot.png")

fig, ax = plt.subplots(figsize=FIGSIZE)
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
ax.set_xlabel(r"$t$ (reduced units)")
ax.set_ylabel(r"$r^{2}(t)$")
ax.legend()
ax.grid(True, alpha=0.25)
save(fig, "msd_plot.png")

fig, ax = plt.subplots(figsize=FIGSIZE)
ax.plot(omega, sdw, color="#2ca02c", lw=1.2)
ax.axhline(0.0, color="0.45", ls="--", lw=0.8)
ax.set_xlabel(r"$\omega$ (reduced units)")
ax.set_ylabel(r"$S(\omega)$")
ax.grid(True, alpha=0.25)
save(fig, "sdw_plot.png")

fig, ax = plt.subplots(figsize=FIGSIZE)
ax.plot(t_d, dinst, color="#9467bd", lw=1.0, label=r"$r^{2}(t)/(6t)$")
ax.axhline(
    d_einstein, color="black", ls="--", lw=1.2, label=rf"Einstein $D={d_einstein:.4f}$"
)
ax.axvspan(tmin, tmax, color="0.85", zorder=0, label="fit window")
ax.set_ylim(0.0, dinst.max() * 1.35)
ax.set_xlabel(r"$t$ (reduced units)")
ax.set_ylabel(r"$D(t)$ (reduced units)")
ax.legend()
ax.grid(True, alpha=0.25)
save(fig, "diffusion_plot.png")
