import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

DATA_FILE = Path(__file__).with_name("ex2_results.dat")
PLOT_FILE = Path(__file__).with_name("ex2_plot.png")

k, pi_aprox = np.loadtxt(DATA_FILE, unpack=True)
error = np.abs(pi_aprox - math.pi)

fig, (ax_pi, ax_err) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

ax_pi.plot(k, pi_aprox, color="C0", label=r"$\pi_{\mathrm{aprox}}$")
ax_pi.axhline(math.pi, color="C3", linestyle="--", linewidth=1.2, label=r"$\pi$")
ax_pi.set_ylabel("Valor de π")
ax_pi.set_title("Convergencia de la serie de Leibniz")
ax_pi.legend()
ax_pi.grid(True, alpha=0.3)

ax_err.plot(k, error, color="C1", label="|π_aprox − π|")
ax_err.set_xlabel("Iteración K")
ax_err.set_ylabel("Error absoluto")
ax_err.set_yscale("log")
ax_err.legend()
ax_err.grid(True, which="both", alpha=0.3)

fig.tight_layout()
fig.savefig(PLOT_FILE, dpi=150)
plt.show()
print(f"Gráfico guardado en {PLOT_FILE}")
