from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

DATA_FILE = Path(__file__).with_name("ex4_results.dat")
PLOT_FILE = Path(__file__).with_name("ex4_plot.png")

N, A_average, A_deviation = np.loadtxt(DATA_FILE, unpack=True)
A_exact = 4.0 * (1.0 + np.log(11.0))

fig, (ax_a, ax_err) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

ax_a.plot(N, A_average, "o-", color="C0", label=r"$A_{\mathrm{average}}$")
ax_a.axhline(A_exact, color="C3", linestyle="--", linewidth=1.2, label=r"$A_{\mathrm{exact}}$")
ax_a.set_ylabel("Área")
ax_a.set_title("Estimación Monte Carlo del área")
ax_a.legend()
ax_a.grid(True, alpha=0.3)

ax_err.plot(N, A_deviation, "o-", color="C1", label=r"$A_{\mathrm{deviation}}$")
ax_err.set_xlabel("Número de puntos N")
ax_err.set_ylabel("Desviación RMS")
ax_err.legend()
ax_err.grid(True, alpha=0.3)

fig.tight_layout()
fig.savefig(PLOT_FILE, dpi=150)
plt.show()
print(f"Gráfico guardado en {PLOT_FILE}")
