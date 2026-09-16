from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

SIGNAL_FILE = Path(__file__).with_name("ex3_data.dat")
RESULT_FILE = Path(__file__).with_name("ex3_results.dat")
PLOT_FILE = Path(__file__).with_name("ex3_plot.png")

t, C = np.loadtxt(SIGNAL_FILE, unpack=True)
w, F = np.loadtxt(RESULT_FILE, unpack=True)

fig, (ax_c, ax_f) = plt.subplots(2, 1, figsize=(8, 7))

ax_c.plot(t, C, color="C0", label=r"$C(t)$")
ax_c.set_xlabel("t")
ax_c.set_ylabel("C(t)")
ax_c.set_title("Señal de entrada")
ax_c.legend()
ax_c.grid(True, alpha=0.3)

ax_f.plot(w, F, color="C1", label=r"$F(\omega)=\int C(t)\cos(\omega t)\,dt$")
ax_f.set_xlabel(r"Frecuencia $\omega$")
ax_f.set_ylabel(r"$F(\omega)$")
ax_f.set_title("Transformada de Fourier en coseno (Simpson)")
ax_f.legend()
ax_f.grid(True, alpha=0.3)

fig.tight_layout()
fig.savefig(PLOT_FILE, dpi=150)
plt.show()
print(f"Gráfico guardado en {PLOT_FILE}")
