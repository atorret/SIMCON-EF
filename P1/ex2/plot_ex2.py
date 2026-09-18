from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

# Rutas de archivos
data_file = Path(__file__).with_name("ex2_results.dat")
plot_file = Path(__file__).with_name("ex2_plot.png")

# Carga de datos y cálculo de error
k, pi_aprox = np.loadtxt(data_file, unpack=True)
error = np.abs(pi_aprox - np.pi)

# Configuración del lienzo
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 6), sharex=True)

# Panel 1: Aproximación de Pi
ax1.plot(k, pi_aprox, color="#1f77b4", lw=1.2, label=r"$\pi_{\mathrm{aprox}}$")
ax1.axhline(np.pi, color="#d62728", ls="--", lw=1.2, label=r"$\pi$ exacto")
ax1.set_title("Convergencia de la serie de Leibniz", fontsize=11)
ax1.set_ylabel(r"Valor de $\pi$")
ax1.grid(True, alpha=0.25)
ax1.legend(loc="upper right")

# Panel 2: Error absoluto
ax2.plot(k, error, color="#ff7f0e", lw=1.2, label=r"$|\pi_{\mathrm{aprox}} - \pi|$")
ax2.set_yscale("log")
ax2.set_xlabel("Iteración (K)")
ax2.set_ylabel("Error")
ax2.grid(True, which="both", alpha=0.25)
ax2.legend(loc="upper right")

# Guardar y mostrar
fig.tight_layout()
fig.savefig(plot_file, dpi=200)
plt.show()