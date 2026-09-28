import matplotlib.pyplot as plt
import numpy as np

# Cargar los ficheros de resultados (2 columnas: tiempo/paso, valor)
t_en, etot = np.loadtxt("energy-leap.dat", unpack=True)
t_tp, temp = np.loadtxt("temp-leap.dat", unpack=True)

# Crear figura con dos subplots (Energía y Temperatura)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

# Gráfico de Energía Total
ax1.plot(t_en, etot, color="#d62728", lw=1.0)
ax1.set_ylabel("Energía Total $E_{tot}$")
ax1.set_title("Evolución temporal del sistema (Lennard-Jones NVE)")
ax1.grid(True, alpha=0.25)

# Gráfico de Temperatura
ax2.plot(t_tp, temp, color="#1f77b4", lw=1.0)
ax2.set_xlabel("Tiempo / Pasos")
ax2.set_ylabel("Temperatura $T$")
ax2.grid(True, alpha=0.25)

fig.tight_layout()
fig.savefig("md_results_plot.png", dpi=200)
plt.show()

print("Gráfico generado: md_results_plot.png")