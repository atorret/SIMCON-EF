from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

# Carga de archivos .dat omitiendo la cabecera de texto con skiprows=1
t_001, x_001, _ = np.loadtxt("euler_001.dat", skiprows=1, unpack=True)
t_01,  x_01,  _ = np.loadtxt("euler_01.dat",  skiprows=1, unpack=True)
t_02,  x_02,  _ = np.loadtxt("euler_02.dat",  skiprows=1, unpack=True)
t_prd, x_prd, _ = np.loadtxt("euler_predictor.dat", skiprows=1, unpack=True)
t_ver, x_ver, _ = np.loadtxt("verlet.dat",    skiprows=1, unpack=True)

# Solución analítica del oscilador armónico: x(t) = x_0 * cos(omega * t)
# omega = sqrt(k / m) = sqrt(2.0 / 0.2) = sqrt(10) ≈ 3.1622777
omega = np.sqrt(2.0 / 0.200)
t_fine = np.linspace(0.0, 6.0, 1000)
x_exact = 0.05 * np.cos(omega * t_fine)

# Configuración de la figura con dos subgráficas en columna
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

# Panel 1: Efecto del cambio de paso de tiempo dt en el método de Euler
ax1.plot(t_fine, x_exact, color="black", linestyle="--", linewidth=1.2, label="Analítica", alpha=0.7)
ax1.plot(t_001,  x_001,  color="#1f77b4", linewidth=1.0, label=r"Euler $\Delta t = 0.001$")
ax1.plot(t_01,   x_01,   color="#ff7f0e", linewidth=1.2, label=r"Euler $\Delta t = 0.01$")
ax1.plot(t_02,   x_02,   color="#2ca02c", linewidth=1.2, label=r"Euler $\Delta t = 0.02$")
ax1.set_ylabel("Posición $x$ (m)")
ax1.set_title("Influencia del paso temporal $\Delta t$ en el método de Euler")
ax1.grid(True, alpha=0.25)
ax1.legend(loc="upper left")

# Panel 2: Comparativa de integradores para el mismo paso dt = 0.02 s
ax2.plot(t_fine, x_exact, color="black", linestyle="--", linewidth=1.2, label="Analítica", alpha=0.7)
ax2.plot(t_02,   x_02,   color="#ff7f0e", linewidth=1.2, label="Euler estándar", alpha=0.8)
ax2.plot(t_prd,  x_prd,  color="#9467bd", linewidth=1.5, label="Euler Predictor")
ax2.plot(t_ver,  x_ver,  color="#d62728", linewidth=1.5, label="Verlet")
ax2.set_xlabel("Tiempo $t$ (s)")
ax2.set_ylabel("Posición $x$ (m)")
ax2.set_title(r"Comparación de métodos numéricos para $\Delta t = 0.02$ s")
ax2.grid(True, alpha=0.25)
ax2.legend(loc="upper left")

# Ajustes de diseño, guardado y despliegue
fig.tight_layout()
fig.savefig("ex5_plot.png", dpi=200)
plt.show()

print("Gráfico generado y guardado como ex5_plot.png")