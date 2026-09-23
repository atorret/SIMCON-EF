from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

# Carga de datos
t_001, x_001, _ = np.loadtxt("euler_001.dat", skiprows=1, unpack=True)
t_01,  x_01,  _ = np.loadtxt("euler_01.dat",  skiprows=1, unpack=True)
t_02,  x_02,  _ = np.loadtxt("euler_02.dat",  skiprows=1, unpack=True)
t_prd, x_prd, _ = np.loadtxt("euler_predictor.dat", skiprows=1, unpack=True)
t_ver, x_ver, _ = np.loadtxt("verlet.dat",    skiprows=1, unpack=True)

# Solución analítica y cálculo de errores
omega = np.sqrt(2.0 / 0.200)
t_fine = np.linspace(0.0, 4.0, 1000)
x_exact = 0.05 * np.cos(omega * t_fine)

err_001 = np.abs(x_001 - 0.05 * np.cos(omega * t_001))
err_01  = np.abs(x_01  - 0.05 * np.cos(omega * t_01))
err_02  = np.abs(x_02  - 0.05 * np.cos(omega * t_02))
err_prd = np.abs(x_prd - 0.05 * np.cos(omega * t_prd))
err_ver = np.abs(x_ver - 0.05 * np.cos(omega * t_ver))

# ==========================================
# FIGURA 1: POSICIONES
# ==========================================
fig1, (ax1_pos, ax2_pos) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

ax1_pos.plot(t_fine, x_exact, "k--", lw=1.2, label="Analítica", alpha=0.7)
ax1_pos.plot(t_001, x_001, color="#1f77b4", label=r"Euler $\Delta t=0.001$")
ax1_pos.plot(t_01,  x_01,  color="#ff7f0e", label=r"Euler $\Delta t=0.01$")
ax1_pos.plot(t_02,  x_02,  color="#2ca02c", label=r"Euler $\Delta t=0.02$")
ax1_pos.set_ylabel("Posición $x$ (m)")
ax1_pos.set_title(r"Influencia de $\Delta t$ en la estabilidad de Euler")
ax1_pos.grid(True, alpha=0.25)
ax1_pos.legend(loc="upper left")

ax2_pos.plot(t_fine, x_exact, "k--", lw=1.2, label="Analítica", alpha=0.7)
ax2_pos.plot(t_02,  x_02,  color="#ff7f0e", label="Euler estándar", alpha=0.8)
ax2_pos.plot(t_prd, x_prd, color="#9467bd", label="Euler Predictor", lw=1.5)
ax2_pos.plot(t_ver, x_ver, color="#d62728", label="Velocity-Verlet", lw=1.5)
ax2_pos.set_xlabel("Tiempo $t$ (s)")
ax2_pos.set_ylabel("Posición $x$ (m)")
ax2_pos.set_title(r"Comparación de métodos para $\Delta t = 0.02$ s")
ax2_pos.grid(True, alpha=0.25)
ax2_pos.legend(loc="upper left")

fig1.tight_layout()
fig1.savefig("ex5_pos_plot.png", dpi=200)

# ==========================================
# FIGURA 2: ERRORES (semilogy)
# ==========================================
fig2, (ax1_err, ax2_err) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

ax1_err.semilogy(t_001, err_001, color="#1f77b4", lw=1.2, label=r"Euler $\Delta t = 0.001$")
ax1_err.semilogy(t_01,  err_01,  color="#ff7f0e", lw=1.2, label=r"Euler $\Delta t = 0.01$")
ax1_err.semilogy(t_02,  err_02,  color="#2ca02c", lw=1.2, label=r"Euler $\Delta t = 0.02$")
ax1_err.set_ylabel("Error absoluto (m)")
ax1_err.set_title(r"Error del método de Euler según el paso temporal $\Delta t$")
ax1_err.grid(True, which="both", alpha=0.25)
ax1_err.legend(loc="lower right")

ax2_err.semilogy(t_02,  err_02,  color="#ff7f0e", lw=1.2, label="Euler estándar", alpha=0.8)
ax2_err.semilogy(t_prd, err_prd, color="#9467bd", lw=1.5, label="Euler Predictor")
ax2_err.semilogy(t_ver, err_ver, color="#d62728", lw=1.5, label="Verlet")
ax2_err.set_xlabel("Tiempo $t$ (s)")
ax2_err.set_ylabel("Error absoluto (m)")
ax2_err.set_title(r"Comparación del error numérico para $\Delta t = 0.02$ s")
ax2_err.grid(True, which="both", alpha=0.25)
ax2_err.legend(loc="lower right")

fig2.tight_layout()
fig2.savefig("ex5_error_plot.png", dpi=200)

# ==========================================
# MOSTRAR AMBAS FIGURAS
# ==========================================
plt.show()

<<<<<<< HEAD
print("Plot generated and saved as ex5_plot.png")
=======
print("Gráficos guardados como ex5_pos_plot.png y ex5_error_plot.png")
>>>>>>> 292a06963409cbae393f448cf7a2a9ec18e331a6
