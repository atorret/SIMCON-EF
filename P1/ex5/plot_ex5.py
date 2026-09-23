from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

# Carga de datos completa (5 columnas: t, x, v, Ek, Ep)
t_001, x_001, _, ek_001, ep_001 = np.loadtxt("euler_001.dat", skiprows=1, unpack=True)
t_01,  x_01,  _, ek_01,  ep_01  = np.loadtxt("euler_01.dat",  skiprows=1, unpack=True)
t_02,  x_02,  _, ek_02,  ep_02  = np.loadtxt("euler_02.dat",  skiprows=1, unpack=True)
t_prd, x_prd, _, ek_prd, ep_prd = np.loadtxt("euler_predictor.dat", skiprows=1, unpack=True)
t_ver, x_ver, _, ek_ver, ep_ver = np.loadtxt("verlet.dat",    skiprows=1, unpack=True)

# Energías totales (E_tot = Ek + Ep)
etot_001 = ek_001 + ep_001
etot_01  = ek_01  + ep_01
etot_02  = ek_02  + ep_02
etot_prd = ek_prd + ep_prd
etot_ver = ek_ver + ep_ver

# Solución analítica y cálculo de errores de posición
omega = np.sqrt(2.0 / 0.200)
t_fine = np.linspace(0.0, 4.0, 1000)
x_exact = 0.05 * np.cos(omega * t_fine)
e_exact = 0.5 * 2.0 * (0.05**2)  # E0 = 1/2 * k * x0^2 = 0.0025 J

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
# FIGURA 3: BALANCES ENERGÉTICOS POR MÉTODO (dt = 0.02 s)
# ==========================================
fig3, (ax1_e, ax2_e, ax3_e) = plt.subplots(3, 1, figsize=(8, 9), sharex=True)

# Subplot 1: Euler estándar
ax1_e.plot(t_02, ek_02, color="#1f77b4", label=r"Cinética $E_k$")
ax1_e.plot(t_02, ep_02, color="#2ca02c", label=r"Potencial $E_p$")
ax1_e.plot(t_02, etot_02, color="#d62728", lw=1.5, label=r"Total $E_{tot}$")
ax1_e.axhline(e_exact, color="black", linestyle="--", lw=1.0, alpha=0.7, label=r"$E_0$ analítica")
ax1_e.set_ylabel("Energía (J)")
ax1_e.set_title(r"Balance energético: Euler estándar ($\Delta t = 0.02$ s)")
ax1_e.grid(True, alpha=0.25)
ax1_e.legend(loc="upper left")

# Subplot 2: Euler Predictor
ax2_e.plot(t_prd, ek_prd, color="#1f77b4", label=r"Cinética $E_k$")
ax2_e.plot(t_prd, ep_prd, color="#2ca02c", label=r"Potencial $E_p$")
ax2_e.plot(t_prd, etot_prd, color="#d62728", lw=1.5, label=r"Total $E_{tot}$")
ax2_e.axhline(e_exact, color="black", linestyle="--", lw=1.0, alpha=0.7, label=r"$E_0$ analítica")
ax2_e.set_ylabel("Energía (J)")
ax2_e.set_title(r"Balance energético: Euler Predictor ($\Delta t = 0.02$ s)")
ax2_e.grid(True, alpha=0.25)
ax2_e.legend(loc="upper left")

# Subplot 3: Velocity-Verlet
ax3_e.plot(t_ver, ek_ver, color="#1f77b4", label=r"Cinética $E_k$")
ax3_e.plot(t_ver, ep_ver, color="#2ca02c", label=r"Potencial $E_p$")
ax3_e.plot(t_ver, etot_ver, color="#d62728", lw=1.5, label=r"Total $E_{tot}$")
ax3_e.axhline(e_exact, color="black", linestyle="--", lw=1.0, alpha=0.7, label=r"$E_0$ analítica")
ax3_e.set_xlabel("Tiempo $t$ (s)")
ax3_e.set_ylabel("Energía (J)")
ax3_e.set_title(r"Balance energético: Velocity-Verlet ($\Delta t = 0.02$ s)")
ax3_e.grid(True, alpha=0.25)
ax3_e.legend(loc="upper left")

fig3.tight_layout()
fig3.savefig("ex5_energy_plot.png", dpi=200)

# ==========================================
# FIGURA 4: COMPARACIÓN DIRECTA DE ENERGÍA TOTAL
# ==========================================
fig4, ax_tot = plt.subplots(figsize=(8, 4.5))

ax_tot.plot(t_02,  etot_02,  color="#ff7f0e", label="Euler estándar", alpha=0.8)
ax_tot.plot(t_prd, etot_prd, color="#9467bd", label="Euler Predictor", lw=1.5)
ax_tot.plot(t_ver, etot_ver, color="#d62728", label="Velocity-Verlet", lw=1.5)
ax_tot.axhline(e_exact, color="black", linestyle="--", lw=1.0, alpha=0.7, label=r"$E_0$ analítica")
ax_tot.set_xlabel("Tiempo $t$ (s)")
ax_tot.set_ylabel(r"Energía Total $E_{tot}$ (J)")
ax_tot.set_title(r"Comparación de la conservación de $E_{tot}$ ($\Delta t = 0.02$ s)")
ax_tot.grid(True, alpha=0.25)
ax_tot.legend(loc="upper left")

fig4.tight_layout()
fig4.savefig("ex5_total_energy_plot.png", dpi=200)

# ==========================================
# MOSTRAR TODAS LAS FIGURAS
# ==========================================
plt.show()

print("Gráficos generados:")
print(" - ex5_pos_plot.png")
print(" - ex5_error_plot.png")
print(" - ex5_energy_plot.png (desglose Ek, Ep, Etot para Euler, Predictor y Verlet)")
print(" - ex5_total_energy_plot.png (comparación directa de Etot)")
