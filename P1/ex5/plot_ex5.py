from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

# Carga de datos completa (6 columnas: t, x, v, Ek, Ep, ET)
t_001, x_001, _, ek_001, ep_001, etot_001 = np.loadtxt("euler_001.dat", skiprows=1, unpack=True)
t_01,  x_01,  _, ek_01,  ep_01,  etot_01  = np.loadtxt("euler_01.dat",  skiprows=1, unpack=True)
t_02,  x_02,  _, ek_02,  ep_02,  etot_02  = np.loadtxt("euler_02.dat",  skiprows=1, unpack=True)
t_prd, x_prd, _, ek_prd, ep_prd, etot_prd = np.loadtxt("euler_predictor.dat", skiprows=1, unpack=True)
t_ver, x_ver, _, ek_ver, ep_ver, etot_ver = np.loadtxt("verlet.dat",    skiprows=1, unpack=True)

# Solución analítica y cálculo de referencias
omega = np.sqrt(2.0 / 0.200)
t_fine = np.linspace(0.0, 4.0, 1000)
x_exact = 0.05 * np.cos(omega * t_fine)
e_exact = 0.5 * 2.0 * (0.05**2)  # E0 = 0.0025 J

# Errores absolutos de posición
err_001 = np.abs(x_001 - 0.05 * np.cos(omega * t_001))
err_01  = np.abs(x_01  - 0.05 * np.cos(omega * t_01))
err_02  = np.abs(x_02  - 0.05 * np.cos(omega * t_02))
err_prd = np.abs(x_prd - 0.05 * np.cos(omega * t_prd))
err_ver = np.abs(x_ver - 0.05 * np.cos(omega * t_ver))

# Errores absolutos de energía (|Etot - E0|)
err_e_001 = np.abs(etot_001 - e_exact)
err_e_01  = np.abs(etot_01  - e_exact)
err_e_02  = np.abs(etot_02  - e_exact)
err_e_prd = np.abs(etot_prd - e_exact)
err_e_ver = np.abs(etot_ver - e_exact)

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
ax2_pos.plot(t_ver, x_ver, color="#d62728", label="Verlet", lw=1.5)
ax2_pos.set_xlabel("Tiempo $t$ (s)")
ax2_pos.set_ylabel("Posición $x$ (m)")
ax2_pos.set_title(r"Comparación de métodos para $\Delta t = 0.02$ s")
ax2_pos.grid(True, alpha=0.25)
ax2_pos.legend(loc="upper left")

fig1.tight_layout()
fig1.savefig("ex5_pos_plot.png", dpi=200)

# ==========================================
# FIGURA 2: ERRORES DE POSICIÓN (semilogy)
# ==========================================
fig2, (ax1_err, ax2_err) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

ax1_err.semilogy(t_001, err_001, color="#1f77b4", lw=1.2, label=r"Euler $\Delta t = 0.001$")
ax1_err.semilogy(t_01,  err_01,  color="#ff7f0e", lw=1.2, label=r"Euler $\Delta t = 0.01$")
ax1_err.semilogy(t_02,  err_02,  color="#2ca02c", lw=1.2, label=r"Euler $\Delta t = 0.02$")
ax1_err.set_ylabel("Error absoluto en $x$ (m)")
ax1_err.set_title(r"Error de posición en Euler según el paso temporal $\Delta t$")
ax1_err.grid(True, which="both", alpha=0.25)
ax1_err.legend(loc="lower right")

ax2_err.semilogy(t_02,  err_02,  color="#ff7f0e", lw=1.2, label="Euler estándar", alpha=0.8)
ax2_err.semilogy(t_prd, err_prd, color="#9467bd", lw=1.5, label="Euler Predictor")
ax2_err.semilogy(t_ver, err_ver, color="#d62728", lw=1.5, label="Verlet")
ax2_err.set_xlabel("Tiempo $t$ (s)")
ax2_err.set_ylabel("Error absoluto en $x$ (m)")
ax2_err.set_title(r"Comparación del error de posición para $\Delta t = 0.02$ s")
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

# Subplot 3: Verlet
ax3_e.plot(t_ver, ek_ver, color="#1f77b4", label=r"Cinética $E_k$")
ax3_e.plot(t_ver, ep_ver, color="#2ca02c", label=r"Potencial $E_p$")
ax3_e.plot(t_ver, etot_ver, color="#d62728", lw=1.5, label=r"Total $E_{tot}$")
ax3_e.axhline(e_exact, color="black", linestyle="--", lw=1.0, alpha=0.7, label=r"$E_0$ analítica")
ax3_e.set_xlabel("Tiempo $t$ (s)")
ax3_e.set_ylabel("Energía (J)")
ax3_e.set_title(r"Balance energético: Verlet ($\Delta t = 0.02$ s)")
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
ax_tot.plot(t_ver, etot_ver, color="#d62728", label="Verlet", lw=1.5)
ax_tot.axhline(e_exact, color="black", linestyle="--", lw=1.0, alpha=0.7, label=r"$E_0$ analítica")
ax_tot.set_xlabel("Tiempo $t$ (s)")
ax_tot.set_ylabel(r"Energía Total $E_{tot}$ (J)")
ax_tot.set_title(r"Comparación de la conservación de $E_{tot}$ ($\Delta t = 0.02$ s)")
ax_tot.grid(True, alpha=0.25)
ax_tot.legend(loc="upper left")

fig4.tight_layout()
fig4.savefig("ex5_total_energy_plot.png", dpi=200)

# ==========================================
# FIGURA 5: ERROR DE CONSERVACIÓN DE ENERGÍA (semilogy)
# ==========================================
fig5, (ax1_ee, ax2_ee) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

# Subplot 1: Efecto del dt en Euler
ax1_ee.semilogy(t_001, err_e_001, color="#1f77b4", lw=1.2, label=r"Euler $\Delta t = 0.001$")
ax1_ee.semilogy(t_01,  err_e_01,  color="#ff7f0e", lw=1.2, label=r"Euler $\Delta t = 0.01$")
ax1_ee.semilogy(t_02,  err_e_02,  color="#2ca02c", lw=1.2, label=r"Euler $\Delta t = 0.02$")
ax1_ee.set_ylabel(r"$|E_{tot} - E_0|$ (J)")
ax1_ee.set_title(r"Error absoluto de energía en Euler según $\Delta t$")
ax1_ee.grid(True, which="both", alpha=0.25)
ax1_ee.legend(loc="lower right")

# Subplot 2: Comparativa de métodos para dt = 0.02 s
ax2_ee.semilogy(t_02,  err_e_02,  color="#ff7f0e", lw=1.2, label="Euler estándar", alpha=0.8)
ax2_ee.semilogy(t_prd, err_e_prd, color="#9467bd", lw=1.5, label="Euler Predictor")
ax2_ee.semilogy(t_ver, err_e_ver, color="#d62728", lw=1.5, label="Verlet")
ax2_ee.set_xlabel("Tiempo $t$ (s)")
ax2_ee.set_ylabel(r"$|E_{tot} - E_0|$ (J)")
ax2_ee.set_title(r"Comparación del error de energía entre métodos ($\Delta t = 0.02$ s)")
ax2_ee.grid(True, which="both", alpha=0.25)
ax2_ee.legend(loc="lower right")

fig5.tight_layout()
fig5.savefig("ex5_energy_error_plot.png", dpi=200)

plt.show()

print("Gráficos generados:")
print(" - ex5_pos_plot.png")
print(" - ex5_error_plot.png")
print(" - ex5_energy_plot.png")
print(" - ex5_total_energy_plot.png")
print(" - ex5_energy_error_plot.png")