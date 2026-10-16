import numpy as np
import matplotlib.pyplot as plt
import os

# Nombres de los archivos generados por Fortran
file_rdf = "../simulation/rdf_output.dat"
file_energy = "../simulation/energy-leap.dat"
file_pressure = "../simulation/pressure.dat"

# ---------------------------------------------------------
# 1. Carga de datos
# ---------------------------------------------------------
# RDF: r (col 0), g(r) (col 1)
data_rdf = np.loadtxt(file_rdf)
r = data_rdf[:, 0]
g_r = data_rdf[:, 1]

# Energía: tiempo (col 0), energía total (col 1)
data_energy = np.loadtxt(file_energy)
time = data_energy[:, 0]
energy = data_energy[:, 1]

# Presión: tiempo (col 0), presión (col 1)
# (Se asume que la longitud de tiempo es la misma que la de energía)
data_pressure = np.loadtxt(file_pressure)
pressure = data_pressure[:, 1]

# ---------------------------------------------------------
# 2. Análisis de Convergencia (Medias acumuladas)
# ---------------------------------------------------------
# La media acumulada nos dice cómo el valor promedio se estabiliza con el tiempo
energy_running_avg = np.cumsum(energy) / np.arange(1, len(energy) + 1)
pressure_running_avg = np.cumsum(pressure) / np.arange(1, len(pressure) + 1)

# Imprimir en consola los resultados finales (promedios en el equilibrio)
# Nota: Para un análisis de equilibrio estricto, normalmente se descarta
# la fase de "termalización" (por ejemplo, el primer 20% de los datos).
term_cut = int(
    0.2 * len(energy)
)  # Descartamos el primer 20% para el promedio final riguroso
mean_energy = np.mean(energy[term_cut:])
mean_pressure = np.mean(pressure[term_cut:])

print("-" * 50)
print("MEAN RESULTS (excluding initial 20% of thermalization):")
print(f"Total Energy Mean (E) : {mean_energy:.4f} (Reduced Units)")
print(f"Pressure Mean (P)       : {mean_pressure:.4f} (Reduced Units)")
print("-" * 50)

# ---------------------------------------------------------
# 3. Creación de los gráficos
# ---------------------------------------------------------
fig, axs = plt.subplots(1, 3, figsize=(16, 5))

# Gráfico 1: Función de Distribución Radial g(r)
axs[0].plot(r, g_r, color="purple", lw=2)
axs[0].axhline(
    y=1.0, color="black", linestyle="--", alpha=0.6
)  # Límite asintótico de g(r)
axs[0].set_title("Radial Distribution Function, g(r)")
axs[0].set_xlabel("Distance r (reduced units)")
axs[0].set_ylabel("g(r)")
axs[0].grid(True, alpha=0.3)
axs[0].set_xlim(0, max(r))

# Gráfico 2: Convergencia de la Energía Total
axs[1].plot(time, energy, color="tab:blue", alpha=0.3, label="Instantaneous Value")
axs[1].plot(time, energy_running_avg, color="navy", lw=2, label="Running Average")
axs[1].axhline(y=mean_energy, color="red", linestyle="--", label="Mean (equilibrium)")
axs[1].set_title("Convergence of Total Energy")
axs[1].set_xlabel("Time (reduced units)")
axs[1].set_ylabel("Total Energy (reduced units)")
axs[1].legend()
axs[1].grid(True, alpha=0.3)

# Gráfico 3: Convergencia de la Presión
axs[2].plot(time, pressure, color="tab:green", alpha=0.3, label="Instantaneous Value")
axs[2].plot(
    time, pressure_running_avg, color="darkgreen", lw=2, label="Running Average"
)
axs[2].axhline(y=mean_pressure, color="red", linestyle="--", label="Mean (equilibrium)")
axs[2].set_title("Convergence of Pressure")
axs[2].set_xlabel("Time (reduced units)")
axs[2].set_ylabel("Pressure P (reduced units)")
axs[2].legend()
axs[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("md_analysis_results.png", dpi=300)
plt.show()
