import numpy as np
import matplotlib.pyplot as plt
import os

# Nombres de los archivos generados por el programa Fortran
energy_file = "../simulation/energy-leap.dat"
temp_file = "../simulation/temp-leap.dat"

# Cargar los datos
# energy_data[:, 0] es el tiempo, energy_data[:, 1] es la energía total
energy_data = np.loadtxt(energy_file)
time_energy = energy_data[:, 0]
total_energy = energy_data[:, 1]

# temp_data[:, 0] es el tiempo, temp_data[:, 1] es la temperatura
temp_data = np.loadtxt(temp_file)
time_temp = temp_data[:, 0]
temperature = temp_data[:, 1]

# Crear la figura con 2 subgráficos (2 filas, 1 columna)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

# Plot de la Energía Total
ax1.plot(time_energy, total_energy, color="tab:blue", linewidth=1.5)
ax1.set_ylabel("Total Energy (reduced units)", fontsize=12)
ax1.set_title("NVE Molecular Dynamics: Energy and Temperature vs Time", fontsize=14)
ax1.grid(True, linestyle="--", alpha=0.7)

# Plot de la Temperatura
ax2.plot(time_temp, temperature, color="tab:red", linewidth=1.5)
ax2.set_ylabel("Temperature (reduced units)", fontsize=12)
ax2.set_xlabel("Time (reduced units)", fontsize=12)
ax2.grid(True, linestyle="--", alpha=0.7)

# Ajustar el diseño para que no se superpongan los textos
plt.tight_layout()

# Guardar la figura como imagen (opcional)
plt.savefig("md_results.png", dpi=300)

# Mostrar el gráfico en pantalla
plt.show()
