# Cap_02/scripts/plot_fdtd_results.py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# Diretórios
out_dir = "../out"

# Carrega dados
voltage_final = np.loadtxt(os.path.join(out_dir, "fdtd_voltage.csv"))
voltage_series = np.loadtxt(os.path.join(out_dir, "fdtd_time_series.csv"), delimiter=",")
spectrum = np.loadtxt(os.path.join(out_dir, "fdtd_spectrum.csv"), delimiter=",")

comparison_voltage = pd.read_csv(os.path.join(out_dir, 'comparison_voltage.csv'))

# Extrai os dados
z_exact = comparison_voltage['z_exact']
V_exact_real = comparison_voltage['Re(V_exact)']
V_exact_imag = comparison_voltage['Im(V_exact)']
z_fdtd = comparison_voltage['z_fdtd']
V_fdtd_real = comparison_voltage['Re(V_fdtd)']
V_fdtd_imag = comparison_voltage['Im(V_fdtd)']

# Cria o gráfico
plt.figure(figsize=(8, 6))

# Linha cheia: parte real exata
plt.plot(z_exact, V_exact_real, '-', label='Real, exact')

# Linha tracejada: parte imaginária exata
plt.plot(z_exact, V_exact_imag, '--', label='Imag, exact')

# Círculo transparente: parte real FDTD
plt.plot(z_fdtd, V_fdtd_real, 'o', mfc='none', label='Real, FDTD')

# Sinal de +: parte imaginária FDTD
plt.plot(z_fdtd, V_fdtd_imag, '+', label='Imag, FDTD')

# Configurações dos eixos e da legenda
plt.xlabel('z (m)')
plt.ylabel('Steady-state voltage (V)')
plt.legend()
plt.grid(True)
plt.tight_layout()

# Salva o gráfico como arquivo PNG
plt.savefig(os.path.join(out_dir, 'comparison_voltage_plot.png'), dpi=300)

# Exibe o gráfico
plt.show()



# Parâmetros
Nz = voltage_final.shape[0]
Nk = voltage_series.shape[0]
dz = 0.25 / (Nz - 1)
z = np.linspace(0, 0.25, Nz)
t = np.arange(Nk)

# --- Gráfico 1: tensão final ---
plt.figure()
plt.plot(z, voltage_final)
plt.title("Tensão final no tempo (último passo)")
plt.xlabel("z (m)")
plt.ylabel("V(z, t_final)")
plt.grid(True)
plt.savefig(os.path.join(out_dir, "voltage_final.png"))

# --- Gráfico 2: mapa de calor V(z,t) ---
plt.figure()
plt.imshow(voltage_series, aspect='auto', origin='lower',
           extent=[0, 0.25, 0, Nk], cmap='jet')
plt.colorbar(label="Tensão (V)")
plt.xlabel("z (m)")
plt.ylabel("Tempo (passos)")
plt.title("Evolução temporal da tensão V(z,t)")
plt.savefig(os.path.join(out_dir, "voltage_heatmap.png"))

# --- Gráfico 3: espectro (FFT componente k=2) ---
real = spectrum[:, 0]
imag = spectrum[:, 1]
magnitude = np.sqrt(real**2 + imag**2)
phase = np.arctan2(imag, real)

plt.figure()
plt.plot(z, magnitude)
plt.title("Módulo da FFT na componente k=2")
plt.xlabel("z (m)")
plt.ylabel("|V_k=2(z)|")
plt.grid(True)
plt.savefig(os.path.join(out_dir, "fft_magnitude.png"))

plt.figure()
plt.plot(z, phase)
plt.title("Fase da FFT na componente k=2")
plt.xlabel("z (m)")
plt.ylabel("Fase (rad)")
plt.grid(True)
plt.savefig(os.path.join(out_dir, "fft_phase.png"))

#plt.show()
