# Cap_02/scripts/plot_fdtd_results.py
import numpy as np
import matplotlib.pyplot as plt
import os

# Diretórios
out_dir = "../out"

# Carrega dados
voltage_final = np.loadtxt(os.path.join(out_dir, "fdtd_voltage.csv"))
voltage_series = np.loadtxt(os.path.join(out_dir, "fdtd_time_series.csv"), delimiter=",")
spectrum = np.loadtxt(os.path.join(out_dir, "fdtd_spectrum.csv"), delimiter=",")


#data = csvread('comparison_voltage.csv', 1, 0); % pula cabeçalho

#plot(data(:,1), data(:,2), '-', data(:,1), data(:,3), '--', ...
#     data(:,4), data(:,5), 'o', data(:,4), data(:,6), '+');
#legend('Real, exact','Imag, exact','Real, FDTD','Imag, FDTD', 'Location', 'Best');
#xlabel('z (m)');
#ylabel('Steady-state voltage (V)');



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
