# Cap_02/scripts/plot_fdtd_results.py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os

# Diretórios
out_dir = "../out"

comparison_voltage = pd.read_csv(os.path.join(out_dir, 'comparison_voltage.csv'))
parametros = pd.read_csv(os.path.join(out_dir, 'simulation_parameters.csv'))

#print("\nParâmetros da Simulação:")
#print(parametros.to_string(index=False))

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

param_text = "\n".join(f"{row['Parameter']}: {row['Value']}" for _, row in parametros.iterrows())
plt.gcf().text(0.9, 0.2, param_text, ha='right', va='bottom', bbox=dict(facecolor='white', alpha=0.8))

# Salva o gráfico como arquivo PNG
plt.savefig(os.path.join(out_dir, 'comparison_voltage_plot.png'))

# Exibe o gráfico
#plt.show()

# Lê a tabela ignorando as duas primeiras linhas
df = pd.read_csv(os.path.join(out_dir, 'erro_relativo.csv'), skiprows=2)

# Captura a norma L2 da primeira linha
with open(os.path.join(out_dir, 'erro_relativo.csv'), 'r') as f:
    first_line = f.readline().strip()
    norma_l2 = float(first_line.split(':')[1].strip())

# Cria o gráfico
plt.figure(figsize=(10, 6))
plt.plot(df['Index'], df['Erro_relativo'], marker='o')
plt.title(f'Erro Relativo ponto a ponto\nNorma L2 Relativa: {norma_l2:.6f}') # Adiciona a norma L2 no topo
plt.xlabel('Index (z)')
plt.ylabel('Erro Relativo (%)')

plt.grid(True)
plt.savefig(os.path.join(out_dir, 'erro_relativo.png'))

#plt.show()
# Lê o arquivo de tensões
voltage_series = pd.read_csv(os.path.join(out_dir, 'voltage_over_time.csv'))

fig, ax = plt.subplots()
line, = ax.plot([], [], lw=2)
ax.set_xlim(1, voltage_series.shape[1]-1)
ax.set_ylim(-1, 1)
ax.set_xlabel('Index (z)')
ax.set_ylabel('Normalized Voltage')

def init():
    line.set_data([], [])
    return line,

def update(frame):
    x = range(1, voltage_series.shape[1])  # posição espacial
    y = voltage_series.iloc[frame, 1:]     # ignora a coluna TimeStep
    line.set_data(x, y)
    ax.set_title(f'Voltage at timestep {voltage_series.iloc[frame, 0]}')
    return line,

ani = FuncAnimation(fig, update, frames=len(voltage_series), init_func=init, blit=True)

# Salva o vídeo
ani.save(os.path.join(out_dir, 'voltage_simulation.mp4'), writer='ffmpeg')
ani.save(os.path.join(out_dir, 'voltage_simulation.gif'), writer='pillow')
plt.close(fig)  # fecha apenas depois de salvar

# Lê o arquivo CSV da corrente
current_series = pd.read_csv(os.path.join(out_dir, 'current_over_time.csv'))

fig, ax = plt.subplots()
line, = ax.plot([], [], lw=2)
ax.set_xlim(1, current_series.shape[1]-1)
ax.set_ylim(-1, 1)  # ajuste se os valores da corrente tiverem outra faixa
ax.set_xlabel('Index (z)')
ax.set_ylabel('Current (A)')

def init():
    line.set_data([], [])
    return line,

def update(frame):
    x = range(1, current_series.shape[1])  # posição espacial
    y = current_series.iloc[frame, 1:]     # ignora a coluna TimeStep
    line.set_data(x, y)
    ax.set_title(f'Current at timestep {int(current_series.iloc[frame, 0])}')
    return line,

ani = FuncAnimation(fig, update, frames=len(current_series), init_func=init, blit=True)

# Salva o GIF
ani.save(os.path.join(out_dir, 'current_simulation.mp4'), writer='ffmpeg')
ani.save(os.path.join(out_dir, 'current_simulation.gif'), writer='pillow')

plt.close(fig)  # fecha apenas depois de salvar