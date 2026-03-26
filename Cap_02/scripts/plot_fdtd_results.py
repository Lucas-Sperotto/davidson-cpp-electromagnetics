import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from pathlib import Path
import csv

OUT_DIR = Path(__file__).resolve().parents[1] / "out"


def read_dict_rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


comparison_rows = read_dict_rows(OUT_DIR / "comparison_voltage.csv")
parametros = read_dict_rows(OUT_DIR / "simulation_parameters.csv")

z_exact = np.array([float(row["z_exact"]) for row in comparison_rows if row["z_exact"]], dtype=float)
V_exact_real = np.array([float(row["Re(V_exact)"]) for row in comparison_rows if row["z_exact"]], dtype=float)
V_exact_imag = np.array([float(row["Im(V_exact)"]) for row in comparison_rows if row["z_exact"]], dtype=float)
z_fdtd = np.array([float(row["z_fdtd"]) for row in comparison_rows if row["z_fdtd"]], dtype=float)
V_fdtd_real = np.array([float(row["Re(V_fdtd)"]) for row in comparison_rows if row["z_fdtd"]], dtype=float)
V_fdtd_imag = np.array([float(row["Im(V_fdtd)"]) for row in comparison_rows if row["z_fdtd"]], dtype=float)

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

param_text = "\n".join(f"{row['Parameter']}: {row['Value']}" for row in parametros)
plt.gcf().text(0.9, 0.2, param_text, ha='right', va='bottom', bbox=dict(facecolor='white', alpha=0.8))

# Salva o gráfico como arquivo PNG
plt.savefig(OUT_DIR / "comparison_voltage_plot.png")

with (OUT_DIR / "erro_relativo.csv").open() as handle:
    first_line = handle.readline().strip()
    handle.readline()
    reader = csv.DictReader(handle)
    erro_rows = list(reader)

norma_l2 = float(first_line.split(":")[1].strip())
erro_index = np.array([float(row["Index"]) for row in erro_rows], dtype=float)
erro_relativo = np.array([float(row["Erro_relativo"]) for row in erro_rows], dtype=float)

# Cria o gráfico
plt.figure(figsize=(10, 6))
plt.plot(erro_index, erro_relativo, marker='o')
plt.title(f'Erro Relativo Percentual Ponto a Ponto\nNorma L2 Relativa: {norma_l2:.6f} (%)')
plt.xlabel("Index (z)")
plt.ylabel("Erro Relativo (%)")

plt.grid(True)
plt.savefig(OUT_DIR / "erro_relativo.png")

voltage_series = np.genfromtxt(OUT_DIR / "voltage_over_time.csv", delimiter=",", skip_header=1)

fig, ax = plt.subplots()
line, = ax.plot([], [], lw=2)
ax.set_xlim(1, voltage_series.shape[1] - 1)
ax.set_ylim(-1, 1)
ax.set_xlabel("Index (z)")
ax.set_ylabel("Normalized Voltage")


def init():
    line.set_data([], [])
    return line,


def update(frame):
    x = range(1, voltage_series.shape[1])
    y = voltage_series[frame, 1:]
    line.set_data(x, y)
    ax.set_title(f"Voltage at timestep {int(voltage_series[frame, 0])}")
    return line,


ani = FuncAnimation(fig, update, frames=len(voltage_series), init_func=init, blit=True)

ani.save(OUT_DIR / "voltage_simulation.mp4", writer="ffmpeg")
ani.save(OUT_DIR / "voltage_simulation.gif", writer="pillow")
plt.close(fig)

current_series = np.genfromtxt(OUT_DIR / "current_over_time.csv", delimiter=",", skip_header=1)

fig, ax = plt.subplots()
line, = ax.plot([], [], lw=2)
ax.set_xlim(1, current_series.shape[1] - 1)
ax.set_ylim(-1, 1)
ax.set_xlabel("Index (z)")
ax.set_ylabel("Current (A)")


def init():
    line.set_data([], [])
    return line,


def update(frame):
    x = range(1, current_series.shape[1])
    y = current_series[frame, 1:]
    line.set_data(x, y)
    ax.set_title(f"Current at timestep {int(current_series[frame, 0])}")
    return line,


ani = FuncAnimation(fig, update, frames=len(current_series), init_func=init, blit=True)

ani.save(OUT_DIR / "current_simulation.mp4", writer="ffmpeg")
ani.save(OUT_DIR / "current_simulation.gif", writer="pillow")

plt.close(fig)
