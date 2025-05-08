# Cap_02/scripts/plot_fdtd_wb_results.py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os

out_dir = "../out"

# Carrega sinais no tempo
vs = np.loadtxt(os.path.join(out_dir, "fdtd_wb_vs.csv"), delimiter=",")
vl = np.loadtxt(os.path.join(out_dir, "fdtd_wb_vl.csv"), delimiter=",")

# Plota sinais Vs(t) e Vl(t)
plt.figure()
plt.plot(vs[:,0], vs[:,1], label=r"$V_+(t)$")
plt.plot(vl[:,0], vl[:,1], label=r"$V_L(t)$", linestyle="--")
plt.title("Tensão da Onda Incidente e Tensão na Carga.")
plt.xlabel("Tempo (s)")
plt.ylabel("Tensão (V)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "wb_time_response.png"))

# Verifica se função de transferência foi salva
tf_path = os.path.join(out_dir, "fdtd_wb_tf.csv")
if os.path.exists(tf_path):
    tf = np.loadtxt(tf_path, delimiter=",")
    plt.figure()
    plt.plot(tf[:,0], tf[:,1])
    plt.title(r"Função de Transferência $|V_L/V_O|$")
    plt.xlabel("Frequência (Hz)")
    plt.ylabel(r"$|V_L/V_O|$")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "wb_transfer_function.png"))

plt.show()

# Lê o arquivo de tensões
voltage_series = pd.read_csv(os.path.join(out_dir, 'WB_voltage_over_time.csv'))

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
ani.save(os.path.join(out_dir, 'WB_voltage_simulation.mp4'), writer='ffmpeg')
ani.save(os.path.join(out_dir, 'WB_voltage_simulation.gif'), writer='pillow')
plt.close(fig)  # fecha apenas depois de salvar
