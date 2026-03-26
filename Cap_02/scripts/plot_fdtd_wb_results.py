import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from pathlib import Path

OUT_DIR = Path(__file__).resolve().parents[1] / "out"

# Carrega sinais no tempo
vs = np.loadtxt(OUT_DIR / "fdtd_wb_vs.csv", delimiter=",")
vl = np.loadtxt(OUT_DIR / "fdtd_wb_vl.csv", delimiter=",")

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
plt.savefig(OUT_DIR / "wb_time_response.png")

# Verifica se função de transferência foi salva
tf_path = OUT_DIR / "fdtd_wb_tf.csv"
if tf_path.exists():
    tf = np.loadtxt(tf_path, delimiter=",")
    plt.figure()
    plt.plot(tf[:,0], tf[:,1])
    plt.title(r"Função de Transferência $|V_L/V_O|$")
    plt.xlabel("Frequência (Hz)")
    plt.ylabel(r"$|V_L/V_O|$")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(OUT_DIR / "wb_transfer_function.png")

plt.close("all")

voltage_series = np.genfromtxt(OUT_DIR / "WB_voltage_over_time.csv", delimiter=",", skip_header=1)

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

ani.save(OUT_DIR / "WB_voltage_simulation.mp4", writer="ffmpeg")
ani.save(OUT_DIR / "WB_voltage_simulation.gif", writer="pillow")
plt.close(fig)
