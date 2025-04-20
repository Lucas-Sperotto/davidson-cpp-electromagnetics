# Cap_02/scripts/plot_fdtd_wb_results.py
import numpy as np
import matplotlib.pyplot as plt
import os

out_dir = "../out"

# Carrega sinais no tempo
vs = np.loadtxt(os.path.join(out_dir, "fdtd_wb_vs.csv"), delimiter=",")
vl = np.loadtxt(os.path.join(out_dir, "fdtd_wb_vl.csv"), delimiter=",")

# Plota sinais Vs(t) e Vl(t)
plt.figure()
plt.plot(vs[:,0], vs[:,1], label="Vs(t)")
plt.plot(vl[:,0], vl[:,1], label="Vl(t)", linestyle="--")
plt.title("Resposta no tempo: Fonte e Carga")
plt.xlabel("Tempo (s)")
plt.ylabel("Tensão (V)")
plt.legend()
plt.grid(True)
plt.savefig(os.path.join(out_dir, "wb_time_response.png"))

# Verifica se função de transferência foi salva
tf_path = os.path.join(out_dir, "fdtd_wb_tf.csv")
if os.path.exists(tf_path):
    tf = np.loadtxt(tf_path, delimiter=",")
    plt.figure()
    plt.plot(tf[:,0], tf[:,1])
    plt.title("Função de Transferência |Vl/Vs|")
    plt.xlabel("Frequência (Hz)")
    plt.ylabel("|Vl/Vs|")
    plt.grid(True)
    plt.savefig(os.path.join(out_dir, "wb_transfer_function.png"))

plt.show()
