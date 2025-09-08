import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# === Lê dados ===
df = pd.read_csv("../out/field_map.csv")        # x,y,phi,Ex,Ey
gp = pd.read_csv("../out/geom_params.csv")      # a_half,b,h,w,Nx,Ny
gs = pd.read_csv("../out/strip_nodes.csv")      # x,y

a_half = float(gp["a_half"].iloc[0])
b      = float(gp["b"].iloc[0])
h      = float(gp["h"].iloc[0])
w      = float(gp["w"].iloc[0])
Nx     = int(gp["Nx"].iloc[0])
Ny     = int(gp["Ny"].iloc[0])

nx1 = Nx + 1
ny1 = Ny + 1

# Constrói grades
# Atenção: df está em ordem varrendo i mais rápido; pivot facilita obter matrizes
phi_mat = df.pivot(index="y", columns="x", values="phi").values
Ex_mat  = df.pivot(index="y", columns="x", values="Ex").values
Ey_mat  = df.pivot(index="y", columns="x", values="Ey").values

# Eixos X,Y (ordenados crescentes)
xs = np.sort(df["x"].unique())
ys = np.sort(df["y"].unique())
X, Y = np.meshgrid(xs, ys)

# === Figura ===
plt.figure(figsize=(9, 4.5))

# 1) Equipotenciais (linhas)
levels = 20  # ajuste se quiser mais/menos linhas
cont = plt.contour(X, Y, phi_mat, levels=levels, linewidths=0.9)
plt.clabel(cont, inline=True, fontsize=8, fmt="%.2f")

# 2) Campo elétrico (setas) - faça uma amostragem para não poluir
step_x = max(1, Nx // 25)  # ajuste densidade
step_y = max(1, Ny // 20)
plt.quiver(
    X[::step_y, ::step_x],
    Y[::step_y, ::step_x],
    Ex_mat[::step_y, ::step_x],
    Ey_mat[::step_y, ::step_x],
    angles='xy', scale_units='xy', scale=None, width=0.002
)

# 3) Interface dielétrica y = h
plt.axhline(h, linestyle="--", linewidth=1.5, label="interface dielétrica (y=h)")

# 4) Microstrip (a partir dos nós prescritos): desenha o retângulo mínimo que os contém
if len(gs) > 0:
    minx, maxx = gs["x"].min(), gs["x"].max()
    miny, maxy = gs["y"].min(), gs["y"].max()
    rect = Rectangle((minx, miny), maxx - minx, maxy - miny,
                     linewidth=2, fill=False, color="C3", label="microstrip (V)")
    plt.gca().add_patch(rect)
    # (Opcional) pontos dos nós prescritos:
    plt.scatter(gs["x"], gs["y"], s=6, color="C3")

# Moldura e rótulos
plt.xlim(0, a_half)
plt.ylim(0, b)
plt.gca().set_aspect('equal', 'box')
plt.xlabel("x (unid.)")
plt.ylabel("y (unid.)")
plt.title("Microstrip: equipotenciais e campo elétrico (E)")

plt.legend(loc="upper right", fontsize=8)
plt.tight_layout()
plt.savefig("../out/microstrip_field.png", dpi=200)
plt.show()
