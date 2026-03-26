from pathlib import Path
import csv

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np


OUT_DIR = Path(__file__).resolve().parents[1] / "out"


def read_single_row_csv(path: Path):
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise RuntimeError(f"Arquivo vazio: {path}")
    return rows[0]


field = np.genfromtxt(OUT_DIR / "field_map.csv", delimiter=",", names=True)
geom = read_single_row_csv(OUT_DIR / "geom_params.csv")
strip = np.genfromtxt(OUT_DIR / "strip_nodes.csv", delimiter=",", names=True)

a_half = float(geom["a_half"])
b = float(geom["b"])
h = float(geom["h"])
nx = int(float(geom["Nx"])) + 1
ny = int(float(geom["Ny"])) + 1

x = field["x"]
y = field["y"]
phi = field["phi"]
ex = field["Ex"]
ey = field["Ey"]

x_grid = x.reshape(nx, ny).T
y_grid = y.reshape(nx, ny).T
phi_grid = phi.reshape(nx, ny).T
ex_grid = ex.reshape(nx, ny).T
ey_grid = ey.reshape(nx, ny).T

plt.figure(figsize=(9, 4.5))

cont = plt.contour(x_grid, y_grid, phi_grid, levels=20, linewidths=0.9)
plt.clabel(cont, inline=True, fontsize=8, fmt="%.2f")

step_x = max(1, (nx - 1) // 25)
step_y = max(1, (ny - 1) // 20)
plt.quiver(
    x_grid[::step_y, ::step_x],
    y_grid[::step_y, ::step_x],
    ex_grid[::step_y, ::step_x],
    ey_grid[::step_y, ::step_x],
    angles="xy",
    scale_units="xy",
    scale=None,
    width=0.002,
)

plt.axhline(h, linestyle="--", linewidth=1.5, label="interface dielétrica (y=h)")

if strip.size > 0:
    strip_x = np.atleast_1d(strip["x"])
    strip_y = np.atleast_1d(strip["y"])
    minx, maxx = np.min(strip_x), np.max(strip_x)
    miny, maxy = np.min(strip_y), np.max(strip_y)
    rect = Rectangle((minx, miny), maxx - minx, maxy - miny,
                     linewidth=2, fill=False, color="C3", label="microstrip (V)")
    plt.gca().add_patch(rect)
    plt.scatter(strip_x, strip_y, s=6, color="C3")

plt.xlim(0, a_half)
plt.ylim(0, b)
plt.gca().set_aspect("equal", "box")
plt.xlabel("x (unid.)")
plt.ylabel("y (unid.)")
plt.title("Microstrip: equipotenciais e campo elétrico (E)")
plt.legend(loc="upper right", fontsize=8)
plt.tight_layout()
plt.savefig(OUT_DIR / "microstrip_field.png", dpi=200)
plt.close()

print(f"Figura salva em {OUT_DIR / 'microstrip_field.png'}")
