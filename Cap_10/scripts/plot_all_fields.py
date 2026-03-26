from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


BASE_DIR = Path(__file__).resolve().parents[1]
OUT_DIR = BASE_DIR / "out"
IMG_DIR = OUT_DIR / "img"


def load_field(filepath: Path):
    data = np.genfromtxt(filepath, delimiter=",", names=True)
    x = data["x"]
    y = data["y"]
    ex = data["Ex"]
    ey = data["Ey"]

    x_unique = np.unique(x)
    y_unique = np.unique(y)
    nx = len(x_unique)
    ny = len(y_unique)

    x_grid = x.reshape(nx, ny).T
    y_grid = y.reshape(nx, ny).T
    ex_grid = ex.reshape(nx, ny).T
    ey_grid = ey.reshape(nx, ny).T
    return x_grid, y_grid, ex_grid, ey_grid


def render_field(filepath: Path, output: Path):
    x_grid, y_grid, ex_grid, ey_grid = load_field(filepath)

    plt.figure(figsize=(8, 6))
    plt.quiver(x_grid, y_grid, ex_grid, ey_grid, scale=1, scale_units="xy", angles="xy")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(filepath.stem)
    plt.axis("equal")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output, dpi=200)
    plt.close()
    print(f"Imagem salva: {output}")


if __name__ == "__main__":
    IMG_DIR.mkdir(parents=True, exist_ok=True)

    patterns = ["field_kc_*.csv", "spur_field*_kc_*.csv"]
    files = []
    for pattern in patterns:
        files.extend(sorted(OUT_DIR.glob(pattern)))

    for filepath in files:
        render_field(filepath, IMG_DIR / f"{filepath.stem}.png")
