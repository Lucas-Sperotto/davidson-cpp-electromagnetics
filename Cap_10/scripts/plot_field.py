from pathlib import Path
import argparse

import matplotlib.pyplot as plt
import numpy as np


OUT_DIR = Path(__file__).resolve().parents[1] / "out"


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


def plot_field_from_csv(filepath: Path, output: Path | None):
    x_grid, y_grid, ex_grid, ey_grid = load_field(filepath)

    plt.figure(figsize=(8, 6))
    plt.quiver(x_grid, y_grid, ex_grid, ey_grid, scale=1, scale_units="xy", angles="xy")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(f"Campo vetorial: {filepath.name}")
    plt.axis("equal")
    plt.grid(True)
    plt.tight_layout()

    if output is not None:
        plt.savefig(output, dpi=200)
        plt.close()
        print(f"Imagem salva em {output}")
    else:
        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Visualiza campos vetoriais salvos em CSV.")
    parser.add_argument("filename", type=str, help="Nome do arquivo CSV na pasta out/")
    parser.add_argument("--save", type=str, default="", help="Salva a figura em vez de abrir janela.")
    args = parser.parse_args()

    filepath = OUT_DIR / args.filename
    if not filepath.exists():
        raise SystemExit(f"Arquivo nao encontrado: {filepath}")

    output = Path(args.save) if args.save else None
    plot_field_from_csv(filepath, output)
