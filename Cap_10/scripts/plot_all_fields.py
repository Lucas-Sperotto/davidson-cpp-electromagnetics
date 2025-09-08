import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def plot_field_from_csv(filepath, save_dir):
    # Lê os dados
    df = pd.read_csv(filepath)

    x = df["x"].values
    y = df["y"].values
    Ex = df["Ex"].values
    Ey = df["Ey"].values

    x_unique = np.unique(x)
    y_unique = np.unique(y)

    X, Y = np.meshgrid(x_unique, y_unique)
    U = Ex.reshape(len(x_unique), len(y_unique)).T
    V = Ey.reshape(len(x_unique), len(y_unique)).T

    # Nome base
    name = os.path.splitext(os.path.basename(filepath))[0]

    # Cria figura
    plt.figure(figsize=(8, 6))
    plt.quiver(X, Y, U, V, scale=1, scale_units='xy', angles='xy')
    #plt.quiver(x, y, Ex, Ey, scale=1, scale_units='xy', angles='xy')
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(f"Campo vetorial: {name}")
    plt.axis("equal")
    plt.grid(True)
    plt.tight_layout()

    # Salva como imagem
    out_file = os.path.join(save_dir, name + ".png")
    plt.savefig(out_file)
    plt.close()
    print(f"Imagem salva: {out_file}")


if __name__ == "__main__":
    data_dir = "../out"
    save_dir = "../out/img"

    os.makedirs(save_dir, exist_ok=True)

    for file in os.listdir(data_dir):
        if file.endswith(".csv") and file.startswith("field_kc_"):
            filepath = os.path.join(data_dir, file)
            plot_field_from_csv(filepath, save_dir)
