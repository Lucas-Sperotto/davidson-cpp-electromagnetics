import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def plot_field_from_csv(filepath):
    # Lê os dados
    df = pd.read_csv(filepath)
    
    # Extrai coordenadas e componentes do campo
    x = df["x"].values
    y = df["y"].values
    Ex = df["Ex"].values
    Ey = df["Ey"].values

    # Determina a grade
    x_unique = np.unique(x)
    y_unique = np.unique(y)

    X, Y = np.meshgrid(x_unique, y_unique)
    U = Ex.reshape(len(x_unique), len(y_unique)).T
    V = Ey.reshape(len(x_unique), len(y_unique)).T

    # Cria o gráfico
    plt.figure(figsize=(8, 6))
    plt.quiver(X, Y, U, V, scale=1, scale_units='xy', angles='xy')
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(f"Campo vetorial: {os.path.basename(filepath)}")
    plt.axis("equal")
    plt.grid(True)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Visualiza campos vetoriais salvos em CSV.")
    parser.add_argument("filename", type=str, help="Nome do arquivo CSV na pasta ../out/")
    args = parser.parse_args()

    filepath = os.path.join("../out", args.filename)
    
    if not os.path.exists(filepath):
        print(f"Arquivo não encontrado: {filepath}")
    else:
        plot_field_from_csv(filepath)
