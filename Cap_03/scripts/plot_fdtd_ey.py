from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

OUT_DIR = Path(__file__).resolve().parent.parent / "out"


def load_signal(path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.genfromtxt(path, delimiter=",", names=True, dtype=float, encoding="utf-8")
    return np.atleast_1d(data["tempo_ns"]), np.atleast_1d(data["Ey_Vpm"])


# Carrega os dados
try:
    tempo_ns, ey = load_signal(OUT_DIR / "ey_point1.csv")
    plt.plot(tempo_ns, ey, label='Sem PML')
except FileNotFoundError:
    print("Arquivo ey_point1.csv não encontrado.")

try:
    tempo_ns, ey = load_signal(OUT_DIR / "ey_point1_pml.csv")
    plt.plot(tempo_ns, ey, label='Com PML')
except FileNotFoundError:
    print("Arquivo ey_point1_pml.csv não encontrado.")

plt.xlabel('Tempo [ns]')
plt.ylabel('$E_y$ [V/m]')
plt.title('Campo $E_y$ no ponto de observação')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(OUT_DIR / 'ey_point1_comparison.png', dpi=300)
plt.close()
