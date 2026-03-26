from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

OUT_DIR = Path(__file__).resolve().parent.parent / "out"


def load_fft(path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.genfromtxt(path, delimiter=",", names=True, dtype=float, encoding="utf-8")
    return np.atleast_1d(data["freq_Hz"]), np.atleast_1d(data["abs_Hz"])

# Carrega os dados da FFT
try:
    freq_hz, abs_hz = load_fft(OUT_DIR / 'hz_center_fft.csv')
    plt.plot(freq_hz / 1e6, abs_hz, label='$|H_z(f)|$')
except FileNotFoundError:
    print("Arquivo hz_center_fft.csv não encontrado.")
    exit()

# Linhas verticais para frequências analíticas (modo TE101 e outros)
f_exact = [249662404.0, 360048742.0, 390255684.0, 426286054.0]  # Hz
for f in f_exact:
    plt.axvline(x=f/1e6, color='k', linestyle='--')

plt.xlabel('Frequência [MHz]')
plt.ylabel('$|H_z|$ (normalizado)')
plt.title('Espectro de Fourier de $H_z$ no centro da cavidade')
plt.grid(True)
plt.tight_layout()
plt.savefig(OUT_DIR / 'hz_center_fft_plot.png', dpi=300)
plt.close()
