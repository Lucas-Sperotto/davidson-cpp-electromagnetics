# plot_fdtd3d_fft.py
import matplotlib.pyplot as plt
import pandas as pd

# Carrega os dados da FFT
try:
    df_fft = pd.read_csv('out/hz_center_fft.csv')
    plt.plot(df_fft['freq_Hz'] / 1e6, df_fft['abs_Hz'], label='$|H_z(f)|$')
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
plt.savefig('out/hz_center_fft_plot.png', dpi=300)
plt.show()
