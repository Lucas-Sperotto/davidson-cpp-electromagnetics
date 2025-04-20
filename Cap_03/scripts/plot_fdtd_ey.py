# plot_fdtd_ey.py
import matplotlib.pyplot as plt
import pandas as pd

# Carrega os dados
try:
    df_std = pd.read_csv('out/ey_point1.csv')
    plt.plot(df_std['tempo_ns'], df_std['Ey_Vpm'], label='Sem PML')
except FileNotFoundError:
    print("Arquivo ey_point1.csv não encontrado.")

try:
    df_pml = pd.read_csv('out/ey_point1_pml.csv')
    plt.plot(df_pml['tempo_ns'], df_pml['Ey_Vpm'], label='Com PML')
except FileNotFoundError:
    print("Arquivo ey_point1_pml.csv não encontrado.")

plt.xlabel('Tempo [ns]')
plt.ylabel('$E_y$ [V/m]')
plt.title('Campo $E_y$ no ponto de observação')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('out/ey_point1_comparison.png', dpi=300)
plt.show()
