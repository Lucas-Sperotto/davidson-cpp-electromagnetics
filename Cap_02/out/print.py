import pandas as pd
import matplotlib.pyplot as plt

# Caminho do arquivo
file_path = 'erro_relativo.csv'

# Lê a tabela ignorando as duas primeiras linhas
df = pd.read_csv(file_path, skiprows=2)

# Captura a norma L2 da primeira linha
with open(file_path, 'r') as f:
    first_line = f.readline().strip()
    norma_l2 = float(first_line.split(':')[1].strip())

# Cria o gráfico
plt.figure(figsize=(10, 6))
plt.plot(df['Index'], df['Erro_relativo'], marker='o')
plt.title('Erro Relativo ponto a ponto')
plt.xlabel('Index (z)')
plt.ylabel('Erro Relativo (%)')

# Adiciona a norma L2 no topo
plt.text(0.5, 1.02, f'Norma L2 Relativa: {norma_l2:.6f}',
         ha='center', va='bottom',
         transform=plt.gca().transAxes,
         fontsize=12, fontweight='bold')

plt.grid(True)
plt.show()
