# Capítulo 2 - Simulações FDTD 1D

Este diretório contém a tradução para C++ dos códigos MATLAB do Capítulo 2 do livro:

> D.B. Davidson, *Computational Electromagnetics for RF and Microwave Engineering*, Cambridge University Press, 2ª edição.

## Estrutura

```
Cap_02/
├── src/
│   ├── fdtd_1D_demo.cpp             # Simulação FDTD 1D com senoide
│   ├── fdtd_1D_WB_demo.cpp          # Simulação FDTD 1D banda larga
│   └── gaussder_norm.hpp           # Função impulso gaussiano derivado normalizado
├── scripts/
│   ├── plot_fdtd_results.py        # Visualização de resultados da senoide
│   └── plot_fdtd_wb_results.py     # Visualização de resultados banda larga
├── out/                            # Saídas geradas pelos simuladores
└── README.md                       # Este arquivo
```

## Requisitos

- C++17 ou superior
- [FFTW3](http://www.fftw.org/) para cálculo da FFT
- Python 3 com `numpy` e `matplotlib` para visualização

## Compilação (exemplo)

```bash
sudo apt install libfftw3-dev
cd Cap_02
mkdir build && cd build
cmake ..
make
```

## Execução e Visualização

Execute os binários gerados dentro de `build/`:

```bash
./fdtd_1D_demo
./fdtd_1D_WB_demo
```

E visualize os resultados com:

```bash
cd ../scripts
python3 plot_fdtd_results.py
python3 plot_fdtd_wb_results.py
```

## Saídas geradas e suas interpretações

As imagens são salvas na pasta `Cap_02/out/`. Veja abaixo algumas delas:

### 📈 `comparison_voltage_plot.png`

> **Tensão no tempo final** da simulação senoidal (`fdtd_1D_demo`). Mostra a distribuição espacial da tensão após a convergência, juntamente com a solução analítica. Para verificar os parâmetros de execução, consulte o arquivo gerado: [Ver arquivo CSV](./out/simulation_parameters.csv).

<p align="center">
  <img src="/Cap_02/out/comparison_voltage_plot.png" alt="comparison_voltage_plot">
</p>

---

### 🌡️ `erro_relativo.png`

> **Erro relativo percentual** ao longo do espaço, correspondente aos valores do gráfico anterior.

<p align="center">
  <img src="/Cap_02/out/erro_relativo.png" alt="erro_relativo">
</p>

---

### ⚡ `voltage_simulation.gif`

> **Animação da evolução temporal da tensão** SV(z,t)S ao longo do espaço.

<p align="center">
  <img src="/Cap_02/out/voltage_simulation.gif" alt="voltage_simulation">
</p>

---

### 🔌 `current_simulation.gif`

> **Animação da evolução temporal da corrente** SI(z,t)S ao longo do espaço.

<p align="center">
  <img src="/Cap_02/out/current_simulation.gif" alt="current_simulation">
</p>

---

### 🕒 `wb_time_response.png`
> Gráfico temporal do pulso aplicado e da resposta na carga (versão wideband).

![wb_time_response](/Cap_02/out/wb_time_response.png)

---

### 📡 `wb_transfer_function.png`
> Função de transferência S|V_L / V_S|S em função da frequência, simulada via pulso e FFT.

![wb_transfer_function](/Cap_02/out/wb_transfer_function.png)

## Ligação com README geral

Este é o módulo correspondente ao **Capítulo 2** do livro. Consulte o [README principal](../README.md) para acessar os demais capítulos.
