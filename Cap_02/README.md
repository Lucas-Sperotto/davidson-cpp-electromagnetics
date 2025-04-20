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
cd Cap_02/src
g++ fdtd_1D_demo.cpp -o fdtd_demo -lfftw3 -lm
./fdtd_demo
```

## Visualização

```bash
cd ../scripts
python3 plot_fdtd_results.py
python3 plot_fdtd_wb_results.py
```

## Ligação com README geral

Este é o módulo correspondente ao **Capítulo 2** do livro. Consulte o [README principal](../README.md) para acessar os demais capítulos.
