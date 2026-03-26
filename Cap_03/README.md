# Capítulo 3 - FDTD 2D e 3D

Este diretório reúne a tradução para C++ dos demos FDTD do Capítulo 3 do livro de Davidson, com foco em:

- espalhamento 2D com fronteira campo total/campo espalhado;
- absorção por PML em 2D;
- análise modal de cavidade 3D via FDTD e FFT.

## Status

Estado atual do capítulo:

- `fdtd_2d_demo.m` -> traduzido em `fdtd_2d_demo.cpp`;
- `fdtd_2d_pml_demo.m` -> traduzido em `fdtd2d_pml.cpp` com adaptacao controlada na injecao scat/tot do PML;
- `fdtd_3D_demo.m` -> traduzido em `fdtd3d.cpp`;
- `gaussder_norm.m` -> incorporado em `gaussder.cpp`;
- `PMLperformance.m` -> traduzido em `scripts/pml_performance.py`;
- `cyl_test.m` -> traduzido em `scripts/cyl_test.py`.

Os demos principais já compilam e rodam, e os scripts auxiliares principais de comparação já possuem equivalente. O capítulo agora entra numa fase de validação fina e refinamento didático.

## Estrutura Atual

```text
Cap_03/
├── include/
│   ├── fdtd2d.hpp
│   ├── fdtd2d_pml.hpp
│   └── fdtd3d.hpp
├── src/
│   ├── fdtd_2d_demo.cpp
│   ├── fdtd2d_pml.cpp
│   ├── fdtd3d.cpp
│   └── gaussder.cpp
├── scripts/
│   ├── plot_fdtd_ey.py
│   ├── plot_fdtd3d_fft.py
│   ├── pml_performance.py
│   └── cyl_test.py
├── out/
│   ├── ey_point1.csv
│   ├── ey_point1_meta.csv
│   ├── ey_point1_pml.csv
│   ├── ey_point1_pml_meta.csv
│   ├── hz_center_fft.csv
│   ├── hz_center_time.csv
│   ├── hz_center_meta.csv
│   ├── ey_point1_comparison.png
│   └── hz_center_fft_plot.png
├── main_fdtd2d.cpp
├── main_fdtd2d_pml.cpp
├── main_fdtd3d.cpp
├── CMakeLists.txt
└── README.md
```

## Compilação

```bash
cd Cap_03
cmake -S . -B build
cmake --build build -j$(nproc)
```

## Executáveis

- `./build/fdtd_2d_demo` - demo 2D com fronteiras ABC.
- `./build/fdtd2d_pml` - demo 2D com PML.
- `./build/fdtd3d` - demo 3D para cavidade retangular PEC.

Os executáveis agora aceitam parâmetros por linha de comando.

Exemplos:

```bash
./build/fdtd_2d_demo --cyl-present --refine 1 --pulse-compress 2
./build/fdtd2d_pml --refine 1 --d-cell 10 --poly-m 3
./build/fdtd3d --refine 2 --seed 12345
```

## Scripts Originais de Referência

Arquivos MATLAB usados como base local:

- `original_matlab/Chapter 3/FDTD_2D/fdtd_2d_demo.m`
- `original_matlab/Chapter 3/FDTD_2D/fdtd_2d_pml_demo.m`
- `original_matlab/Chapter 3/FDTD_2D/PMLperformance.m`
- `original_matlab/Chapter 3/FDTD_2D/cyl_test.m`
- `original_matlab/Chapter 3/FDTD_3D/fdtd_3D_demo.m`

## Saídas

Arquivos gerados em `out/`:

- `ey_point1.csv` - histórico temporal do campo elétrico para o caso 2D com ABC;
- `ey_point1_meta.csv` - metadados da rodada 2D com ABC;
- `ey_point1_pml.csv` - histórico temporal do campo elétrico para o caso 2D com PML;
- `ey_point1_pml_meta.csv` - metadados da rodada 2D com PML;
- `hz_center_fft.csv` - espectro do campo magnético no centro da cavidade 3D;
- `hz_center_time.csv` - histórico temporal do campo magnético no centro da cavidade;
- `hz_center_meta.csv` - metadados da rodada 3D;
- `ey_point1_comparison.png` - comparação entre ABC e PML;
- `hz_center_fft_plot.png` - gráfico do espectro 3D.

## Visualização

```bash
cd Cap_03
python3 scripts/plot_fdtd_ey.py
python3 scripts/plot_fdtd3d_fft.py
python3 scripts/pml_performance.py --thick out/ey_point1_pml_thick.csv --thick-meta out/ey_point1_pml_thick_meta.csv --reference out/ey_point1.csv
python3 scripts/cyl_test.py --old out/cyl_old.csv --new out/cyl_new.csv --new-longer out/cyl_new_longer.csv
```

## Limitações Atuais

- O demo 2D com ABC ainda contém trechos explicitamente marcados no código como não testados ou sujeitos a sinal espúrio nos cantos da interface campo total/campo espalhado.
- O demo PML atual preserva a estrutura, os parâmetros e os pós-processamentos do original, mas a injeção scat/tot nas componentes split de `H_z` segue marcada para validação adicional; por enquanto, a tradução mantém a formulação estável já validada no port C++.
- Os scripts auxiliares já foram portados, mas ainda precisam ser usados com um conjunto de rodadas de referência mais representativo.

## Próximos Passos Naturais

- consolidar a validação numérica de `fdtd_2d_demo.cpp` e `fdtd2d_pml.cpp`;
- documentar a teoria do capítulo com o mesmo nível didático do `Cap_02`;
- validar numericamente os resultados C++ contra os dados MATLAB de referência.
