# Capítulo 9 - FEM 1D para Linha de Transmissão

Este diretório reúne a tradução para C++ dos códigos MATLAB do Capítulo 9 do livro de Davidson para a solução 1D da equação de onda em uma linha de transmissão usando elementos finitos lineares de primeira ordem.

## Status

Estado atual do capítulo:

- `FEM_1D_1st.m` -> traduzido em `FEM_1D_1st.cpp`;
- `FEM_1D_solver.m` -> traduzido em `FEM_1D_solver.cpp/.h`;
- `FEM_pp.m` -> traduzido em `FEM_pp.cpp/.h`.

O núcleo numérico está fiel ao fluxo MATLAB original. As principais adaptações desta versão em C++ são:

- interface por linha de comando, com suporte tanto a argumentos nomeados quanto ao modo posicional legado;
- gravação padronizada de resultados em `Cap_09/out/`;
- geração explícita de CSVs de perfis, convergência, nós e metadados de execução.

## Teoria

O problema modelado é uma linha de transmissão uniforme, com parâmetros distribuídos `L` e `C`, excitada em uma extremidade e terminada em circuito aberto na outra. O capítulo usa a analogia entre a linha TEM e a equação de onda 1D para mostrar como o método dos elementos finitos pode ser aplicado em um caso simples e bem controlado.

Assumindo:

- indutância por unidade de comprimento `L`;
- capacitância por unidade de comprimento `C`;
- frequência `f`;
- velocidade de propagação `c = 1/sqrt(LC)`;
- comprimento da linha `ell = lambda/2`, com `lambda = c/f`.

O capítulo discretiza a linha em `N_elem` elementos lineares e monta o sistema matricial global a partir de duas contribuições:

- uma matriz de rigidez associada ao termo espacial;
- uma matriz de massa associada ao termo harmônico em frequência.

Depois disso:

1. resolve-se o sistema reduzido com a condição prescrita em uma extremidade;
2. reconstrói-se a solução nodal completa;
3. interpola-se a solução FEM em pontos internos para comparação com a expressão analítica usada no script MATLAB;
4. calcula-se o erro RMS e, com refinamentos sucessivos, a curva de convergência.

## Estrutura

```text
Cap_09/
├── CMakeLists.txt
├── FEM_1D_1st.cpp
├── FEM_1D_solver.cpp
├── FEM_1D_solver.h
├── FEM_pp.cpp
├── FEM_pp.h
├── plot_fem.py
├── out/
│   ├── fem_profile_stage_#.csv
│   ├── fem_profile_stage_#.png
│   ├── fem_nodes_stage_#.csv
│   ├── fem_convergence.csv
│   ├── fem_run_metadata.csv
│   └── convergence.png
└── README.md
```

## Arquivos MATLAB de Referência

Arquivos usados como base local para a tradução:

- `original_matlab/Chapter 9/FEM_1D_1st.m`
- `original_matlab/Chapter 9/FEM_1D_solver.m`
- `original_matlab/Chapter 9/FEM_pp.m`

## Compilação

```bash
cd Cap_09
cmake -S . -B build
cmake --build build -j$(nproc)
```

Executável gerado:

```bash
./build/fem_tl
```

## Execução

Forma recomendada, com argumentos nomeados:

```bash
./build/fem_tl --n-elem-init 2 --num-meshes 4 --freq 1.0 --inductance 1.0 --capacitance 1.0 --vin 1.0
```

Forma legada, mantida por compatibilidade:

```bash
./build/fem_tl 2 4 1.0 1.0 1.0 1.0
```

Parâmetros:

- `--n-elem-init` - número inicial de elementos, padrão `2`;
- `--num-meshes` - número de etapas de refinamento, padrão `1`;
- `--freq` - frequência em Hz, padrão `1.0`;
- `--inductance` - indutância por unidade de comprimento, padrão `1.0`;
- `--capacitance` - capacitância por unidade de comprimento, padrão `1.0`;
- `--vin` - tensão prescrita no último nó, padrão `1.0`.

Ajuda rápida:

```bash
./build/fem_tl --help
```

## Saídas

Os arquivos são gerados em `Cap_09/out/`:

- `fem_profile_stage_#.csv` - perfil interpolado da solução FEM e solução analítica em cada estágio;
- `fem_nodes_stage_#.csv` - valores nodais da solução FEM em cada estágio;
- `fem_convergence.csv` - curva `h/lambda` versus erro RMS;
- `fem_run_metadata.csv` - parâmetros da rodada e inclinação log-log quando houver mais de uma malha.

## Visualização

```bash
cd Cap_09
python3 plot_fem.py
```

Arquivos gerados:

- `out/convergence.png` - curva log-log de convergência;
- `out/fem_profile_stage_#.png` - comparação FEM versus solução analítica em cada estágio.

## Interpretação dos Resultados

- `fem_profile_stage_#.png` mostra como a solução interpolada por elementos finitos acompanha a solução analítica usada no exemplo do livro.
- `fem_convergence.csv` e `convergence.png` mostram a queda do erro RMS conforme a malha é refinada.
- `fem_run_metadata.csv` ajuda a documentar a rodada e a guardar a inclinação observada da curva de convergência.

## Limitações e Observações

- A expressão analítica de comparação foi preservada do script MATLAB original, inclusive na forma como é escrita no exemplo.
- Este capítulo ainda pode crescer em documentação teórica, mas o fluxo numérico principal já está bem estabilizado.
- Como o foco aqui é fidelidade ao material original, as matrizes ainda são montadas de forma explícita, em vez de usar uma formulação mais enxuta ou esparsa.

## Ligação com o Repositório

Este é o módulo correspondente ao **Capítulo 9**. Consulte o [README principal](../README.md) e a [matriz de tradução](../TRANSLATION_MATRIX.md) para acompanhar o restante do projeto.
