# Capítulo 11 - FEM Vetorial 3D

Este diretório reúne a tradução para C++ dos códigos MATLAB do Capítulo 11 do livro de Davidson. O capítulo estuda modos próprios 3D em cavidades PEC usando tetraedros, elementos de Whitney/CTLN, elementos LTQN e alguns utilitários de malha e quadratura associados.

## Status

Estado atual do capítulo:

- `Eigen3D_v0.m` -> traduzido em `Eigen3D_v0.cpp`;
- `Eigen3D_CTLN.m` -> traduzido em `Eigen3D_CTLN.cpp`;
- `Eigen3D_LTQN.m` -> traduzido em `Eigen3D_LTQN.cpp`;
- `FETD_FDTD.m` -> traduzido em `FETD_FDTD.cpp`;
- `test_tet_quad.m` -> traduzido em `test_tet_quad.cpp`;
- `brick_mesh.m`, `read_gmsh2.m`, `edgemake3D.m`, `facemake3D.m`, `free_dof3D.m`, `free_dof3D_LTQN.m`, `free_nodes3D.m`, `renumber_dof.m`, `renumber_dof_LTQN.m`, `sandt3D.m`, `sandt3D_LTQN.m`, `LTQN3D.m`, `curl_LTQN3D.m`, `tet_quad.m`, `avg_mesh_length.m`, `eig_err3D.m`, `TE_cavity_modes3D.m` e `cavity_modes3D.m` -> traduzidos em `src/`.

O núcleo numérico do capítulo já está funcional. As principais adaptações desta versão em C++ são:

- substituição da geração via `delaunay3` por uma decomposição interna determinística de cada brick em seis tetraedros;
- interface por linha de comando para os scripts principais;
- gravação padronizada de CSVs em `Cap_11/out/`;
- geração de figuras por script Python, em vez de `plot(...)` e `tetramesh(...)` do MATLAB.

## Teoria

O capítulo estende a formulação vetorial do Capítulo 10 para 3D. A cavidade retangular PEC é discretizada com tetraedros, e o problema de autovalores generalizado é montado a partir das matrizes:

- `S`, associada ao termo de rotacional;
- `T`, associada ao termo de massa.

No caso `CTLN`, os graus de liberdade vivem nas arestas e correspondem às funções de Whitney de primeira ordem. No caso `LTQN`, entram também funções adicionais de aresta e de face, elevando a ordem do espaço aproximante.

As rotinas auxiliares do capítulo fazem exatamente o que seus nomes sugerem no MATLAB:

- `brick_mesh` gera os nós de uma malha regular hexaédrica;
- `edgemake3D` e `facemake3D` constroem topologia global;
- `free_dof3D`, `free_dof3D_LTQN` e `free_nodes3D` aplicam as condições de contorno PEC/simetria;
- `sandt3D` e `sandt3D_LTQN` calculam as matrizes elementares;
- `tet_quad`, `LTQN3D` e `curl_LTQN3D` apoiam a formulação de segunda ordem;
- `eig_err3D` compara os primeiros autovalores com o caso de referência da cavidade retangular.

## Estrutura

```text
Cap_11/
├── CMakeLists.txt
├── Eigen3D_v0.cpp
├── Eigen3D_CTLN.cpp
├── Eigen3D_LTQN.cpp
├── FETD_FDTD.cpp
├── test_tet_quad.cpp
├── include/
│   ├── *.hpp
├── src/
│   ├── *.cpp
├── scripts/
│   └── plot_cap11.py
├── out/
│   ├── eigdata3D_*.csv
│   ├── eigen3d_*_summary.csv
│   ├── fetd_fdtd_*.csv
│   ├── test_tet_quad_summary.csv
│   └── *.png
└── README.md
```

## Arquivos MATLAB de Referência

Arquivos usados como base local para a tradução:

- `original_matlab/Chapter 11/Eigen3D_v0.m`
- `original_matlab/Chapter 11/Eigen3D_CTLN.m`
- `original_matlab/Chapter 11/Eigen3D_LTQN.m`
- `original_matlab/Chapter 11/FETD_FDTD.m`
- `original_matlab/Chapter 11/LTQN3D.m`
- `original_matlab/Chapter 11/TE_cavity_modes3D.m`
- `original_matlab/Chapter 11/avg_mesh_length.m`
- `original_matlab/Chapter 11/brick_mesh.m`
- `original_matlab/Chapter 11/cavity_modes3D.m`
- `original_matlab/Chapter 11/curl_LTQN3D.m`
- `original_matlab/Chapter 11/edgemake3D.m`
- `original_matlab/Chapter 11/eig_err3D.m`
- `original_matlab/Chapter 11/facemake3D.m`
- `original_matlab/Chapter 11/free_dof3D.m`
- `original_matlab/Chapter 11/free_dof3D_LTQN.m`
- `original_matlab/Chapter 11/free_nodes3D.m`
- `original_matlab/Chapter 11/read_gmsh2.m`
- `original_matlab/Chapter 11/renumber_dof.m`
- `original_matlab/Chapter 11/renumber_dof_LTQN.m`
- `original_matlab/Chapter 11/sandt3D.m`
- `original_matlab/Chapter 11/sandt3D_LTQN.m`
- `original_matlab/Chapter 11/test_tet_quad.m`
- `original_matlab/Chapter 11/tet_quad.m`

## Compilação

```bash
cd Cap_11
cmake -S . -B build
cmake --build build -j$(nproc)
```

Executáveis gerados:

```bash
./build/Eigen3D_v0
./build/Eigen3D_CTLN
./build/Eigen3D_LTQN
./build/FETD_FDTD
./build/test_tet_quad
```

## Execução

Versão inicial com malha interna:

```bash
./build/Eigen3D_v0 --nx 2 --ny 1 --nz 2
```

CTLN com malha interna:

```bash
./build/Eigen3D_CTLN
```

CTLN lendo malha externa Gmsh:

```bash
./build/Eigen3D_CTLN --internal-mesh 0 --mesh-file "../original_matlab/Chapter 11/Gmsh files/box_20.msh"
```

LTQN:

```bash
./build/Eigen3D_LTQN
./build/Eigen3D_LTQN --element-order 2 --nx 2 --ny 1 --nz 2
```

Funções CTLN do exemplo FETD-FDTD:

```bash
./build/FETD_FDTD
```

Teste da quadratura tetraédrica:

```bash
./build/test_tet_quad
```

## Saídas

Os arquivos são gerados em `Cap_11/out/`:

- `eigdata3D_v0_modes.csv` e `eigen3d_v0_summary.csv`;
- `eigdata3D_ctlN_modes.csv` e `eigen3d_ctln_summary.csv`;
- `eigdata3D_ltqn_modes.csv` e `eigen3d_ltqn_summary.csv`;
- `fetd_fdtd_ctlN_funcs.csv` e `fetd_fdtd_summary.csv`;
- `test_tet_quad_summary.csv`.

## Visualização

```bash
cd Cap_11
python3 scripts/plot_cap11.py
```

Figuras geradas:

- `out/fetd_fdtd_ctlN_funcs.png`
- `out/test_tet_quad_errors.png`
- `out/eigen3d_v0_modes.png`
- `out/eigen3d_v0_modes_relerr.png`
- `out/eigen3d_ctln_modes.png`
- `out/eigen3d_ctln_modes_relerr.png`
- `out/eigen3d_ltqn_modes.png`
- `out/eigen3d_ltqn_modes_relerr.png`

## Resultados de Referência

Na validação inicial desta tradução:

- `Eigen3D_CTLN` com malha interna padrão (`2 x 1 x 2` bricks) gerou `24` tetraedros e `9` DOFs livres;
- `Eigen3D_LTQN --element-order 2 --nx 2 --ny 1 --nz 2` gerou `82` DOFs e apresentou os primeiros modos físicos próximos de `5.2495`, `7.0475`, `7.4577`, `7.5660`;
- `test_tet_quad` confirmou a regra de 11 pontos com erros relativos numéricos da ordem de máquina;
- `Eigen3D_CTLN --internal-mesh 0 --mesh-file "../original_matlab/Chapter 11/Gmsh files/box_20.msh"` rodou com sucesso, validando também a leitura do caminho Gmsh.

## Limitações e Observações

- A geração interna da malha foi adaptada: em vez de `delaunay3`, cada brick é dividido deterministicamente em seis tetraedros. Isso mantém a ideia do capítulo, mas não reproduz exatamente a conectividade do MATLAB em todos os casos.
- O arquivo `cavity_modes3D.m` original já está marcado como incompleto, e a tradução preserva esse estado.
- O `Eigen3D_v0` mantém os defaults originais de malha mais densa, mas a validação prática aqui foi feita em malha reduzida para manter o custo razoável.

## Ligação com o Repositório

Este é o módulo correspondente ao **Capítulo 11**. Consulte o [README principal](../README.md) e a [matriz de tradução](../TRANSLATION_MATRIX.md) para acompanhar o restante do projeto.
