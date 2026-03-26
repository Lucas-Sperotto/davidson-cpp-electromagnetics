# Capítulo 6 - MoM 3D com RWG e RCS de Esfera

Este diretório reúne a tradução para C++ dos códigos MATLAB do Capítulo 6 do livro de Davidson. O capítulo combina dois blocos:

- a solução analítica da RCS monostática de uma esfera PEC;
- a formulação de espalhamento 3D por EFIE de potencial misto com funções de base RWG sobre uma placa PEC triangularizada.

## Status

Estado atual do capítulo:

- `MoM3D_demo.m` -> traduzido em `MoM3D_demo.cpp`;
- `converge_sphereRCS.m` -> traduzido em `converge_sphereRCS.cpp`;
- `sphereRCS.m` -> traduzido em `src/sphereRCS.cpp`;
- `ComputeRho_c.m`, `FillVVector.m`, `FillZMatrixByEdge.m`, `FillZMatrixByFace.m`, `GL_quad_rule.m`, `Int_pq.m`, `PostProcMoM.m`, `edge_conx_elem.m`, `edgemake_MoM.m`, `find_local_dofs.m`, `hat.m`, `intg_sing_SGF.m`, `outside_edge.m`, `renumber_RWG.m`, `simplex_area.m`, `tri_area3D.m`, `tri_quad.m` e `trimesh3D.m` -> traduzidos em `src/`.

As principais adaptações desta versão em C++ são:

- interface por linha de comando no lugar dos `input(...)` do MATLAB;
- gravação de CSVs em `Cap_06/out/` para inspeção e plot;
- geração das figuras por script Python, em vez de `plot(...)` e `semilogy(...)` diretos;
- exportação explícita de malha, correntes e cortes para facilitar validação.

## Teoria

Na parte analítica, `sphereRCS` implementa a série modal da RCS monostática de uma esfera PEC. O script `converge_sphereRCS` repete a avaliação para vários truncamentos `N` e reproduz a ideia da Fig. 6.5 do livro.

Na parte numérica, `MoM3D_demo` segue a formulação de Rao, Wilton e Glisson para a EFIE de potencial misto em superfícies trianguladas. Os graus de liberdade vivem nas arestas internas da malha e as integrais de acoplamento são montadas de duas formas:

- por par de arestas, em `FillZMatrixByEdge`;
- por par de faces, em `FillZMatrixByFace`.

O capítulo também preserva o tratamento especial do termo singular via `intg_sing_SGF`, baseado na transformação descrita por Khayat e Wilton.

## Estrutura

```text
Cap_06/
├── CMakeLists.txt
├── MoM3D_demo.cpp
├── converge_sphereRCS.cpp
├── include/
│   ├── *.hpp
├── src/
│   ├── *.cpp
├── scripts/
│   └── plot_cap06.py
├── out/
│   ├── mom3d_*.csv
│   ├── sphere_rcs_*.csv
│   └── *.png
└── README.md
```

## Arquivos MATLAB de Referência

Arquivos usados como base local para a tradução:

- `original_matlab/Chapter 6/Mixed Potential EFIE RWG/MoM3D_demo.m`
- `original_matlab/Chapter 6/Mixed Potential EFIE RWG/ComputeRho_c.m`
- `original_matlab/Chapter 6/Mixed Potential EFIE RWG/FillVVector.m`
- `original_matlab/Chapter 6/Mixed Potential EFIE RWG/FillZMatrixByEdge.m`
- `original_matlab/Chapter 6/Mixed Potential EFIE RWG/FillZMatrixByFace.m`
- `original_matlab/Chapter 6/Mixed Potential EFIE RWG/GL_quad_rule.m`
- `original_matlab/Chapter 6/Mixed Potential EFIE RWG/Int_pq.m`
- `original_matlab/Chapter 6/Mixed Potential EFIE RWG/PostProcMoM.m`
- `original_matlab/Chapter 6/Mixed Potential EFIE RWG/edge_conx_elem.m`
- `original_matlab/Chapter 6/Mixed Potential EFIE RWG/edgemake_MoM.m`
- `original_matlab/Chapter 6/Mixed Potential EFIE RWG/find_local_dofs.m`
- `original_matlab/Chapter 6/Mixed Potential EFIE RWG/hat.m`
- `original_matlab/Chapter 6/Mixed Potential EFIE RWG/intg_sing_SGF.m`
- `original_matlab/Chapter 6/Mixed Potential EFIE RWG/outside_edge.m`
- `original_matlab/Chapter 6/Mixed Potential EFIE RWG/renumber_RWG.m`
- `original_matlab/Chapter 6/Mixed Potential EFIE RWG/simplex_area.m`
- `original_matlab/Chapter 6/Mixed Potential EFIE RWG/tri_area3D.m`
- `original_matlab/Chapter 6/Mixed Potential EFIE RWG/tri_quad.m`
- `original_matlab/Chapter 6/Mixed Potential EFIE RWG/trimesh3D.m`
- `original_matlab/Chapter 6/Sphere RCS/converge_sphereRCS.m`
- `original_matlab/Chapter 6/Sphere RCS/sphereRCS.m`

## Compilação

```bash
cd Cap_06
cmake -S . -B build
cmake --build build -j$(nproc)
```

Executáveis gerados:

```bash
./build/MoM3D_demo
./build/converge_sphereRCS
```

## Execução

Convergência da RCS da esfera:

```bash
./build/converge_sphereRCS
```

Caso RWG82 equivalente à Fig. 6:

```bash
./build/MoM3D_demo
```

Caso RWG82 equivalente à Fig. 5 com integral singular:

```bash
./build/MoM3D_demo --prob-type 5 --sing 1
```

Teste mais leve:

```bash
./build/MoM3D_demo --prob-type 6 --xmesh 2 --ymesh 3
```

## Saídas

Os arquivos são gerados em `Cap_06/out/`:

- `sphere_rcs_convergence.csv` e `sphere_rcs_summary.csv`;
- `mom3d_summary.csv`;
- `mom3d_currents.csv`;
- `mom3d_cuts.csv`;
- `mom3d_mesh_nodes.csv`;
- `mom3d_mesh_elements.csv`;
- `mom3d_postproc_summary.csv`.

## Visualização

```bash
cd Cap_06
python3 scripts/plot_cap06.py
```

Figuras geradas:

- `out/sphere_rcs_convergence.png`
- `out/mom3d_cuts.png`

## Limitações e Observações

- A implementação mantém a malha especializada do original, restrita a uma placa retangular em `z = 0`.
- A excitação incidente segue restrita ao caso normal com campo elétrico dirigido em `x`, como no MATLAB.
- O tratamento singular foi traduzido de forma direta; ele é o trecho numericamente mais sensível do capítulo e merece validação fina caso a caso.

## Ligação com o Repositório

Este é o módulo correspondente ao **Capítulo 6**. Consulte o [README principal](../README.md) e a [matriz de tradução](../TRANSLATION_MATRIX.md) para acompanhar o restante do projeto.
