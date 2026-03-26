# Translation Matrix

Matriz de acompanhamento entre os arquivos MATLAB de referência e seus equivalentes atuais no repositório C++.

Legenda de status:

- `Traduzido` - existe equivalente funcional no repositório;
- `Traduzido com adaptação` - existe equivalente, mas com mudanças estruturais relevantes;
- `Pendente` - ainda não há equivalente claro no repositório.

## Capítulo 2

| MATLAB original | Equivalente atual | Status | Observação |
| --- | --- | --- | --- |
| `fdtd_1D_demo.m` | `Cap_02/src/fdtd_1D_demo.cpp` | Traduzido | Demo principal do capítulo |
| `fdtd_1D_WB_demo.m` | `Cap_02/src/fdtd_1D_WB_demo.cpp` | Traduzido | Demo banda larga |
| `gaussder_norm.m` | `Cap_02/src/gaussder.cpp` / `Cap_02/include/gaussder.hpp` | Traduzido | Função isolada como módulo próprio |

## Capítulo 3

| MATLAB original | Equivalente atual | Status | Observação |
| --- | --- | --- | --- |
| `FDTD_2D/fdtd_2d_demo.m` | `Cap_03/src/fdtd_2d_demo.cpp` | Traduzido com adaptação | Alguns parâmetros estão fixos e há trechos ainda marcados como não testados |
| `FDTD_2D/fdtd_2d_pml_demo.m` | `Cap_03/src/fdtd2d_pml.cpp` | Traduzido com adaptação | Estrutura e parâmetros preservados; a injeção scat/tot do PML segue em formulação estável adaptada |
| `FDTD_3D/fdtd_3D_demo.m` | `Cap_03/src/fdtd3d.cpp` | Traduzido com adaptação | O cenário principal foi preservado |
| `FDTD_2D/gaussder_norm.m` | `Cap_03/src/gaussder.cpp` | Traduzido | Precisa apenas de padronização estrutural |
| `FDTD_2D/PMLperformance.m` | `Cap_03/scripts/pml_performance.py` | Traduzido com adaptação | Traduzido como script Python para comparação entre rodadas |
| `FDTD_2D/cyl_test.m` | `Cap_03/scripts/cyl_test.py` | Traduzido com adaptação | Traduzido como script Python de comparação temporal |

## Capítulo 4

| MATLAB original | Equivalente atual | Status | Observação |
| --- | --- | --- | --- |
| `MoM_2D_TM.m` | `Cap_04/MoM_2D_TM.cpp` | Traduzido com adaptação | Interface interativa substituída por CLI e saídas CSV |
| `MoM_TM_solver.m` | `Cap_04/src/mom_tm_solver.cpp` / `Cap_04/include/mom_tm_solver.hpp` | Traduzido | Solver principal do espalhamento TM |
| `cyl_TM_echo_width.m` | `Cap_04/src/cyl_tm_echo_width.cpp` / `Cap_04/include/cyl_tm_echo_width.hpp` | Traduzido | Solução analítica para eco TM de cilindro PEC |
| `static_mom.m` | `Cap_04/static_mom.cpp` | Traduzido com adaptação | Script convertido para executável com saída CSV |
| `thin_dipole.m` | `Cap_04/thin_dipole.cpp` | Traduzido com adaptação | Entradas interativas convertidas para CLI e saídas CSV |

## Capítulo 9

| MATLAB original | Equivalente atual | Status | Observação |
| --- | --- | --- | --- |
| `FEM_1D_1st.m` | `Cap_09/FEM_1D_1st.cpp` | Traduzido | Programa principal |
| `FEM_1D_solver.m` | `Cap_09/FEM_1D_solver.cpp` / `Cap_09/FEM_1D_solver.h` | Traduzido | Solver modularizado |
| `FEM_pp.m` | `Cap_09/FEM_pp.cpp` / `Cap_09/FEM_pp.h` | Traduzido | Pós-processamento modularizado |

## Capítulo 10

| MATLAB original | Equivalente atual | Status | Observação |
| --- | --- | --- | --- |
| `Static2D.m` | `Cap_10/Static2D.cpp` | Traduzido | Solver principal |
| `Eigen2D.m` | `Cap_10/Eigen2D.cpp` | Traduzido | Solver principal com Whitney |
| `Eigen2D_LTQN.m` | `Cap_10/Eigen2D_LTQN.cpp` | Traduzido | Solver principal com LTQN |
| `sandt.m` | `Cap_10/src/sandt.cpp` | Traduzido | Matriz local Whitney |
| `sandt_LTQN.m` | `Cap_10/src/sandt_LTQN.cpp` | Traduzido | Matriz local LTQN |
| `LTQN.m` | `Cap_10/src/LTQN.cpp` | Traduzido | Base LTQN |
| `curl_LTQN.m` | `Cap_10/src/curl_LTQN.cpp` | Traduzido | Rotacional LTQN |
| `whitney.m` | `Cap_10/src/whitney.cpp` | Traduzido | Base de Whitney |
| `plot_field.m` | `Cap_10/src/plot_field.cpp` | Traduzido | Pós-processamento de campo |
| `TEeig_err.m` | `Cap_10/src/TEeig_err.cpp` | Traduzido | Avaliação de erro |
| `trimesh.m` | `Cap_10/src/trimesh.cpp` | Traduzido | Geração de malha |
| `edgemake.m` | `Cap_10/src/edgemake.cpp` | Traduzido | Geração de arestas |
| `find_local_dofs.m` | `Cap_10/src/find_local_dofs.cpp` | Traduzido | Topologia local |
| `free_dof.m` | `Cap_10/src/free_dof.cpp` | Traduzido | Contorno PEC em arestas |
| `free_nodes.m` | `Cap_10/src/free_nodes.cpp` | Traduzido | Contorno em nós |
| `free_nodes_mstrip.m` | `Cap_10/src/free_nodes_mstrip.cpp` | Traduzido | Contorno do microstrip |
| `prescr_nodes_mstrip.m` | `Cap_10/src/prescr_nodes_mstrip.cpp` | Traduzido | Potenciais prescritos |
| `renumber_dof.m` | `Cap_10/src/renumber_dof.cpp` | Traduzido | Renumeração de DOFs |
| `renumber_dof_LTQN.m` | `Cap_10/src/renumber_dof_LTQN.cpp` | Traduzido | Renumeração LTQN |
| `simplex2D.m` | `Cap_10/src/simplex2D.cpp` | Traduzido | Coordenadas baricêntricas |
| `tri_quad.m` | `Cap_10/src/tri_quad.cpp` | Traduzido | Quadratura em triângulo |
| `nodal_1st.m` | `Cap_10/src/nodal_1st.cpp` | Traduzido | Base nodal |
| `s_nodal.m` | `Cap_10/src/s_nodal.cpp` | Traduzido | Matriz nodal |

## Uso Prático da Matriz

Próximo uso recomendado desta matriz:

1. associar cada linha a um caso de validação numérica;
2. anexar, no futuro, um campo de evidência com CSV, figura ou comparação MATLAB versus C++;
3. usar o status `Pendente` para ordenar a continuação do trabalho.
