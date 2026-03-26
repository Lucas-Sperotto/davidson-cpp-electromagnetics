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

## Capítulo 6

| MATLAB original | Equivalente atual | Status | Observação |
| --- | --- | --- | --- |
| `Mixed Potential EFIE RWG/MoM3D_demo.m` | `Cap_06/MoM3D_demo.cpp` | Traduzido com adaptação | Prompt interativo substituído por CLI e saídas CSV |
| `Mixed Potential EFIE RWG/ComputeRho_c.m` | `Cap_06/src/ComputeRho_c.cpp` | Traduzido | Vetores `rho_c^+` e `rho_c^-` |
| `Mixed Potential EFIE RWG/FillVVector.m` | `Cap_06/src/FillVVector.cpp` | Traduzido | Montagem do vetor de excitação |
| `Mixed Potential EFIE RWG/FillZMatrixByEdge.m` | `Cap_06/src/FillZMatrixByEdge.cpp` | Traduzido | Montagem por par de arestas |
| `Mixed Potential EFIE RWG/FillZMatrixByFace.m` | `Cap_06/src/FillZMatrixByFace.cpp` | Traduzido | Montagem por par de faces |
| `Mixed Potential EFIE RWG/GL_quad_rule.m` | `Cap_06/src/GL_quad_rule.cpp` | Traduzido | Quadratura de Gauss-Legendre em linha |
| `Mixed Potential EFIE RWG/Int_pq.m` | `Cap_06/src/Int_pq.cpp` | Traduzido | Integrais sobre triângulo fonte |
| `Mixed Potential EFIE RWG/PostProcMoM.m` | `Cap_06/src/PostProcMoM.cpp` | Traduzido com adaptação | Pós-processamento exportado em CSV |
| `Mixed Potential EFIE RWG/edge_conx_elem.m` | `Cap_06/src/edge_conx_elem.cpp` | Traduzido | Conectividade aresta-elemento |
| `Mixed Potential EFIE RWG/edgemake_MoM.m` | `Cap_06/src/edgemake_MoM.cpp` | Traduzido | Geração das arestas globais |
| `Mixed Potential EFIE RWG/find_local_dofs.m` | `Cap_06/src/find_local_dofs.cpp` | Traduzido | Numeração local das arestas por DOF |
| `Mixed Potential EFIE RWG/hat.m` | `Cap_06/src/hat.cpp` | Traduzido | Vetor unitário |
| `Mixed Potential EFIE RWG/intg_sing_SGF.m` | `Cap_06/src/intg_sing_SGF.cpp` | Traduzido | Integração singular por transformação |
| `Mixed Potential EFIE RWG/outside_edge.m` | `Cap_06/src/outside_edge.cpp` | Traduzido | Marca arestas de contorno |
| `Mixed Potential EFIE RWG/renumber_RWG.m` | `Cap_06/src/renumber_RWG.cpp` | Traduzido | Renumeração dos DOFs RWG |
| `Mixed Potential EFIE RWG/simplex_area.m` | `Cap_06/src/simplex_area.cpp` | Traduzido | Coordenadas simplex |
| `Mixed Potential EFIE RWG/tri_area3D.m` | `Cap_06/src/tri_area3D.cpp` | Traduzido | Área de triângulo 3D |
| `Mixed Potential EFIE RWG/tri_quad.m` | `Cap_06/src/tri_quad.cpp` | Traduzido | Quadratura de Dunavant |
| `Mixed Potential EFIE RWG/trimesh3D.m` | `Cap_06/src/trimesh3D.cpp` | Traduzido | Malha triangular da placa PEC |
| `Sphere RCS/converge_sphereRCS.m` | `Cap_06/converge_sphereRCS.cpp` | Traduzido com adaptação | Script convertido para executável com saída CSV |
| `Sphere RCS/sphereRCS.m` | `Cap_06/src/sphereRCS.cpp` | Traduzido | Solução analítica modal da RCS |

## Capítulo 7

| MATLAB original | Equivalente atual | Status | Observação |
| --- | --- | --- | --- |
| `scalar_pot.m` | `Cap_07/scalar_pot.cpp` | Traduzido com adaptação | Gera CSVs para as figuras das integrais de Sommerfeld |
| `V_pot_eps.m` | `Cap_07/V_pot_eps.cpp` | Traduzido | Estudo do potencial para diferentes permissividades |
| `V_pot_height.m` | `Cap_07/V_pot_height.cpp` | Traduzido | Estudo do potencial para diferentes alturas |
| `MoM_Som.m` | `Cap_07/MoM_Som.cpp` | Traduzido com adaptação | Prompt interativo substituído por CLI |
| `D_TE.m` | `Cap_07/src/D_TE.cpp` | Traduzido | Denominador TE |
| `D_TM.m` | `Cap_07/src/D_TM.cpp` | Traduzido | Denominador TM |
| `F.m` | `Cap_07/src/F.cpp` | Traduzido | Integrando original |
| `F_reg1.m` | `Cap_07/src/F_reg1.cpp` | Traduzido | Região 1 |
| `F_reg2.m` | `Cap_07/src/F_reg2.cpp` | Traduzido | Região 2 com singularidade extraída |
| `F_reg3.m` | `Cap_07/src/F_reg3.cpp` | Traduzido | Região 3 com termo estático extraído |
| `F_static.m` | `Cap_07/src/F_static.cpp` | Traduzido | Termo estático |
| `V_int.m` | `Cap_07/src/V_int.cpp` | Traduzido | Integração do potencial escalar |
| `root_D_TM.m` | `Cap_07/src/root_D_TM.cpp` | Traduzido | Busca do polo TM |

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

## Capítulo 11

| MATLAB original | Equivalente atual | Status | Observação |
| --- | --- | --- | --- |
| `Eigen3D_v0.m` | `Cap_11/Eigen3D_v0.cpp` | Traduzido com adaptação | Geração interna de tetraedros feita por divisão determinística do brick |
| `Eigen3D_CTLN.m` | `Cap_11/Eigen3D_CTLN.cpp` | Traduzido com adaptação | Suporta malha interna e leitura Gmsh |
| `Eigen3D_LTQN.m` | `Cap_11/Eigen3D_LTQN.cpp` | Traduzido com adaptação | CLI substitui prompts e suporta ordem 1 e 2 |
| `FETD_FDTD.m` | `Cap_11/FETD_FDTD.cpp` | Traduzido com adaptação | Saída em CSV no lugar do plot direto |
| `LTQN3D.m` | `Cap_11/src/LTQN3D.cpp` | Traduzido | Base LTQN 3D |
| `TE_cavity_modes3D.m` | `Cap_11/src/TE_cavity_modes3D.cpp` | Traduzido | Modos TE analíticos |
| `avg_mesh_length.m` | `Cap_11/src/avg_mesh_length.cpp` | Traduzido | Métrica média da malha |
| `brick_mesh.m` | `Cap_11/src/brick_mesh.cpp` | Traduzido | Geração da malha hexaédrica regular |
| `cavity_modes3D.m` | `Cap_11/src/cavity_modes3D.cpp` | Traduzido | Mantém o estado incompleto do original |
| `curl_LTQN3D.m` | `Cap_11/src/curl_LTQN3D.cpp` | Traduzido | Rotacional LTQN 3D |
| `edgemake3D.m` | `Cap_11/src/edgemake3D.cpp` | Traduzido | Geração de arestas 3D |
| `eig_err3D.m` | `Cap_11/src/eig_err3D.cpp` | Traduzido | Erro dos autovalores |
| `facemake3D.m` | `Cap_11/src/facemake3D.cpp` | Traduzido | Geração de faces 3D |
| `free_dof3D.m` | `Cap_11/src/free_dof3D.cpp` | Traduzido | Contorno CTLN |
| `free_dof3D_LTQN.m` | `Cap_11/src/free_dof3D_LTQN.cpp` | Traduzido | Contorno LTQN |
| `free_nodes3D.m` | `Cap_11/src/free_nodes3D.cpp` | Traduzido | Nós livres para o núcleo nulo |
| `read_gmsh2.m` | `Cap_11/src/read_gmsh2.cpp` | Traduzido | Leitor Gmsh v2 ASCII |
| `renumber_dof.m` | `Cap_11/src/renumber_dof3D.cpp` | Traduzido | Renumeração CTLN |
| `renumber_dof_LTQN.m` | `Cap_11/src/renumber_dof3D_LTQN.cpp` | Traduzido | Renumeração LTQN |
| `sandt3D.m` | `Cap_11/src/sandt3D.cpp` | Traduzido | Matrizes locais 3D de Whitney |
| `sandt3D_LTQN.m` | `Cap_11/src/sandt3D_LTQN.cpp` | Traduzido | Matrizes locais 3D LTQN |
| `test_tet_quad.m` | `Cap_11/test_tet_quad.cpp` | Traduzido com adaptação | Resultado exportado em CSV |
| `tet_quad.m` | `Cap_11/src/tet_quad.cpp` | Traduzido | Regra de quadratura tetraédrica |

## Uso Prático da Matriz

Próximo uso recomendado desta matriz:

1. associar cada linha a um caso de validação numérica;
2. anexar, no futuro, um campo de evidência com CSV, figura ou comparação MATLAB versus C++;
3. usar o status `Pendente` para ordenar a continuação do trabalho.
