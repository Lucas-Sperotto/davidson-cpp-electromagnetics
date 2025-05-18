# Capítulo 10 — Códigos Traduzidos para C++

Este diretório contém as versões em C++ dos códigos MATLAB apresentados no **Capítulo 10** do livro:

**D.B. Davidson**, *Computational Electromagnetics for RF and Microwave Engineering*, Cambridge University Press, 2ª edição (2009).

🔗 [Acesse o livro na Cambridge University Press](https://www.cambridge.org/br/universitypress/subjects/engineering/rf-and-microwave-engineering/computational-electromagnetics-rf-and-microwave-engineering-2nd-edition?format=HB&isbn=9780521518918)

---
### Main Programs
- `Static2D.cpp` – Computes the quasi-TEM mode in boxed microstrip
- `Eigen2D.cpp` – Computes TE eigenvalues and eigenmodes of hollow rectangular waveguide

### Support Functions
- `free_nodes.cpp`, `free_nodes_mstrip.cpp` – Identify free nodes based on PEC boundaries
- `free_dof.cpp` – Identifies free edge DOFs
- `find_local_dofs.cpp` – Locates local edge indices per DOF
- `renumber_dof.cpp`, `renumber_dof_LTQN.cpp` – Renumber DOFs for Whitney and LTQN elements
- `sandt.cpp`, `s_nodal.cpp`, `sandt_LTQN.cpp` – Compute local stiffness/mass matrices
- `simplex2D.cpp` – Computes barycentric coordinates
- `edgemake.cpp` – Builds global edge list from mesh
- `whitney.cpp` – Whitney 1-form basis functions
- `LTQN.cpp`, `curl_LTQN.cpp` – LTQN basis functions and their curls
- `plot_field.cpp` – Evaluates and outputs vector field from DOFs
- `TEeig_err.cpp` – Evaluates relative error of eigenvalues
- `trimesh.cpp` – Generates structured triangular mesh

All auxiliary functions are modular and reusable.

## Arquivos principais

- `Static2D.cpp`: calcula o modo quasi-TEM em microfita encapsulada (boxed microstrip)
- `Eigen2D.cpp`: calcula os autovalores e modos próprios TE de um guia de onda retangular oco

## Arquivos auxiliares

As funções auxiliares implementam geração de malha, montagem de matrizes, avaliação de campo, elementos de Whitney e LTQN, entre outros:

- Geração de malha: `trimesh.cpp`, `edgemake.cpp`
- Geometria e topologia: `simplex2D.cpp`, `find_local_dofs.cpp`, `renumber_dof.cpp`, `renumber_dof_LTQN.cpp`
- Elementos: `whitney.cpp`, `LTQN.cpp`, `curl_LTQN.cpp`
- Matrizes locais: `sandt.cpp`, `s_nodal.cpp`, `sandt_LTQN.cpp`
- Condições de contorno: `free_nodes.cpp`, `free_nodes_mstrip.cpp`, `free_dof.cpp`, `prescr_nodes_mstrip.cpp`
- Pós-processamento: `plot_field.cpp`, `TEeig_err.cpp`

---

## Observação

Todos os arquivos foram fielmente traduzidos para manter a lógica original dos scripts MATLAB de Davidson. Os nomes de funções e estrutura de dados foram modernizados e modularizados para melhor integração em projetos C++.

Este repositório é de uso educacional e não é afiliado oficialmente ao autor ou à editora.

