# Davidson C++ Electromagnetics

This repository contains C++ translations of selected MATLAB scripts by Prof. D. B. Davidson, as described in his book:

**"Computational Electromagnetics for RF and Microwave Engineering"**, 2nd Edition, Cambridge University Press.

🔗 [View book on Cambridge University Press](https://www.cambridge.org/br/universitypress/subjects/engineering/rf-and-microwave-engineering/computational-electromagnetics-rf-and-microwave-engineering-2nd-edition?format=HB&isbn=9780521518918)

These C++ implementations aim to preserve the educational value of the original MATLAB versions while providing an accessible path for those working in C++ environments.

## 🔍 Chapters

- [`Cap_10/`](Cap_10/) – C++ translations of Chapter 10 codes ([Capítulo 10 README](Cap_10/README.md))

## 📂 Cap_10 – Chapter 10 Files

The files in this folder are C++ translations of the original MATLAB scripts discussed in Chapter 10 of Davidson's book. They include:

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

## 📄 Disclaimer

This project is unofficial and is not affiliated with Prof. D. B. Davidson or Cambridge University Press. The translations are provided for educational and research purposes only.

Dependências: C++ sudo apt install build-essential gfortran -y
Cmake sudo apt install cmake -y
pkg-config sudo apt install pkg-config -y
python sudo apt install python3 python3-pip -y
numpy sudo apt install python3-numpy -y
matplotlib sudo apt install python3-matplotlib -y
pandas sudo apt install python3-pandas -y
pillow sudo apt install python3-pil