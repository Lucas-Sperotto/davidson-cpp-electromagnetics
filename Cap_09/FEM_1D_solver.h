#ifndef CAP_09_FEM_1D_SOLVER_H
#define CAP_09_FEM_1D_SOLVER_H

#include <vector>

// Solver FEM 1D (linha de transmissao) com elementos lineares e BC: V(N)=V_in.
// Retorna V_tot com N_elem+1 nos.
std::vector<double> FEM_1D_solver(int N_elem, double h, double omega,
                                  double L, double C, double V_in);

#endif // CAP_09_FEM_1D_SOLVER_H
