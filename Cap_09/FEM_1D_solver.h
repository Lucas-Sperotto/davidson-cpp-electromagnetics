#include <vector>

// Solver FEM 1D (linha de transmissão) com elementos lineares e BC: V(N)=V_in
// Retorna V_tot (N_elem+1 nos)
std::vector<double> FEM_1D_solver(int N_elem, double h, double omega,
                                  double L, double C, double V_in);