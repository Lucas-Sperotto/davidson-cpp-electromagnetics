#ifndef CAP_09_FEM_PP_H
#define CAP_09_FEM_PP_H

#include <vector>

std::vector<double> linspace(double a, double b, int N);

// Pos-processamento: interpola V nos N_int pontos em [0, ell], retorna V_fem e escreve z_out.
std::vector<double> FEM_postprocess(const std::vector<double> &V_tot,
                                    int N_elem, double ell, double h,
                                    int N_int, std::vector<double> &z_out);

#endif // CAP_09_FEM_PP_H
