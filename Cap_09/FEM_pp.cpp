#include <vector>
#include <cmath>
#include "FEM_pp.h"

// Linspace [a,b] com N pontos
std::vector<double> linspace(double a, double b, int N)
{
    std::vector<double> x(N);
    if (N == 1)
    {
        x[0] = a;
        return x;
    }
    double step = (b - a) / (N - 1);
    for (int i = 0; i < N; ++i)
        x[i] = a + step * (double)i;
    return x;
}

// Pós-processamento (interpolação FEM linear) em N_int pontos ao longo de [0, ell]
// Retorna V_fem e escreve z_out
std::vector<double> FEM_postprocess(const std::vector<double> &V_tot,
                                    int N_elem, double ell, double h,
                                    int N_int, std::vector<double> &z_out)
{
    z_out = linspace(0.0, ell, N_int);
    std::vector<double> V_fem(N_int, 0.0);
    for (int i = 0; i < N_int; ++i)
    {
        for (int j = 0; j < N_elem; ++j)
        {
            double z_l = j * h;
            double z_r = (j + 1) * h;
            if (z_l <= z_out[i] && z_out[i] <= z_r)
            {
                double phi_l = (z_r - z_out[i]) / h;
                double phi_r = (z_out[i] - z_l) / h;
                V_fem[i] = V_tot[j] * phi_l + V_tot[j + 1] * phi_r;
            }
        }
    }
    return V_fem;
}