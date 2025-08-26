#include <vector>
#include <Eigen/Dense>
#include "FEM_1D_solver.h"

// Solver FEM 1D (transmission line) com elementos lineares e BC: V(N)=V_in
// Retorna V_tot (N_elem+1 nós)
std::vector<double> FEM_1D_solver(int N_elem, double h, double omega,
                                  double L, double C, double V_in)
{
    const int N = N_elem + 1;                          // nós
    const double alpha = 1.0 / (h * L);                // escala de S
    const double gamma = -omega * omega * h * C / 6.0; // escala de T
    // Sistema reduzido (primeiros N-1 nós). Vetor força F reduzido.
    const int n = N - 1;

    Eigen::MatrixXd S = Eigen::MatrixXd::Zero(N, N);
    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(N, N);

    for (int i = 0; i < N; ++i)
    {
        S(i, i) = 2.0;
        T(i, i) = 4.0;
        if (i > 0)
        {
            S(i, i - 1) = -1.0;
            T(i, i - 1) = 1.0;
        }
        if (i < (N - 1))
        {
            S(i, i + 1) = -1.0;
            T(i, i + 1) = 1.0;
        }
    }

    S(0, 0) = 1.0; // Set first and last diagonal elements
    S(N - 1, N - 1) = 1.0;
    S *= alpha; // Scale entries

    T(0, 0) = 2.0;
    T(N - 1, N - 1) = 2.0;
    T *= gamma;

    // Monta A (n x n) e F (n)
    Eigen::MatrixXd M =  S + T;

    Eigen::VectorXd F = Eigen::VectorXd::Zero(n); 
    F(n - 1) = -M(N - 2, N - 1) * V_in;               // Voltage at source end.

    Eigen::MatrixXd A = M.topLeftCorner(n, n);

    Eigen::PartialPivLU<Eigen::MatrixXd> lu(A);
    if (lu.matrixLU().diagonal().array().abs().minCoeff() < 1e-15)
    {
        throw std::runtime_error("Sistema quase singular.");
    }
    Eigen::VectorXd V_reduced = lu.solve(F);
    if ((A * V_reduced - F).norm() > 1e-10 * F.norm())
    {
        throw std::runtime_error("Solução instável ou divergente.");
    }

    // Reconstrói V_tot com condição de contorno em N-1
    std::vector<double> V_tot(N, 0.0);
    for (int i = 0; i < n; ++i)
        V_tot[i] = V_reduced(i);
    V_tot[N - 1] = V_in;
    return V_tot;
}