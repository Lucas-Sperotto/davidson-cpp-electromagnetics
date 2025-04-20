// src/fdtd2d_pml.cpp
#include "../include/fdtd2d_pml.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "gaussder.cpp"

const double c = 2.997925e8;
const double eps_0 = 8.854e-12;
const double mu_0 = 4 * M_PI * 1e-7;
const double eta_0 = std::sqrt(mu_0 / eps_0);

void run_fdtd2d_pml_simulation() {
    int refine = 1;
    int N_x = refine * 200;
    int N_y = refine * 200;
    int M = refine * 350;
    int L = N_x / 2;
    int d_cell = 10;

    double delta_s = 0.005 / refine;
    double delta_t = delta_s / (c * std::sqrt(2.0));
    double sigma = 1.0e-10;
    double m_offset = 4 * sigma;
    double Peak = 1.0;

    std::vector<std::vector<double>> E_x(N_x, std::vector<double>(N_y, 0.0));
    std::vector<std::vector<double>> E_y(N_x, std::vector<double>(N_y, 0.0));
    std::vector<std::vector<double>> H_z(N_x, std::vector<double>(N_y, 0.0));
    std::vector<std::vector<double>> H_zx(N_x, std::vector<double>(N_y, 0.0));
    std::vector<std::vector<double>> H_zy(N_x, std::vector<double>(N_y, 0.0));

    std::vector<std::vector<double>> sigma_x(N_x, std::vector<double>(N_y, 0.0));
    std::vector<std::vector<double>> sigma_y(N_x, std::vector<double>(N_y, 0.0));
    std::vector<std::vector<double>> sigma_ast_x(N_x, std::vector<double>(N_y, 0.0));
    std::vector<std::vector<double>> sigma_ast_y(N_x, std::vector<double>(N_y, 0.0));

    double sigma_max = 0.8 * 4 / (eta_0 * delta_s);
    int poly_m = 3;

    for (int jj = 0; jj < d_cell; ++jj) {
        double dist = d_cell - jj;
        double val = sigma_max * std::pow(dist / d_cell, poly_m);
        for (int i = 0; i < N_x; ++i) {
            sigma_y[i][jj] = val;
            sigma_y[i][N_y - jj - 1] = val;
        }
        dist -= 0.5;
        double val_ast = eta_0 * eta_0 * sigma_max * std::pow(dist / d_cell, poly_m);
        for (int i = 0; i < N_x; ++i) {
            sigma_ast_y[i][jj] = val_ast;
            sigma_ast_y[i][N_y - jj - 1] = val_ast;
        }
    }
    for (int ii = 0; ii < d_cell; ++ii) {
        double dist = d_cell - ii;
        double val = sigma_max * std::pow(dist / d_cell, poly_m);
        for (int j = 0; j < N_y; ++j) {
            sigma_x[ii][j] = val;
            sigma_x[N_x - ii - 1][j] = val;
        }
        dist -= 0.5;
        double val_ast = eta_0 * eta_0 * sigma_max * std::pow(dist / d_cell, poly_m);
        for (int j = 0; j < N_y; ++j) {
            sigma_ast_x[ii][j] = val_ast;
            sigma_ast_x[N_x - ii - 1][j] = val_ast;
        }
    }

    std::vector<std::vector<double>> C_aEx(N_x, std::vector<double>(N_y));
    std::vector<std::vector<double>> C_bEx(N_x, std::vector<double>(N_y));
    std::vector<std::vector<double>> C_aEy(N_x, std::vector<double>(N_y));
    std::vector<std::vector<double>> C_bEy(N_x, std::vector<double>(N_y));
    std::vector<std::vector<double>> D_aHzx(N_x, std::vector<double>(N_y));
    std::vector<std::vector<double>> D_bHzx(N_x, std::vector<double>(N_y));
    std::vector<std::vector<double>> D_aHzy(N_x, std::vector<double>(N_y));
    std::vector<std::vector<double>> D_bHzy(N_x, std::vector<double>(N_y));

    for (int i = 0; i < N_x; ++i) {
        for (int j = 0; j < N_y; ++j) {
            C_aEx[i][j] = (1 - sigma_y[i][j] * delta_t / (2 * eps_0)) /
                          (1 + sigma_y[i][j] * delta_t / (2 * eps_0));
            C_bEx[i][j] = delta_t / (eps_0 * delta_s) /
                          (1 + sigma_y[i][j] * delta_t / (2 * eps_0));
            C_aEy[i][j] = (1 - sigma_x[i][j] * delta_t / (2 * eps_0)) /
                          (1 + sigma_x[i][j] * delta_t / (2 * eps_0));
            C_bEy[i][j] = delta_t / (eps_0 * delta_s) /
                          (1 + sigma_x[i][j] * delta_t / (2 * eps_0));
            D_aHzx[i][j] = (1 - sigma_ast_x[i][j] * delta_t / (2 * mu_0)) /
                           (1 + sigma_ast_x[i][j] * delta_t / (2 * mu_0));
            D_bHzx[i][j] = delta_t / (mu_0 * delta_s) /
                           (1 + sigma_ast_x[i][j] * delta_t / (2 * mu_0));
            D_aHzy[i][j] = (1 - sigma_ast_y[i][j] * delta_t / (2 * mu_0)) /
                           (1 + sigma_ast_y[i][j] * delta_t / (2 * mu_0));
            D_bHzy[i][j] = delta_t / (mu_0 * delta_s) /
                           (1 + sigma_ast_y[i][j] * delta_t / (2 * mu_0));
        }
    }

    int point1_x = 190;
    int point1_y = N_y / 2;
    std::vector<double> E_y_point1(M);

    for (int m = 1; m < M; ++m) {
        // Atualizar H
        for (int i = 0; i < N_x - 1; ++i) {
            for (int j = 0; j < N_y - 1; ++j) {
                H_zx[i][j] = D_aHzx[i][j] * H_zx[i][j] - D_bHzx[i][j] * (E_y[i + 1][j] - E_y[i][j]);
                H_zy[i][j] = D_aHzy[i][j] * H_zy[i][j] + D_bHzy[i][j] * (E_x[i][j + 1] - E_x[i][j]);
            }
        }
        for (int i = 0; i < N_x; ++i) {
            for (int j = 0; j < N_y; ++j) {
                H_z[i][j] = H_zx[i][j] + H_zy[i][j];
            }
        }

        // Atualizar E
        for (int i = 1; i < N_x; ++i) {
            for (int j = 1; j < N_y; ++j) {
                E_x[i][j] = C_aEx[i][j] * E_x[i][j] + C_bEx[i][j] * (H_z[i][j] - H_z[i][j - 1]);
                E_y[i][j] = C_aEy[i][j] * E_y[i][j] - C_bEy[i][j] * (H_z[i][j] - H_z[i - 1][j]);
            }
        }

        // Fonte
        double inc = Peak * gaussder_norm((m - 1) * delta_t - (L - 1) * delta_s / c, m_offset, sigma);
        E_y[L][point1_y] += inc;

        // Registro
        E_y_point1[m] = E_y[point1_x][point1_y];
    }

    std::ofstream file("out/ey_point1_pml.csv");
    file << "tempo_ns,Ey_Vpm\n";
    for (int m = 0; m < M; ++m) {
        double time_ns = m * delta_t * 1e9;
        file << time_ns << "," << E_y_point1[m] << "\n";
    }
    file.close();
    std::cout << "Simulação com PML finalizada. Resultados em out/ey_point1_pml.csv" << std::endl;
}