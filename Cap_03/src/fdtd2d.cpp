// src/fdtd2d.cpp
#include "../include/fdtd2d.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <filesystem>
#include "gaussder.cpp"

// Constantes físicas
const double c = 2.997925e8;
const double eps_0 = 8.854e-12;
const double mu_0 = 4 * M_PI * 1e-7;
const double eta_0 = std::sqrt(mu_0 / eps_0);

void run_fdtd2d_simulation() {

    const std::string out_dir = PROJECT_OUT_DIR;
    std::filesystem::create_directories(out_dir);

    int refine = 1;
    int N_x = refine * 200;
    int N_y = refine * 200;
    int M = refine * 512;
    int L = N_x / 2;
    int L1 = 20;

    double delta_s = 0.005 * 2 / refine;
    double R = 1.0;
    double delta_t = R * delta_s / (c * std::sqrt(2.0));
    double sigma = 1.0e-10;
    double m_offset = 4 * sigma;
    double Peak = 1.0;

    std::vector<double> E_y_point1(M);
    int point1_x = N_x / 4;
    int point1_y = N_y / 2;

    std::vector<std::vector<double>> H_z(N_x, std::vector<double>(N_y, 0.0));
    std::vector<std::vector<double>> E_x(N_x + 1, std::vector<double>(N_y + 1, 0.0));
    std::vector<std::vector<double>> E_y(N_x + 1, std::vector<double>(N_y + 1, 0.0));

    for (int m = 1; m < M; ++m) {
        // Campo magnético H_z
        for (int i = 0; i < N_x; ++i) {
            for (int j = 0; j < N_y; ++j) {
                H_z[i][j] += delta_t / (mu_0 * delta_s) * (
                    E_x[i][j + 1] - E_x[i][j] +
                    E_y[i][j] - E_y[i + 1][j]);
            }
        }

        // Campo elétrico E_x
        for (int i = 0; i < N_x; ++i) {
            for (int j = 1; j < N_y; ++j) {
                E_x[i][j] += delta_t / (eps_0 * delta_s) * (
                    H_z[i][j] - H_z[i][j - 1]);
            }
        }

        // Campo elétrico E_y
        for (int i = 1; i < N_x; ++i) {
            for (int j = 0; j < N_y; ++j) {
                E_y[i][j] -= delta_t / (eps_0 * delta_s) * (
                    H_z[i][j] - H_z[i - 1][j]);
            }
        }

        // Fonte Gaussiana em E_y
        double inc = Peak * gaussder_norm((m - 1) * delta_t - (L - 1) * delta_s / c, m_offset, sigma);
        E_y[L][point1_y] += inc;

        // Salvar valor em ponto específico
        E_y_point1[m] = E_y[point1_x][point1_y];
    }

    std::cout << "Simulação finalizada. Valores de E_y armazenados." << std::endl;

    // Exportar resultados para CSV
    std::ofstream file(out_dir + "/ey_point1.csv");
    file << "tempo_ns,Ey_Vpm\n";
    for (int m = 0; m < M; ++m) {
        double time_ns = m * delta_t * 1e9;
        file << time_ns << "," << E_y_point1[m] << "\n";
    }
    file.close();
    std::cout << "Arquivo salvo em: " << out_dir << "/ey_point1.csv" << std::endl;
}

