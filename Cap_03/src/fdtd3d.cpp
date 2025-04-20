// src/fdtd3d.cpp
#include "../include/fdtd3d.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <complex>

const double eps_0 = 8.854e-12;
const double mu_0 = 4 * M_PI * 1e-7;
const double c = 1.0 / std::sqrt(eps_0 * mu_0);

void run_fdtd3d_simulation() {
    int refine = 2;
    int N_x = 8 * refine;
    int N_y = 4 * refine;
    int N_z = 6 * refine;

    double L_x = 1.0;
    double L_y = 0.5;
    double L_z = 0.75;

    double dx = L_x / N_x;
    double dy = L_y / N_y;
    double dz = L_z / N_z;

    double delta_s = std::min({dx, dy, dz});
    double delta_t = delta_s / (std::sqrt(3.0) * c);

    int N_t = 1024 * refine;

    std::vector<std::vector<std::vector<double>>> E_x(N_x, std::vector<std::vector<double>>(N_y + 1, std::vector<double>(N_z + 1, 0.0)));
    std::vector<std::vector<std::vector<double>>> E_y(N_x + 1, std::vector<std::vector<double>>(N_y, std::vector<double>(N_z + 1, 0.0)));
    std::vector<std::vector<std::vector<double>>> E_z(N_x + 1, std::vector<std::vector<double>>(N_y + 1, std::vector<double>(N_z, 0.0)));
    std::vector<std::vector<std::vector<double>>> H_x(N_x + 1, std::vector<std::vector<double>>(N_y, std::vector<double>(N_z, 0.0)));
    std::vector<std::vector<std::vector<double>>> H_y(N_x, std::vector<std::vector<double>>(N_y + 1, std::vector<double>(N_z, 0.0)));
    std::vector<std::vector<std::vector<double>>> H_z(N_x, std::vector<std::vector<double>>(N_y, std::vector<double>(N_z + 1, 0.0)));

    for (int i = 0; i < N_x; ++i)
        for (int j = 0; j < N_y; ++j)
            for (int k = 0; k < N_z + 1; ++k)
                H_z[i][j][k] = ((double) rand() / RAND_MAX) - 0.5;

    std::vector<double> t(N_t);
    std::vector<double> H_z_t(N_t);

    for (int n = 0; n < N_t; ++n) {
        t[n] = n * delta_t;
        H_z_t[n] = H_z[N_x / 2][N_y / 2][N_z / 2];

        for (int i = 0; i < N_x + 1; ++i)
            for (int j = 0; j < N_y; ++j)
                for (int k = 0; k < N_z; ++k)
                    H_x[i][j][k] += delta_t / mu_0 * (
                        (j < N_y - 1 && k < N_z - 1 ? E_y[i][j][k + 1] - E_y[i][j][k] : 0.0) / dz -
                        (j < N_y - 1 && k < N_z - 1 ? E_z[i][j + 1][k] - E_z[i][j][k] : 0.0) / dy);

        for (int i = 0; i < N_x; ++i)
            for (int j = 0; j < N_y + 1; ++j)
                for (int k = 0; k < N_z; ++k)
                    H_y[i][j][k] += delta_t / mu_0 * (
                        (i < N_x - 1 && k < N_z - 1 ? E_z[i + 1][j][k] - E_z[i][j][k] : 0.0) / dx -
                        (i < N_x - 1 && k < N_z - 1 ? E_x[i][j][k + 1] - E_x[i][j][k] : 0.0) / dz);

        for (int i = 0; i < N_x; ++i)
            for (int j = 0; j < N_y; ++j)
                for (int k = 0; k < N_z + 1; ++k)
                    H_z[i][j][k] += delta_t / mu_0 * (
                        (i < N_x - 1 && j < N_y - 1 ? E_x[i][j + 1][k] - E_x[i][j][k] : 0.0) / dy -
                        (i < N_x - 1 && j < N_y - 1 ? E_y[i + 1][j][k] - E_y[i][j][k] : 0.0) / dx);

        for (int i = 0; i < N_x; ++i)
            for (int j = 1; j < N_y; ++j)
                for (int k = 1; k < N_z; ++k)
                    E_x[i][j][k] += delta_t / eps_0 * (
                        (k < N_z ? H_z[i][j][k] - H_z[i][j - 1][k] : 0.0) / dy -
                        (j < N_y ? H_y[i][j][k] - H_y[i][j][k - 1] : 0.0) / dz);

        for (int i = 1; i < N_x; ++i)
            for (int j = 0; j < N_y; ++j)
                for (int k = 1; k < N_z; ++k)
                    E_y[i][j][k] += delta_t / eps_0 * (
                        (i < N_x ? H_x[i][j][k] - H_x[i][j][k - 1] : 0.0) / dz -
                        (k < N_z ? H_z[i][j][k] - H_z[i - 1][j][k] : 0.0) / dx);

        for (int i = 1; i < N_x; ++i)
            for (int j = 1; j < N_y; ++j)
                for (int k = 0; k < N_z; ++k)
                    E_z[i][j][k] += delta_t / eps_0 * (
                        (j < N_y ? H_y[i][j][k] - H_y[i - 1][j][k] : 0.0) / dx -
                        (i < N_x ? H_x[i][j][k] - H_x[i][j - 1][k] : 0.0) / dy);
    }

    std::ofstream file("out/hz_center_fft.csv");
    file << "freq_Hz,abs_Hz\n";
    double delta_f = 1.0 / (N_t * delta_t);
    for (int k = 1; k < N_t / 2; ++k) {
        std::complex<double> sum(0.0, 0.0);
        for (int n = 0; n < N_t; ++n) {
            double angle = -2.0 * M_PI * k * n / N_t;
            sum += H_z_t[n] * std::exp(std::complex<double>(0, angle));
        }
        file << k * delta_f << "," << std::abs(sum) << "\n";
    }
    file.close();
    std::cout << "FFT de H_z salva em: out/hz_center_fft.csv" << std::endl;
}
