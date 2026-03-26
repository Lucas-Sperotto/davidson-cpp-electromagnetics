// src/fdtd2d.cpp
#include "../include/fdtd2d.hpp"
#include "../include/gaussder.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace
{
using Matrix = std::vector<std::vector<double>>;

constexpr double pi = 3.14159265358979323846;
constexpr double c = 2.99792458e8;
constexpr double eps_0 = 8.854187817e-12;
constexpr double mu_0 = 4.0 * pi * 1.0e-7;
const double eta_0 = std::sqrt(mu_0 / eps_0);

Matrix make_matrix(int nx, int ny, double value = 0.0)
{
    return Matrix(nx, std::vector<double>(ny, value));
}
} // namespace

void run_fdtd2d_simulation(const Fdtd2DConfig &config)
{
    if (config.refine <= 0)
        throw std::invalid_argument("refine deve ser maior que zero.");
    if (config.pulse_compress <= 0.0)
        throw std::invalid_argument("pulse_compress deve ser maior que zero.");

    const std::filesystem::path out_dir = PROJECT_OUT_DIR;
    std::filesystem::create_directories(out_dir);

    const bool cyl_present = config.cyl_present;
    const int refine = config.refine;
    const double pulse_compress = config.pulse_compress;

    const int N_x = refine * 200;
    const int N_y = refine * 100;
    const int M = refine * 512;

    // Mantemos os parametros "MATLAB-style" para que as formulas do texto
    // possam ser conferidas diretamente contra o original.
    const int L_mat = static_cast<int>(std::round(static_cast<double>(N_x) / 2.0));
    const int L1_mat = 20;
    const int L_idx = L_mat - 1;
    const int lower_h_idx = L1_mat - 1;
    const int upper_h_idx = N_y - L1_mat - 1;
    const int right_h_idx = N_x - L1_mat - 1;
    const int right_ey_idx = N_x - L1_mat;
    const int upper_ex_idx = N_y - L1_mat;

    const double delta_s = 0.005 * 2.0 / refine;
    const double delta_t = delta_s / (c * std::sqrt(2.0));
    const double sigma = 1.0e-10 / pulse_compress;
    const double f_max = (1.0 / sigma) / (2.0 * pi);
    const double m_offset = 4.0 * sigma;
    const double peak = 1.0;

    const double radius = 0.03;
    const int N_centre_x_mat = static_cast<int>(std::round(0.75 * N_x));
    const int N_centre_y_mat = static_cast<int>(std::round(0.5 * N_y));

    const int point1_x_mat = N_x / 4;
    const int point1_y_mat = N_y / 2;
    const int point1_x_idx = point1_x_mat - 1;
    const int point1_y_idx = point1_y_mat - 1;
    const int point2_x_mat = static_cast<int>(std::round(static_cast<double>(N_x + L_mat) / 2.0));
    const int point2_y_mat = N_y / 2;
    const int point2_x_idx = point2_x_mat - 1;
    const int point2_y_idx = point2_y_mat - 1;

    if ((N_centre_x_mat - L_mat) * delta_s <= radius)
        throw std::runtime_error("Scattered/total field nao esta inteiramente a esquerda do alvo.");

    Matrix C_Ex = make_matrix(N_x + 1, N_y + 1, delta_t / (eps_0 * delta_s));
    Matrix C_Ey = make_matrix(N_x + 1, N_y + 1, delta_t / (eps_0 * delta_s));
    const double D_Hz = delta_t / (mu_0 * delta_s);

    if (cyl_present)
    {
        for (int ii = 0; ii <= N_x; ++ii)
        {
            const double ii_mat = static_cast<double>(ii + 1);
            for (int jj = 0; jj <= N_y; ++jj)
            {
                const double jj_mat = static_cast<double>(jj + 1);
                const double ex_dist = std::sqrt(
                    std::pow((ii_mat - 0.5 - (N_centre_x_mat - 1)) * delta_s, 2.0) +
                    std::pow((jj_mat - 1.0 - (N_centre_y_mat - 1)) * delta_s, 2.0));
                const double ey_dist = std::sqrt(
                    std::pow((ii_mat - 1.0 - (N_centre_x_mat - 1)) * delta_s, 2.0) +
                    std::pow((jj_mat - 0.5 - (N_centre_y_mat - 1)) * delta_s, 2.0));

                if (ex_dist <= radius)
                    C_Ex[ii][jj] = 0.0;
                if (ey_dist <= radius)
                    C_Ey[ii][jj] = 0.0;
            }
        }
    }

    std::vector<double> H_z_point1(M, 0.0);
    std::vector<double> E_y_point1(M, 0.0);
    std::vector<double> H_z_point2(M, 0.0);
    std::vector<double> E_y_point2(M, 0.0);

    Matrix H_z_nmin1 = make_matrix(N_x, N_y);
    Matrix E_x_nmin1 = make_matrix(N_x + 1, N_y + 1);
    Matrix E_y_nmin1 = make_matrix(N_x + 1, N_y + 1);
    Matrix H_z_n = make_matrix(N_x, N_y);
    Matrix E_x_n = make_matrix(N_x + 1, N_y + 1);
    Matrix E_y_n = make_matrix(N_x + 1, N_y + 1);

    std::vector<double> E_y_nmin1_inc_front(N_y + 1, 0.0);
    std::vector<double> E_y_nmin1_inc_back(N_y + 1, 0.0);
    std::vector<double> E_x_nmin1_inc_top(N_x + 1, 0.0);
    std::vector<double> E_x_nmin1_inc_bottom(N_x + 1, 0.0);
    std::vector<double> H_z_n_inc(N_y, 0.0);
    std::vector<double> H_z_n_inc2(N_x, 0.0);

    for (int m_mat = 2; m_mat <= M; ++m_mat)
    {
        const int m_idx = m_mat - 1;

        for (int i = 0; i < N_x; ++i)
        {
            for (int j = 0; j < N_y; ++j)
            {
                H_z_n[i][j] = H_z_nmin1[i][j] +
                              D_Hz * (E_x_nmin1[i][j + 1] - E_x_nmin1[i][j] +
                                      E_y_nmin1[i][j] - E_y_nmin1[i + 1][j]);
            }
        }

        const double incident_front =
            peak * gaussder_norm((m_mat - 1.0) * delta_t - (L_mat - 1.0) * delta_s / c,
                                 m_offset, sigma);
        const double incident_back =
            peak * gaussder_norm((m_mat - 1.0) * delta_t - (N_x - L1_mat) * delta_s / c,
                                 m_offset, sigma);

        for (int j = 0; j <= N_y; ++j)
        {
            E_y_nmin1_inc_front[j] = incident_front;
            E_y_nmin1_inc_back[j] = incident_back;
        }

        for (int j_mat = L1_mat; j_mat <= N_y - L1_mat; ++j_mat)
        {
            const int j = j_mat - 1;
            H_z_n[L_idx][j] =
                H_z_nmin1[L_idx][j] +
                D_Hz * (E_x_nmin1[L_idx][j + 1] - E_x_nmin1[L_idx][j] +
                        E_y_nmin1[L_idx][j] + E_y_nmin1_inc_front[j] -
                        E_y_nmin1[L_idx + 1][j]);

            H_z_n[right_h_idx][j] =
                H_z_nmin1[right_h_idx][j] +
                D_Hz * (E_x_nmin1[right_h_idx][j + 1] - E_x_nmin1[right_h_idx][j] +
                        E_y_nmin1[right_h_idx][j] - E_y_nmin1_inc_back[j] -
                        E_y_nmin1[right_h_idx + 1][j]);
        }

        for (int i = 0; i <= N_x; ++i)
        {
            E_x_nmin1_inc_top[i] = 0.0;
            E_x_nmin1_inc_bottom[i] = 0.0;
        }

        for (int i_mat = L_mat; i_mat <= N_x - L1_mat; ++i_mat)
        {
            const int i = i_mat - 1;
            H_z_n[i][upper_h_idx] =
                H_z_nmin1[i][upper_h_idx] +
                D_Hz * (E_x_nmin1[i][upper_h_idx + 1] + E_x_nmin1_inc_top[i] -
                        E_x_nmin1[i][upper_h_idx] + E_y_nmin1[i][upper_h_idx] -
                        E_y_nmin1[i + 1][upper_h_idx]);

            H_z_n[i][lower_h_idx] =
                H_z_nmin1[i][lower_h_idx] +
                D_Hz * (E_x_nmin1[i][lower_h_idx + 1] - E_x_nmin1[i][lower_h_idx] -
                        E_x_nmin1_inc_bottom[i] + E_y_nmin1[i][lower_h_idx] -
                        E_y_nmin1[i + 1][lower_h_idx]);
        }

        H_z_n[L_idx][lower_h_idx] =
            H_z_nmin1[L_idx][lower_h_idx] +
            D_Hz * (E_x_nmin1[L_idx][lower_h_idx + 1] - E_x_nmin1_inc_bottom[L_idx] -
                    E_x_nmin1[L_idx][lower_h_idx] + E_y_nmin1[L_idx][lower_h_idx] +
                    E_y_nmin1_inc_front[lower_h_idx] - E_y_nmin1[L_idx + 1][lower_h_idx]);

        H_z_n[L_idx][upper_h_idx] =
            H_z_nmin1[L_idx][upper_h_idx] +
            D_Hz * (E_x_nmin1[L_idx][upper_h_idx + 1] + E_x_nmin1_inc_top[L_idx] -
                    E_x_nmin1[L_idx][upper_h_idx] + E_y_nmin1[L_idx][upper_h_idx] +
                    E_y_nmin1_inc_front[upper_h_idx] - E_y_nmin1[L_idx + 1][upper_h_idx]);

        for (int i = 0; i < N_x; ++i)
        {
            for (int j = 1; j < N_y; ++j)
            {
                E_x_n[i][j] =
                    E_x_nmin1[i][j] + C_Ex[i][j] * (H_z_n[i][j] - H_z_n[i][j - 1]);
            }
        }

        for (int i = 1; i < N_x; ++i)
        {
            for (int j = 0; j < N_y; ++j)
            {
                E_y_n[i][j] =
                    E_y_nmin1[i][j] - C_Ey[i][j] * (H_z_n[i][j] - H_z_n[i - 1][j]);
            }
        }

        const double incident_h_left =
            (peak / eta_0) *
            gaussder_norm((m_mat - 0.5) * delta_t - (L_mat - 0.5) * delta_s / c,
                          m_offset, sigma);
        const double incident_h_right =
            (peak / eta_0) *
            gaussder_norm((m_mat - 0.5) * delta_t - (N_x - L1_mat - 0.5) * delta_s / c,
                          m_offset, sigma);

        for (int j = 0; j < N_y; ++j)
            H_z_n_inc[j] = incident_h_left;

        for (int j_mat = L1_mat; j_mat <= N_y - L1_mat; ++j_mat)
        {
            const int j = j_mat - 1;
            E_y_n[L_idx][j] =
                E_y_nmin1[L_idx][j] -
                C_Ey[L_idx][j] * (H_z_n[L_idx][j] - H_z_n_inc[j] - H_z_n[L_idx - 1][j]);
        }

        for (int j = 0; j < N_y; ++j)
            H_z_n_inc[j] = incident_h_right;

        for (int j_mat = L1_mat; j_mat <= N_y - L1_mat; ++j_mat)
        {
            const int j = j_mat - 1;
            E_y_n[right_ey_idx][j] =
                E_y_nmin1[right_ey_idx][j] -
                C_Ey[right_ey_idx][j] *
                    (H_z_n[right_ey_idx][j] + H_z_n_inc[j] - H_z_n[right_h_idx][j]);
        }

        for (int ii_mat = 1; ii_mat <= N_x; ++ii_mat)
        {
            const int i = ii_mat - 1;
            H_z_n_inc2[i] =
                (peak / eta_0) *
                gaussder_norm((m_mat - 0.5) * delta_t - (ii_mat - 0.5) * delta_s / c,
                              m_offset, sigma);
        }

        for (int i_mat = L_mat; i_mat <= N_x - L1_mat; ++i_mat)
        {
            const int i = i_mat - 1;
            E_x_n[i][upper_ex_idx] =
                E_x_nmin1[i][upper_ex_idx] +
                C_Ex[i][upper_ex_idx] *
                    (H_z_n[i][upper_ex_idx] - H_z_n[i][upper_h_idx] + H_z_n_inc2[i]);

            E_x_n[i][lower_h_idx] =
                E_x_nmin1[i][lower_h_idx] +
                C_Ex[i][lower_h_idx] *
                    (H_z_n[i][lower_h_idx] - H_z_n_inc2[i] - H_z_n[i][lower_h_idx - 1]);
        }

        for (int j = 0; j <= N_y; ++j)
        {
            E_y_n[0][j] = E_y_nmin1[0][j] * (1.0 - c * delta_t / delta_s) +
                          (c * delta_t / delta_s) * E_y_nmin1[1][j];
            E_y_n[N_x][j] = E_y_nmin1[N_x][j] * (1.0 - c * delta_t / delta_s) +
                            (c * delta_t / delta_s) * E_y_nmin1[N_x - 1][j];
        }

        for (int i = 0; i <= N_x; ++i)
        {
            E_x_n[i][0] = E_x_nmin1[i][0] * (1.0 - c * delta_t / delta_s) +
                          (c * delta_t / delta_s) * E_x_nmin1[i][1];
            E_x_n[i][N_y] = E_x_nmin1[i][N_y] * (1.0 - c * delta_t / delta_s) +
                            (c * delta_t / delta_s) * E_x_nmin1[i][N_y - 1];
        }

        H_z_point1[m_idx] = H_z_n[point1_x_idx][point1_y_idx];
        H_z_point2[m_idx] = H_z_n[point2_x_idx][point2_y_idx];
        E_y_point1[m_idx] = E_y_n[point1_x_idx][point1_y_idx];
        E_y_point2[m_idx] = E_y_n[point2_x_idx][point2_y_idx];

        H_z_nmin1 = H_z_n;
        E_x_nmin1 = E_x_n;
        E_y_nmin1 = E_y_n;
    }

    std::ofstream file(out_dir / config.output_filename);
    file.precision(15);
    file << std::scientific;
    file << "tempo_ns,Ey_Vpm\n";
    for (int m = 0; m < M; ++m)
    {
        const double time_ns = m * delta_t * 1.0e9;
        file << time_ns << "," << E_y_point1[m] << "\n";
    }
    file.close();

    std::ofstream meta(out_dir / config.metadata_filename);
    meta << "key,value\n";
    meta << "cyl_present," << (cyl_present ? 1 : 0) << "\n";
    meta << "refine," << refine << "\n";
    meta << "pulse_compress," << pulse_compress << "\n";
    meta << "N_x," << N_x << "\n";
    meta << "N_y," << N_y << "\n";
    meta << "M," << M << "\n";
    meta << "L," << L_mat << "\n";
    meta << "L1," << L1_mat << "\n";
    meta << "delta_s," << delta_s << "\n";
    meta << "delta_t," << delta_t << "\n";
    meta << "sigma," << sigma << "\n";
    meta << "f_max," << f_max << "\n";
    meta << "m_offset," << m_offset << "\n";
    meta << "radius," << radius << "\n";
    meta << "point1_x_mat," << point1_x_mat << "\n";
    meta << "point1_y_mat," << point1_y_mat << "\n";
    meta << "point2_x_mat," << point2_x_mat << "\n";
    meta << "point2_y_mat," << point2_y_mat << "\n";
    meta.close();

    std::cout << "Simulacao 2D com ABC finalizada. Resultados em "
              << (out_dir / config.output_filename) << std::endl;
}
