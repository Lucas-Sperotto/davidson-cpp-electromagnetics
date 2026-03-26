// src/fdtd2d_pml.cpp
#include "../include/fdtd2d_pml.hpp"
#include "../include/gaussder.hpp"

#include <algorithm>
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
constexpr double c = 2.997925e8;
constexpr double eps_0 = 8.854e-12;
constexpr double mu_0 = 4.0 * pi * 1.0e-7;
const double eta_0 = std::sqrt(mu_0 / eps_0);

Matrix make_matrix(int nx, int ny, double value = 0.0)
{
    return Matrix(nx, std::vector<double>(ny, value));
}

void write_matrix_csv(const std::filesystem::path &path, const Matrix &matrix)
{
    std::ofstream file(path);
    file.precision(17);
    file << std::scientific;
    for (std::size_t i = 0; i < matrix.size(); ++i)
    {
        for (std::size_t j = 0; j < matrix[i].size(); ++j)
        {
            if (j)
                file << ',';
            file << matrix[i][j];
        }
        file << '\n';
    }
}

void write_vector_csv(const std::filesystem::path &path, const std::vector<double> &values)
{
    std::ofstream file(path);
    file.precision(17);
    file << std::scientific;
    for (double value : values)
        file << value << '\n';
}
} // namespace

void run_fdtd2d_pml_simulation(const Fdtd2DPmlConfig &config)
{
    if (config.refine <= 0)
        throw std::invalid_argument("refine deve ser maior que zero.");
    if (config.pulse_compress <= 0.0)
        throw std::invalid_argument("pulse_compress deve ser maior que zero.");
    if (config.d_cell <= 0)
        throw std::invalid_argument("d_cell deve ser maior que zero.");
    if (config.poly_m < 0)
        throw std::invalid_argument("poly_m deve ser nao negativo.");

    const std::filesystem::path out_dir = PROJECT_OUT_DIR;
    std::filesystem::create_directories(out_dir);

    const bool cyl_present = config.cyl_present;
    const bool line_source = cyl_present ? false : config.line_source;
    const int refine = config.refine;
    const int d_cell = config.d_cell;
    const int poly_m = config.poly_m;
    const bool snapshot_requested = config.snapshot_step > 0;
    const std::string snapshot_prefix =
        config.snapshot_prefix.empty() ? "fdtd2d_pml_snapshot" : config.snapshot_prefix;

    const int N_x = refine * 200;
    const int N_y = refine * 200;
    const int M = refine * 350;
    const int L_mat = static_cast<int>(std::round(static_cast<double>(N_x) / 2.0));
    const int L_idx = L_mat - 1;

    const double delta_s = line_source ? 0.005 : 0.005 / refine;
    const double delta_t = delta_s / (c * std::sqrt(2.0));
    const double sigma = 1.0e-10 / config.pulse_compress;
    const double f_max = (1.0 / sigma) / (2.0 * pi);
    const double m_offset = 4.0 * sigma;
    const double peak = 1.0;

    const double radius = 0.03;
    const int N_centre_x_mat = static_cast<int>(std::round(0.75 * N_x));
    const int N_centre_y_mat = static_cast<int>(std::round(0.5 * N_y));

    if ((N_centre_x_mat - L_mat) * delta_s <= radius)
        throw std::runtime_error("Scattered/total field nao esta inteiramente a esquerda do alvo.");

    Matrix sigma_x = make_matrix(N_x, N_y);
    Matrix sigma_y = make_matrix(N_x, N_y);
    Matrix sigma_ast_x = make_matrix(N_x, N_y);
    Matrix sigma_ast_y = make_matrix(N_x, N_y);

    const double sigma_max = 0.8 * (poly_m + 1) / (eta_0 * delta_s);

    for (int jj_mat = 1; jj_mat <= d_cell; ++jj_mat)
    {
        const double dist = d_cell - (jj_mat - 1.0);
        const double value = sigma_max * std::pow(dist / d_cell, poly_m);
        const int low_j = jj_mat - 1;
        const int high_j = N_y - jj_mat;
        for (int i = 0; i < N_x; ++i)
        {
            sigma_y[i][low_j] = value;
            sigma_y[i][high_j] = value;
        }
    }

    for (int jj_mat = 1; jj_mat <= d_cell; ++jj_mat)
    {
        const double dist = d_cell - (jj_mat - 1.0) - 0.5;
        const double value = eta_0 * eta_0 * sigma_max * std::pow(dist / d_cell, poly_m);
        const int low_j = jj_mat - 1;
        const int high_j = N_y - jj_mat - 1;
        for (int i = 0; i < N_x; ++i)
        {
            sigma_ast_y[i][low_j] = value;
            if (high_j >= 0)
                sigma_ast_y[i][high_j] = value;
        }
    }

    for (int ii_mat = 1; ii_mat <= d_cell; ++ii_mat)
    {
        const double dist = d_cell - (ii_mat - 1.0);
        const double value = sigma_max * std::pow(dist / d_cell, poly_m);
        const int low_i = ii_mat - 1;
        const int high_i = N_x - ii_mat;
        for (int j = 0; j < N_y; ++j)
        {
            sigma_x[low_i][j] = value;
            sigma_x[high_i][j] = value;
        }
    }

    for (int ii_mat = 1; ii_mat <= d_cell; ++ii_mat)
    {
        const double dist = d_cell - (ii_mat - 1.0) - 0.5;
        const double value = eta_0 * eta_0 * sigma_max * std::pow(dist / d_cell, poly_m);
        const int low_i = ii_mat - 1;
        const int high_i = N_x - ii_mat - 1;
        for (int j = 0; j < N_y; ++j)
        {
            sigma_ast_x[low_i][j] = value;
            if (high_i >= 0)
                sigma_ast_x[high_i][j] = value;
        }
    }

    Matrix C_aEx = make_matrix(N_x, N_y);
    Matrix C_bEx = make_matrix(N_x, N_y);
    Matrix C_aEy = make_matrix(N_x, N_y);
    Matrix C_bEy = make_matrix(N_x, N_y);
    Matrix D_aHzx = make_matrix(N_x, N_y);
    Matrix D_bHzx = make_matrix(N_x, N_y);
    Matrix D_aHzy = make_matrix(N_x, N_y);
    Matrix D_bHzy = make_matrix(N_x, N_y);

    for (int i = 0; i < N_x; ++i)
    {
        for (int j = 0; j < N_y; ++j)
        {
            C_aEx[i][j] = (1.0 - sigma_y[i][j] * delta_t / (2.0 * eps_0)) /
                          (1.0 + sigma_y[i][j] * delta_t / (2.0 * eps_0));
            C_bEx[i][j] = delta_t / (eps_0 * delta_s) /
                          (1.0 + sigma_y[i][j] * delta_t / (2.0 * eps_0));
            C_aEy[i][j] = (1.0 - sigma_x[i][j] * delta_t / (2.0 * eps_0)) /
                          (1.0 + sigma_x[i][j] * delta_t / (2.0 * eps_0));
            C_bEy[i][j] = delta_t / (eps_0 * delta_s) /
                          (1.0 + sigma_x[i][j] * delta_t / (2.0 * eps_0));
            D_aHzx[i][j] = (1.0 - sigma_ast_x[i][j] * delta_t / (2.0 * mu_0)) /
                           (1.0 + sigma_ast_x[i][j] * delta_t / (2.0 * mu_0));
            D_bHzx[i][j] = delta_t / (mu_0 * delta_s) /
                           (1.0 + sigma_ast_x[i][j] * delta_t / (2.0 * mu_0));
            D_aHzy[i][j] = (1.0 - sigma_ast_y[i][j] * delta_t / (2.0 * mu_0)) /
                           (1.0 + sigma_ast_y[i][j] * delta_t / (2.0 * mu_0));
            D_bHzy[i][j] = delta_t / (mu_0 * delta_s) /
                           (1.0 + sigma_ast_y[i][j] * delta_t / (2.0 * mu_0));
        }
    }

    if (cyl_present)
    {
        for (int ii = 0; ii < N_x; ++ii)
        {
            const double ii_mat = static_cast<double>(ii + 1);
            for (int jj = 0; jj < N_y; ++jj)
            {
                const double jj_mat = static_cast<double>(jj + 1);
                const double ex_dist = std::sqrt(
                    std::pow((ii_mat - 0.5 - (N_centre_x_mat - 1)) * delta_s, 2.0) +
                    std::pow((jj_mat - 1.0 - (N_centre_y_mat - 1)) * delta_s, 2.0));
                const double ey_dist = std::sqrt(
                    std::pow((ii_mat - 1.0 - (N_centre_x_mat - 1)) * delta_s, 2.0) +
                    std::pow((jj_mat - 0.5 - (N_centre_y_mat - 1)) * delta_s, 2.0));

                if (ex_dist <= radius)
                    C_bEx[ii][jj] = 0.0;
                if (ey_dist <= radius)
                    C_bEy[ii][jj] = 0.0;
            }
        }
    }

    int point1_x_mat = 190;
    if (refine == 2)
        point1_x_mat = 290;
    else if (refine != 1)
        point1_x_mat = std::min(N_x - d_cell, std::max(1, static_cast<int>(std::round(0.95 * N_x))));

    const int point1_y_mat = static_cast<int>(std::round(static_cast<double>(N_y) / 2.0));
    const int point2_x_mat = static_cast<int>(std::round(static_cast<double>(N_x + L_mat) / 2.0));
    const int point2_y_mat = static_cast<int>(std::round(static_cast<double>(N_y) / 2.0));
    const int point1_x_idx = point1_x_mat - 1;
    const int point1_y_idx = point1_y_mat - 1;
    const int point2_x_idx = point2_x_mat - 1;
    const int point2_y_idx = point2_y_mat - 1;
    const int source_x_idx = (N_x / 2) - 1;
    const int source_y_idx = (N_y / 2) - 1;

    std::vector<double> H_z_point1(M, 0.0);
    std::vector<double> H_z_point2(M, 0.0);
    std::vector<double> E_y_point1(M, 0.0);
    std::vector<double> E_y_point2(M, 0.0);

    Matrix H_zx_nmin1 = make_matrix(N_x, N_y);
    Matrix H_zy_nmin1 = make_matrix(N_x, N_y);
    Matrix E_x_nmin1 = make_matrix(N_x, N_y);
    Matrix E_y_nmin1 = make_matrix(N_x, N_y);

    Matrix H_zx_n = make_matrix(N_x, N_y);
    Matrix H_zy_n = make_matrix(N_x, N_y);
    Matrix H_z_n = make_matrix(N_x, N_y);
    Matrix E_x_n = make_matrix(N_x, N_y);
    Matrix E_y_n = make_matrix(N_x, N_y);
    std::vector<double> E_y_nmin1_inc_profile(N_y, 0.0);
    std::vector<double> H_z_n_inc_profile(N_y, 0.0);

    for (int m_mat = 2; m_mat <= M; ++m_mat)
    {
        const int m_idx = m_mat - 1;

        for (int i = 0; i < N_x - 1; ++i)
        {
            for (int j = 0; j < N_y - 1; ++j)
            {
                H_zx_n[i][j] = D_aHzx[i][j] * H_zx_nmin1[i][j] -
                               D_bHzx[i][j] * (E_y_nmin1[i + 1][j] - E_y_nmin1[i][j]);
                H_zy_n[i][j] = D_aHzy[i][j] * H_zy_nmin1[i][j] +
                               D_bHzy[i][j] * (E_x_nmin1[i][j + 1] - E_x_nmin1[i][j]);
            }
        }

        if (line_source)
        {
            for (int i = 0; i < N_x; ++i)
            {
                for (int j = 0; j < N_y; ++j)
                    H_z_n[i][j] = H_zx_n[i][j] + H_zy_n[i][j];
            }
            H_z_n[source_x_idx][source_y_idx] = gaussder_norm((m_mat - 1.0) * delta_t, m_offset, sigma);
        }
        else
        {
            const double incident_ey =
                peak * gaussder_norm((m_mat - 1.0) * delta_t - (L_mat - 1.0) * delta_s / c,
                                     m_offset, sigma);
            std::fill(E_y_nmin1_inc_profile.begin(), E_y_nmin1_inc_profile.end(), incident_ey);

            for (int j = 0; j < N_y; ++j)
            {
                H_zx_n[L_idx][j] = 0.0;
                H_zy_n[L_idx][j] = 0.0;
            }

            for (int j = 0; j < N_y - 1; ++j)
            {
                H_zx_n[L_idx][j] =
                    D_aHzx[L_idx][j] * H_zx_nmin1[L_idx][j] -
                    D_bHzx[L_idx][j] * (E_x_nmin1[L_idx][j + 1] - E_x_nmin1[L_idx][j]);

                H_zy_n[L_idx][j] =
                    D_aHzy[L_idx][j] * H_zy_nmin1[L_idx][j] +
                    D_bHzy[L_idx][j] *
                        (E_y_nmin1[L_idx + 1][j] - E_y_nmin1[L_idx][j] - incident_ey);
            }

            for (int i = 0; i < N_x; ++i)
            {
                for (int j = 0; j < N_y; ++j)
                    H_z_n[i][j] = H_zx_n[i][j] + H_zy_n[i][j];
            }
        }

        for (int i = 1; i < N_x; ++i)
        {
            for (int j = 1; j < N_y; ++j)
            {
                E_x_n[i][j] = C_aEx[i][j] * E_x_nmin1[i][j] +
                              C_bEx[i][j] * (H_z_n[i][j] - H_z_n[i][j - 1]);
                E_y_n[i][j] = C_aEy[i][j] * E_y_nmin1[i][j] -
                              C_bEy[i][j] * (H_z_n[i][j] - H_z_n[i - 1][j]);
            }
        }

        if (!line_source)
        {
            const double incident_h =
                (peak / eta_0) *
                gaussder_norm((m_mat - 0.5) * delta_t - (L_mat - 0.5) * delta_s / c,
                              m_offset, sigma);
            std::fill(H_z_n_inc_profile.begin(), H_z_n_inc_profile.end(), incident_h);

            for (int j = 1; j < N_y; ++j)
            {
                E_y_n[L_idx][j] =
                    C_aEy[L_idx][j] * E_y_nmin1[L_idx][j] -
                    C_bEy[L_idx][j] * (H_z_n[L_idx][j] - incident_h - H_z_n[L_idx - 1][j]);
            }
        }
        else
        {
            std::fill(E_y_nmin1_inc_profile.begin(), E_y_nmin1_inc_profile.end(), 0.0);
            std::fill(H_z_n_inc_profile.begin(), H_z_n_inc_profile.end(), 0.0);
        }

        for (int j = 0; j < N_y; ++j)
        {
            E_y_n[0][j] = 0.0;
            E_y_n[N_x - 1][j] = 0.0;
        }
        for (int i = 0; i < N_x; ++i)
        {
            E_x_n[i][0] = 0.0;
            E_x_n[i][N_y - 1] = 0.0;
        }

        H_z_point1[m_idx] = H_z_n[point1_x_idx][point1_y_idx];
        H_z_point2[m_idx] = H_z_n[point2_x_idx][point2_y_idx];
        E_y_point1[m_idx] = E_y_n[point1_x_idx][point1_y_idx];
        E_y_point2[m_idx] = E_y_n[point2_x_idx][point2_y_idx];

        if (snapshot_requested && m_mat == config.snapshot_step)
        {
            write_matrix_csv(out_dir / (snapshot_prefix + "_H_zx_n.csv"), H_zx_n);
            write_matrix_csv(out_dir / (snapshot_prefix + "_H_zy_n.csv"), H_zy_n);
            write_matrix_csv(out_dir / (snapshot_prefix + "_H_z_n.csv"), H_z_n);
            write_matrix_csv(out_dir / (snapshot_prefix + "_E_x_n.csv"), E_x_n);
            write_matrix_csv(out_dir / (snapshot_prefix + "_E_y_n.csv"), E_y_n);
            write_matrix_csv(out_dir / (snapshot_prefix + "_H_zx_nmin1.csv"), H_zx_nmin1);
            write_matrix_csv(out_dir / (snapshot_prefix + "_H_zy_nmin1.csv"), H_zy_nmin1);
            write_matrix_csv(out_dir / (snapshot_prefix + "_E_x_nmin1.csv"), E_x_nmin1);
            write_matrix_csv(out_dir / (snapshot_prefix + "_E_y_nmin1.csv"), E_y_nmin1);
            write_vector_csv(out_dir / (snapshot_prefix + "_E_y_nmin1_inc.csv"), E_y_nmin1_inc_profile);
            write_vector_csv(out_dir / (snapshot_prefix + "_H_z_n_inc.csv"), H_z_n_inc_profile);
        }

        H_zx_nmin1 = H_zx_n;
        H_zy_nmin1 = H_zy_n;
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
    meta << "line_source," << (line_source ? 1 : 0) << "\n";
    meta << "refine," << refine << "\n";
    meta << "pulse_compress," << config.pulse_compress << "\n";
    meta << "d_cell," << d_cell << "\n";
    meta << "poly_m," << poly_m << "\n";
    meta << "N_x," << N_x << "\n";
    meta << "N_y," << N_y << "\n";
    meta << "M," << M << "\n";
    meta << "L," << L_mat << "\n";
    meta << "delta_s," << delta_s << "\n";
    meta << "delta_t," << delta_t << "\n";
    meta << "sigma," << sigma << "\n";
    meta << "sigma_max," << sigma_max << "\n";
    meta << "f_max," << f_max << "\n";
    meta << "m_offset," << m_offset << "\n";
    meta << "radius," << radius << "\n";
    meta << "point1_x_mat," << point1_x_mat << "\n";
    meta << "point1_y_mat," << point1_y_mat << "\n";
    meta << "point2_x_mat," << point2_x_mat << "\n";
    meta << "point2_y_mat," << point2_y_mat << "\n";
    meta << "snapshot_step," << config.snapshot_step << "\n";
    meta.close();

    std::cout << "Simulacao 2D com PML finalizada. Resultados em "
              << (out_dir / config.output_filename) << std::endl;
}
