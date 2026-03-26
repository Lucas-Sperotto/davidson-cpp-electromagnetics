#include "include/chapter7_functions.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <vector>

namespace fs = std::filesystem;

int main()
{
    // Translation of V_pot_height.m:
    // scalar potential versus distance for different normalized heights.
    const fs::path out_dir = PROJECT_OUT_DIR;
    fs::create_directories(out_dir);

    const double c = 2.997925e8;
    const double freq = 10e9;
    const double eps_r_prime = 10.0;
    const double tan_delta = 0.0;
    const std::complex<double> eps_r =
        eps_r_prime * (1.0 - std::complex<double>(0.0, tan_delta));
    const double k_0 = 2.0 * std::numbers::pi * freq / c;
    const std::complex<double> k = k_0 * std::sqrt(eps_r);
    const double rho_max = 10.0 / std::abs(k);
    const double rho_min = 0.01 / std::abs(k);
    const int N_rho = 100;
    const double delta_rho = std::exp(std::log(rho_max / rho_min) / N_rho);
    const std::vector<double> b_values = {0.8, 0.5, 0.2};

    std::vector<double> rho_values(N_rho, 0.0);
    std::vector<std::vector<std::complex<double>>> potentials(
        b_values.size(), std::vector<std::complex<double>>(N_rho, 0.0));

    for (std::size_t jj = 0; jj < b_values.size(); ++jj)
    {
        SommerfeldContext context;
        context.k_0 = k_0;
        context.k = k;
        context.eps_r = eps_r;
        context.eps_r_prime = eps_r_prime;
        context.d = b_values[jj] * std::numbers::pi /
                    (2.0 * k_0 * std::sqrt(eps_r_prime - 1.0));
        context.lambda_p =
            root_D_TM(k_0, context.k, context.eps_r, context.d, 50);

        double rho = rho_min;
        for (int kk = 0; kk < N_rho; ++kk)
        {
            rho_values[kk] = rho;
            context.rho = rho;
            potentials[jj][kk] = V_int(context);
            rho *= delta_rho;
        }
    }

    const double norm = std::abs(potentials.front().front());

    std::ofstream out_mag(out_dir / "v_pot_height_mag.csv");
    out_mag << "b,rho_m,k_rho,abs_v,normalized_abs_v,re_v,im_v\n";
    out_mag << std::scientific << std::setprecision(10);
    for (std::size_t jj = 0; jj < b_values.size(); ++jj)
    {
        for (int kk = 0; kk < N_rho; ++kk)
        {
            out_mag << b_values[jj] << "," << rho_values[kk] << ","
                    << rho_values[kk] * std::abs(k) << ","
                    << std::abs(potentials[jj][kk]) << ","
                    << std::abs(potentials[jj][kk]) / norm << ","
                    << potentials[jj][kk].real() << ","
                    << potentials[jj][kk].imag() << "\n";
        }
    }

    std::ofstream out_phase(out_dir / "v_pot_height_phase.csv");
    out_phase << "b,rho_m,k_rho,phase_deg\n";
    out_phase << std::scientific << std::setprecision(10);
    for (std::size_t jj = 0; jj < b_values.size(); ++jj)
    {
        for (int kk = 0; kk < N_rho; ++kk)
        {
            out_phase << b_values[jj] << "," << rho_values[kk] << ","
                      << rho_values[kk] * std::abs(k) << ","
                      << std::arg(potentials[jj][kk]) * 180.0 / std::numbers::pi << "\n";
        }
    }

    std::ofstream summary(out_dir / "v_pot_height_summary.csv");
    summary << "key,value\n";
    summary << "freq_hz," << freq << "\n";
    summary << "eps_r_prime," << eps_r_prime << "\n";
    summary << "tan_delta," << tan_delta << "\n";
    summary << "k_0," << k_0 << "\n";
    summary << "k_abs," << std::abs(k) << "\n";
    summary << "rho_min," << rho_min << "\n";
    summary << "rho_max," << rho_max << "\n";
    summary << "N_rho," << N_rho << "\n";
    summary << "norm," << norm << "\n";

    std::cout << "Arquivos gerados: v_pot_height_mag.csv, v_pot_height_phase.csv e v_pot_height_summary.csv\n";
    return 0;
}
