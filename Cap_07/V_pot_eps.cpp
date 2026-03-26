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
    // Translation of V_pot_eps.m:
    // scalar potential versus distance for different relative permittivities.
    const fs::path out_dir = PROJECT_OUT_DIR;
    fs::create_directories(out_dir);

    const double c = 2.997925e8;
    const double freq = 10e9;
    const double lambda_0 = c / freq;
    const double tan_delta = 0.0;
    const double k_0 = 2.0 * std::numbers::pi * freq / c;
    const double d = 0.05 * lambda_0;
    const double rho_max = 10.0 / k_0;
    const double rho_min = 0.01 / k_0;
    const int N_rho = 50;
    const double delta_rho = std::exp(std::log(rho_max / rho_min) / N_rho);
    const std::vector<double> rel_eps = {1.01, 2.2, 4.34, 9.6};

    std::vector<double> rho_values(N_rho, 0.0);
    std::vector<std::vector<std::complex<double>>> potentials(
        rel_eps.size(), std::vector<std::complex<double>>(N_rho, 0.0));

    for (std::size_t jj = 0; jj < rel_eps.size(); ++jj)
    {
        SommerfeldContext context;
        context.k_0 = k_0;
        context.eps_r_prime = rel_eps[jj];
        context.eps_r =
            rel_eps[jj] * (1.0 - std::complex<double>(0.0, tan_delta));
        context.k = k_0 * std::sqrt(context.eps_r);
        context.d = d;
        context.lambda_p =
            root_D_TM(k_0, context.k, context.eps_r, d, 50);

        double rho = rho_min;
        for (int kk = 0; kk < N_rho; ++kk)
        {
            rho_values[kk] = rho;
            context.rho = rho;
            potentials[jj][kk] = V_int(context);
            rho *= delta_rho;
        }
    }

    const double norm = std::abs(potentials.back().front());

    std::ofstream out(out_dir / "v_pot_eps.csv");
    out << "eps_r_prime,rho_m,k0_rho,abs_v,normalized_abs_v,re_v,im_v\n";
    out << std::scientific << std::setprecision(10);
    for (std::size_t jj = 0; jj < rel_eps.size(); ++jj)
    {
        for (int kk = 0; kk < N_rho; ++kk)
        {
            out << rel_eps[jj] << "," << rho_values[kk] << ","
                << rho_values[kk] * k_0 << ","
                << std::abs(potentials[jj][kk]) << ","
                << std::abs(potentials[jj][kk]) / norm << ","
                << potentials[jj][kk].real() << ","
                << potentials[jj][kk].imag() << "\n";
        }
    }

    std::ofstream summary(out_dir / "v_pot_eps_summary.csv");
    summary << "key,value\n";
    summary << "freq_hz," << freq << "\n";
    summary << "lambda_0," << lambda_0 << "\n";
    summary << "tan_delta," << tan_delta << "\n";
    summary << "k_0," << k_0 << "\n";
    summary << "d," << d << "\n";
    summary << "rho_min," << rho_min << "\n";
    summary << "rho_max," << rho_max << "\n";
    summary << "N_rho," << N_rho << "\n";
    summary << "norm," << norm << "\n";

    std::cout << "Arquivos gerados: v_pot_eps.csv e v_pot_eps_summary.csv\n";
    return 0;
}
