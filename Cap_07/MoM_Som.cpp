#include "include/chapter7_functions.hpp"
#include "include/chapter7_common.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace
{
struct Config
{
    int max_mode_index = 0;
    int num_int = 0;
    double freq_res = 10e9;
    double eps_r_prime = 2.55;
    double tan_delta = 0.0;
    double d_scale = 0.12;
    double L_scale = 0.39;
    double freq_start_scale = 0.7;
    double freq_stop_scale = 1.1;
    double freq_step_scale = 0.025;
};

std::string require_value(int &index, int argc, char **argv)
{
    if (index + 1 >= argc)
        throw std::runtime_error(std::string("Falta valor para ") + argv[index]);
    return argv[++index];
}

void print_help()
{
    std::cout
        << "Uso: ./build/MoM_Som --max-mode-index M --num-int N [opcoes]\n"
        << "  --max-mode-index M       numero maximo de modos (obrigatorio)\n"
        << "  --num-int N              numero de pontos de integracao, par (obrigatorio)\n"
        << "  --freq-res VALOR         frequencia de projeto em Hz (padrao 10e9)\n"
        << "  --eps-r-prime VALOR      permissividade relativa do substrato (padrao 2.55)\n"
        << "  --tan-delta VALOR        tangente de perdas (padrao 0.0)\n"
        << "  --d-scale VALOR          d/lambda_res (padrao 0.12)\n"
        << "  --L-scale VALOR          L/lambda_res (padrao 0.39)\n"
        << "  --freq-start-scale VALOR escala inicial do sweep em freq_res (padrao 0.7)\n"
        << "  --freq-stop-scale VALOR  escala final do sweep em freq_res (padrao 1.1)\n"
        << "  --freq-step-scale VALOR  passo do sweep em freq_res (padrao 0.025)\n"
        << "  --help                   mostra esta ajuda\n";
}

bool is_finite_complex(const std::complex<double> &value)
{
    return std::isfinite(value.real()) && std::isfinite(value.imag());
}

double condition_number(const Eigen::MatrixXcd &matrix)
{
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd(matrix);
    const auto singular_values = svd.singularValues();
    if (singular_values.size() == 0)
        return 0.0;

    const double sigma_max = singular_values.maxCoeff();
    const double sigma_min = singular_values.minCoeff();
    if (sigma_min <= 0.0)
        return std::numeric_limits<double>::infinity();
    return sigma_max / sigma_min;
}
} // namespace

int main(int argc, char **argv)
{
    // Translation of MoM_Som.m:
    // input impedance of a thin printed dipole over a grounded substrate.
    const fs::path out_dir = PROJECT_OUT_DIR;
    fs::create_directories(out_dir);

    Config config;
    try
    {
        for (int i = 1; i < argc; ++i)
        {
            const std::string arg = argv[i];
            if (arg == "--help" || arg == "-h")
            {
                print_help();
                return 0;
            }
            if (arg == "--max-mode-index")
                config.max_mode_index = std::stoi(require_value(i, argc, argv));
            else if (arg == "--num-int")
                config.num_int = std::stoi(require_value(i, argc, argv));
            else if (arg == "--freq-res")
                config.freq_res = std::stod(require_value(i, argc, argv));
            else if (arg == "--eps-r-prime")
                config.eps_r_prime = std::stod(require_value(i, argc, argv));
            else if (arg == "--tan-delta")
                config.tan_delta = std::stod(require_value(i, argc, argv));
            else if (arg == "--d-scale")
                config.d_scale = std::stod(require_value(i, argc, argv));
            else if (arg == "--L-scale")
                config.L_scale = std::stod(require_value(i, argc, argv));
            else if (arg == "--freq-start-scale")
                config.freq_start_scale = std::stod(require_value(i, argc, argv));
            else if (arg == "--freq-stop-scale")
                config.freq_stop_scale = std::stod(require_value(i, argc, argv));
            else if (arg == "--freq-step-scale")
                config.freq_step_scale = std::stod(require_value(i, argc, argv));
            else
                throw std::runtime_error("Opcao desconhecida: " + arg);
        }
    }
    catch (const std::exception &ex)
    {
        std::cerr << "Erro: " << ex.what() << "\n";
        print_help();
        return 1;
    }

    if (config.max_mode_index <= 0 || config.num_int <= 0 || (config.num_int % 2) != 0)
    {
        std::cerr << "Erro: --max-mode-index deve ser positivo e --num-int deve ser positivo e par.\n";
        return 1;
    }

    const double c = 2.997925e8;
    const double eps_0 = 8.854e-12;
    const double mu_0 = std::numbers::pi * 4e-7;
    const double lambda_res = c / config.freq_res;
    const double d = config.d_scale * lambda_res;
    const double L = config.L_scale * lambda_res;
    const std::complex<double> eps_r =
        config.eps_r_prime * (1.0 - std::complex<double>(0.0, config.tan_delta));

    const double del_x = L / static_cast<double>(config.num_int);

    std::vector<double> x_pot(config.num_int, 0.0);
    for (int index = 0; index < config.num_int; ++index)
        x_pot[index] = -L / 2.0 + del_x / 2.0 + index * del_x;

    std::vector<double> freq_vec;
    for (double factor = config.freq_start_scale;
         factor <= config.freq_stop_scale + 1e-12;
         factor += config.freq_step_scale)
    {
        freq_vec.push_back(factor * config.freq_res);
    }

    std::ofstream impedance(out_dir / "mom_som_impedance.csv");
    impedance << "freq_hz,freq_ghz,re_zin,im_zin,abs_zin,lambda_p,cond_z\n";
    impedance << std::scientific << std::setprecision(10);

    std::ofstream currents(out_dir / "mom_som_currents.csv");
    currents << "freq_hz,x_over_L,re_i,im_i,abs_i\n";
    currents << std::scientific << std::setprecision(10);

    std::ofstream vpot_debug(out_dir / "mom_som_vpot_table.csv");
    vpot_debug << "freq_hz,rho_m,re_vpot,im_vpot,abs_vpot\n";
    vpot_debug << std::scientific << std::setprecision(10);

    for (double freq : freq_vec)
    {
        const double omega = 2.0 * std::numbers::pi * freq;
        const double lambda_0 = c / freq;
        const double k_0 = 2.0 * std::numbers::pi * freq / c;

        SommerfeldContext context;
        context.k_0 = k_0;
        context.eps_r_prime = config.eps_r_prime;
        context.eps_r = eps_r;
        context.k = k_0 * std::sqrt(eps_r);
        context.d = d;
        context.lambda_p =
            root_D_TM(k_0, context.k, context.eps_r, d, 50);

        std::vector<double> rho_pot(config.num_int + 1, 0.0);
        std::vector<std::complex<double>> vpot_table(config.num_int + 1, 0.0);

        context.rho = std::abs(del_x / 4.0);
        rho_pot[0] = context.rho;
        vpot_table[0] = V_int(context);
        if (!is_finite_complex(vpot_table[0]))
            throw std::runtime_error("V_int retornou valor nao finito na primeira entrada da tabela.");
        vpot_debug << freq << "," << rho_pot[0] << ","
                   << vpot_table[0].real() << "," << vpot_table[0].imag() << ","
                   << std::abs(vpot_table[0]) << "\n";
        for (int jj = 1; jj < config.num_int; ++jj)
        {
            context.rho = std::abs(x_pot[jj] - x_pot[0]);
            rho_pot[jj] = context.rho;
            vpot_table[jj] = V_int(context);
            if (!is_finite_complex(vpot_table[jj]))
            {
                throw std::runtime_error(
                    "V_int retornou valor nao finito na tabela de interpolacao para freq = " +
                    std::to_string(freq) + " Hz.");
            }
            vpot_debug << freq << "," << rho_pot[jj] << ","
                       << vpot_table[jj].real() << "," << vpot_table[jj].imag() << ","
                       << std::abs(vpot_table[jj]) << "\n";
        }
        context.rho = L;
        rho_pot[config.num_int] = context.rho;
        vpot_table[config.num_int] = V_int(context);
        if (!is_finite_complex(vpot_table[config.num_int]))
            throw std::runtime_error("V_int retornou valor nao finito no ponto extra da tabela.");
        vpot_debug << freq << "," << rho_pot[config.num_int] << ","
                   << vpot_table[config.num_int].real() << ","
                   << vpot_table[config.num_int].imag() << ","
                   << std::abs(vpot_table[config.num_int]) << "\n";

        std::vector<double> x(config.num_int, 0.0);
        std::vector<double> xp(config.num_int, 0.0);
        for (int jj = 0; jj < config.num_int; ++jj)
        {
            x[jj] = -L / 2.0 + del_x / 3.0 + jj * del_x;
            xp[jj] = -L / 2.0 + 2.0 * del_x / 3.0 + jj * del_x;
        }

        std::vector<std::vector<std::complex<double>>> Apot(
            config.num_int, std::vector<std::complex<double>>(config.num_int, 0.0));
        std::vector<std::vector<std::complex<double>>> Vpot(
            config.num_int, std::vector<std::complex<double>>(config.num_int, 0.0));

        for (int jj = 0; jj < config.num_int; ++jj)
        {
            for (int kk = 0; kk < config.num_int; ++kk)
            {
                const double rho = std::abs(x[jj] - xp[kk]);
                const double R_0 = rho;
                const double R_1 = std::sqrt(rho * rho + (2.0 * d) * (2.0 * d));
                Apot[jj][kk] =
                    std::exp(-std::complex<double>(0.0, 1.0) * k_0 * R_0) / R_0 -
                    std::exp(-std::complex<double>(0.0, 1.0) * k_0 * R_1) / R_1;
                Vpot[jj][kk] = interp1_linear(rho_pot, vpot_table, rho);
                if (!is_finite_complex(Vpot[jj][kk]))
                    throw std::runtime_error("Interpolacao do potencial escalar produziu NaN/Inf.");
            }
        }

        Eigen::MatrixXcd A_mat =
            Eigen::MatrixXcd::Zero(config.max_mode_index, config.max_mode_index);
        Eigen::MatrixXcd V_mat =
            Eigen::MatrixXcd::Zero(config.max_mode_index, config.max_mode_index);

        std::vector<std::complex<double>> amn_inner(config.num_int, 0.0);
        std::vector<std::complex<double>> vmn_inner(config.num_int, 0.0);
        std::vector<std::complex<double>> amn_outer(config.num_int, 0.0);
        std::vector<std::complex<double>> vmn_outer(config.num_int, 0.0);

        for (int m_index = 0; m_index < config.max_mode_index; ++m_index)
        {
            const int mm = 2 * (m_index + 1) - 1;
            for (int n_index = 0; n_index < config.max_mode_index; ++n_index)
            {
                const int nn = 2 * (n_index + 1) - 1;
                for (int jj = 0; jj < config.num_int; ++jj)
                {
                    for (int kk = 0; kk < config.num_int; ++kk)
                    {
                        amn_inner[kk] =
                            std::cos(nn * std::numbers::pi * xp[kk] / L) * Apot[jj][kk];
                        vmn_inner[kk] =
                            std::sin(nn * std::numbers::pi * xp[kk] / L) * Vpot[jj][kk];
                    }
                    amn_outer[jj] =
                        std::cos(mm * std::numbers::pi * x[jj] / L) *
                        del_x * trapz_uniform(amn_inner);
                    vmn_outer[jj] =
                        std::sin(mm * std::numbers::pi * x[jj] / L) *
                        del_x * trapz_uniform(vmn_inner);
                }

                A_mat(m_index, n_index) = del_x * trapz_uniform(amn_outer);
                V_mat(m_index, n_index) =
                    (static_cast<double>(mm * nn) * std::numbers::pi * std::numbers::pi) /
                    (L * L) * del_x * trapz_uniform(vmn_outer);
            }
        }

        A_mat *= std::complex<double>(0.0, 1.0) * omega * (mu_0 / (4.0 * std::numbers::pi));
        V_mat *= 1.0 / (std::complex<double>(0.0, 1.0) * omega) *
                 (1.0 / (2.0 * std::numbers::pi * eps_0));

        Eigen::VectorXcd b_vec =
            Eigen::VectorXcd::Ones(config.max_mode_index);

        const Eigen::MatrixXcd Z_mat = A_mat + V_mat;
        const double cond_z = condition_number(Z_mat);
        for (int row = 0; row < Z_mat.rows(); ++row)
        {
            for (int col = 0; col < Z_mat.cols(); ++col)
            {
                if (!is_finite_complex(Z_mat(row, col)))
                    throw std::runtime_error("A matriz Z contem NaN/Inf.");
            }
        }

        const Eigen::VectorXcd I_vec = Z_mat.fullPivLu().solve(b_vec);
        for (int index = 0; index < I_vec.size(); ++index)
        {
            if (!is_finite_complex(I_vec(index)))
            {
                throw std::runtime_error(
                    "A solucao MoM produziu NaN/Inf para freq = " +
                    std::to_string(freq) + " Hz.");
            }
        }

        std::vector<double> x_s(config.num_int, 0.0);
        std::vector<std::complex<double>> I_s(config.num_int, 0.0);
        for (int jj = 0; jj < config.num_int; ++jj)
        {
            x_s[jj] = -L / 2.0 + del_x / 2.0 + jj * del_x;
            for (int m_index = 0; m_index < config.max_mode_index; ++m_index)
            {
                const int mm = 2 * (m_index + 1) - 1;
                I_s[jj] += I_vec(m_index) *
                           std::cos(mm * std::numbers::pi * x_s[jj] / L);
            }
        }

        // The MATLAB script samples I_s(num_int/2), which in 1-based indexing
        // means the point immediately to the left of the strip centre.
        const int centre_index = config.num_int / 2 - 1;
        const std::complex<double> Zin = 1.0 / I_s[centre_index];

        impedance << freq << "," << freq / 1e9 << ","
                  << Zin.real() << "," << Zin.imag() << ","
                  << std::abs(Zin) << "," << context.lambda_p << ","
                  << cond_z << "\n";

        for (int jj = 0; jj < config.num_int; ++jj)
        {
            currents << freq << "," << x_s[jj] / L << ","
                     << I_s[jj].real() << "," << I_s[jj].imag() << ","
                     << std::abs(I_s[jj]) << "\n";
        }

        (void)lambda_0;
    }

    std::ofstream summary(out_dir / "mom_som_summary.csv");
    summary << "key,value\n";
    summary << "freq_res," << config.freq_res << "\n";
    summary << "eps_r_prime," << config.eps_r_prime << "\n";
    summary << "tan_delta," << config.tan_delta << "\n";
    summary << "d," << d << "\n";
    summary << "L," << L << "\n";
    summary << "max_mode_index," << config.max_mode_index << "\n";
    summary << "num_int," << config.num_int << "\n";
    summary << "freq_start_scale," << config.freq_start_scale << "\n";
    summary << "freq_stop_scale," << config.freq_stop_scale << "\n";
    summary << "freq_step_scale," << config.freq_step_scale << "\n";

    std::cout << "Arquivos gerados: mom_som_impedance.csv, mom_som_currents.csv e mom_som_summary.csv\n";
    return 0;
}
