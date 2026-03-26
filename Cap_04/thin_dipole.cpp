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
// Translation of thin_dipole.m:
// thin-wire MoM with piecewise sinusoidal basis functions, including
// magnetic frill and delta-gap feed models.
struct Config
{
    int n_seg = 80;
    double length_lambda = 0.48;
    double radius_lambda = 0.005;
    double b = 0.0;
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
        << "Uso: ./build/thin_dipole [opcoes]\n"
        << "  --n-seg N            numero de segmentos (deve ser par, padrao 80)\n"
        << "  --length L           comprimento do fio em comprimentos de onda (padrao 0.48)\n"
        << "  --radius A           raio do fio em comprimentos de onda (padrao 0.005)\n"
        << "  --b VALOR            raio externo da franja magnetica (padrao 3.5*a)\n"
        << "  --help               mostra esta ajuda\n";
}

Eigen::MatrixXcd toeplitz_from_first_row(const Eigen::VectorXcd &row)
{
    const int n = static_cast<int>(row.size());
    Eigen::MatrixXcd matrix(n, n);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
            matrix(i, j) = row(std::abs(j - i));
    }
    return matrix;
}
} // namespace

int main(int argc, char **argv)
{
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
            if (arg == "--n-seg")
                config.n_seg = std::stoi(require_value(i, argc, argv));
            else if (arg == "--length")
                config.length_lambda = std::stod(require_value(i, argc, argv));
            else if (arg == "--radius")
                config.radius_lambda = std::stod(require_value(i, argc, argv));
            else if (arg == "--b")
                config.b = std::stod(require_value(i, argc, argv));
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

    if (config.n_seg <= 0 || (config.n_seg % 2) != 0 || config.length_lambda <= 0.0 || config.radius_lambda <= 0.0)
    {
        std::cerr << "Erro: n_seg deve ser par e positivo; comprimento e raio devem ser positivos.\n";
        return 1;
    }

    if (config.b <= 0.0)
        config.b = 3.5 * config.radius_lambda;

    const std::complex<double> imag_unit(0.0, 1.0);
    const double beta = 2.0 * std::numbers::pi; // lambda = 1 m
    const double dz = config.length_lambda / static_cast<double>(config.n_seg);
    const double beta_dz = beta * dz;
    const int N = config.n_seg - 1;

    Eigen::VectorXcd c(N);
    for (int n = 1; n <= N; ++n)
    {
        const double r_np1 = std::sqrt(config.radius_lambda * config.radius_lambda + std::pow((n - 2) * dz, 2));
        const double r_n = std::sqrt(config.radius_lambda * config.radius_lambda + std::pow((n - 1) * dz, 2));
        const double r_nm1 = std::sqrt(config.radius_lambda * config.radius_lambda + std::pow(n * dz, 2));
        const std::complex<double> A = std::exp(-imag_unit * beta * r_np1) / r_np1;
        const std::complex<double> B =
            std::exp(-imag_unit * beta * r_n) * std::sin(2.0 * beta_dz) / (r_n * std::sin(beta_dz));
        const std::complex<double> C = std::exp(-imag_unit * beta * r_nm1) / r_nm1;
        c(n - 1) = -imag_unit * 30.0 * (A - B + C) / std::sin(beta_dz);
    }

    const Eigen::MatrixXcd Z = toeplitz_from_first_row(c);
    const int centre_seg = (N - 1) / 2;

    Eigen::VectorXcd E = Eigen::VectorXcd::Zero(N);
    for (int m = centre_seg; m < N; ++m)
    {
        const double zz = (m - centre_seg) * dz;
        const double r1 = std::sqrt(config.radius_lambda * config.radius_lambda + zz * zz);
        const double r2 = std::sqrt(config.b * config.b + zz * zz);
        E(m) = (std::exp(-imag_unit * beta * r1) / r1 - std::exp(-imag_unit * beta * r2) / r2) /
               (2.0 * std::log(config.b / config.radius_lambda));
        E(N - 1 - m) = E(m);
    }

    Eigen::VectorXcd V = -E.conjugate();
    const Eigen::VectorXcd I_mag_frill = Z.partialPivLu().solve(V);
    const std::complex<double> Zin_mag_frill = 1.0 / I_mag_frill(centre_seg);

    // Delta-gap feed model.
    Eigen::VectorXcd V_delta = Eigen::VectorXcd::Zero(N);
    V_delta(centre_seg) = -1.0 / dz;
    const Eigen::VectorXcd I_delta_gap = Z.partialPivLu().solve(V_delta);
    const std::complex<double> Zin_delta_gap = 1.0 / I_delta_gap(centre_seg);

    std::ofstream current(out_dir / "thin_dipole_current.csv");
    current << "z_over_lambda,abs_i_delta_gap,abs_i_mag_frill,re_i_delta_gap,im_i_delta_gap,re_i_mag_frill,im_i_mag_frill\n";
    current << std::scientific << std::setprecision(10);
    for (int i = 0; i < N; ++i)
    {
        const double z = -config.length_lambda / 2.0 + dz / 2.0 + i * dz;
        current << z << ","
                << std::abs(I_delta_gap(i)) << ","
                << std::abs(I_mag_frill(i)) << ","
                << I_delta_gap(i).real() << ","
                << I_delta_gap(i).imag() << ","
                << I_mag_frill(i).real() << ","
                << I_mag_frill(i).imag() << "\n";
    }

    const double gamma_frill =
        20.0 * std::log10(std::abs((Zin_mag_frill - std::complex<double>(75.0, 0.0)) /
                                   (Zin_mag_frill + std::complex<double>(75.0, 0.0))));
    const double gamma_delta =
        20.0 * std::log10(std::abs((Zin_delta_gap - std::complex<double>(75.0, 0.0)) /
                                   (Zin_delta_gap + std::complex<double>(75.0, 0.0))));

    std::ofstream summary(out_dir / "thin_dipole_summary.csv");
    summary << "key,value\n";
    summary << "n_seg," << config.n_seg << "\n";
    summary << "length_lambda," << config.length_lambda << "\n";
    summary << "radius_lambda," << config.radius_lambda << "\n";
    summary << "b," << config.b << "\n";
    summary << "dz," << dz << "\n";
    summary << "zin_mag_frill_real," << Zin_mag_frill.real() << "\n";
    summary << "zin_mag_frill_imag," << Zin_mag_frill.imag() << "\n";
    summary << "zin_delta_gap_real," << Zin_delta_gap.real() << "\n";
    summary << "zin_delta_gap_imag," << Zin_delta_gap.imag() << "\n";
    summary << "gamma_frill_db," << gamma_frill << "\n";
    summary << "gamma_delta_db," << gamma_delta << "\n";

    std::cout << "Arquivos gerados: thin_dipole_current.csv e thin_dipole_summary.csv\n";
    return 0;
}
