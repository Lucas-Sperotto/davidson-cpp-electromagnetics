#include <Eigen/Dense>

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <string>

namespace fs = std::filesystem;

namespace
{
// Translation of static_mom.m:
// compute the charge distribution along a straight wire.
struct Config
{
    double L = 1.0;
    double a = 0.001;
    int N = 5;
    double voltage = 1.0;
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
        << "Uso: ./build/static_mom [opcoes]\n"
        << "  --length VALOR      comprimento do fio em metros (padrao 1.0)\n"
        << "  --radius VALOR      raio do fio em metros (padrao 0.001)\n"
        << "  --segments N        numero de segmentos (padrao 5)\n"
        << "  --voltage VALOR     potencial do fio em volts (padrao 1.0)\n"
        << "  --help              mostra esta ajuda\n";
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
            if (arg == "--length")
                config.L = std::stod(require_value(i, argc, argv));
            else if (arg == "--radius")
                config.a = std::stod(require_value(i, argc, argv));
            else if (arg == "--segments")
                config.N = std::stoi(require_value(i, argc, argv));
            else if (arg == "--voltage")
                config.voltage = std::stod(require_value(i, argc, argv));
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

    if (config.L <= 0.0 || config.a <= 0.0 || config.N <= 0)
    {
        std::cerr << "Erro: comprimento, raio e numero de segmentos devem ser positivos.\n";
        return 1;
    }

    // Define geometry.
    const double delta = config.L / config.N;
    const double eps_0 = 8.854e-12;
    const double voltage = config.voltage; // Potential of wire [V]

    Eigen::VectorXd V_vec = Eigen::VectorXd::Zero(config.N);
    Eigen::VectorXd I_vec = Eigen::VectorXd::Zero(config.N);
    Eigen::MatrixXd Z_mat = Eigen::MatrixXd::Zero(config.N, config.N);

    // Set up matrix equation.
    for (int mm = 0; mm < config.N; ++mm)
    {
        for (int nn = 0; nn < config.N; ++nn)
        {
            const double l_m = std::sqrt(std::pow((mm - nn) * delta, 2) + config.a * config.a);
            const double d_plus_mn = l_m + delta / 2.0;
            const double d_min_mn = l_m - delta / 2.0;
            if (mm == nn)
            {
                Z_mat(mm, nn) = 2.0 * std::log((delta / 2.0 + std::sqrt(config.a * config.a + (delta / 2.0) * (delta / 2.0))) / config.a);
            }
            else if (std::abs(mm - nn) <= 2)
            {
                Z_mat(mm, nn) = std::log((d_plus_mn + std::sqrt(d_plus_mn * d_plus_mn + config.a * config.a)) /
                                         (d_min_mn + std::sqrt(d_min_mn * d_min_mn + config.a * config.a)));
            }
            else
            {
                Z_mat(mm, nn) = std::log(d_plus_mn / d_min_mn);
            }
        }
    }

    V_vec = 4.0 * std::numbers::pi * eps_0 * voltage * Eigen::VectorXd::Ones(config.N);
    I_vec = Z_mat.partialPivLu().solve(V_vec);

    std::ofstream charge(out_dir / "static_mom_charge.csv");
    charge << "z_m,line_charge_C_per_m,line_charge_pC_per_m\n";
    charge << std::scientific << std::setprecision(10);
    for (int i = 0; i < config.N; ++i)
    {
        const double z_axis = delta / 2.0 + i * delta;
        charge << z_axis << "," << I_vec(i) << "," << I_vec(i) / 1e-12 << "\n";
    }

    std::ofstream summary(out_dir / "static_mom_summary.csv");
    summary << "key,value\n";
    summary << "L," << config.L << "\n";
    summary << "a," << config.a << "\n";
    summary << "N," << config.N << "\n";
    summary << "Delta," << delta << "\n";
    summary << "eps_0," << eps_0 << "\n";
    summary << "voltage," << voltage << "\n";

    std::cout << "Arquivos gerados: static_mom_charge.csv e static_mom_summary.csv\n";
    return 0;
}
