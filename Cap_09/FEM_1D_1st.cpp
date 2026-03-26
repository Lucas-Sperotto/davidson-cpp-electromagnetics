#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include <numbers>
#include <filesystem>
#include <sstream>
#include <string>
#include <stdexcept>
#include "FEM_1D_solver.h"
#include "FEM_pp.h"

namespace fs = std::filesystem;

// Salva CSV genérico
static void save_csv(const std::string &filename,
                     const std::vector<std::vector<double>> &cols,
                     const std::vector<std::string> &headers)
{
    std::ofstream f(filename);
    f << std::setprecision(10) << std::fixed;
    if (!headers.empty())
    {
        for (size_t i = 0; i < headers.size(); ++i)
        {
            f << headers[i];
            if (i + 1 < headers.size())
                f << ",";
        }
        f << "\n";
    }
    if (cols.empty())
        return;
    size_t rows = cols[0].size();
    for (size_t r = 0; r < rows; ++r)
    {
        for (size_t c = 0; c < cols.size(); ++c)
        {
            f << cols[c][r];
            if (c + 1 < cols.size())
                f << ",";
        }
        f << "\n";
    }
}

namespace
{
std::string require_value(int &index, int argc, char **argv)
{
    if (index + 1 >= argc)
        throw std::runtime_error(std::string("Falta valor para ") + argv[index]);
    return argv[++index];
}

void print_help()
{
    std::cout
        << "Uso: ./build/fem_tl [opcoes]\n"
        << "  --n-elem-init N    numero inicial de elementos (padrao 2)\n"
        << "  --num-meshes N     numero de refinamentos (padrao 1)\n"
        << "  --freq F           frequencia em Hz (padrao 1.0)\n"
        << "  --inductance L     indutancia por unidade de comprimento (padrao 1.0)\n"
        << "  --capacitance C    capacitancia por unidade de comprimento (padrao 1.0)\n"
        << "  --vin V            tensao prescrita no ultimo no (padrao 1.0)\n"
        << "  --help             mostra esta ajuda\n"
        << "\n"
        << "Compatibilidade: tambem aceita a forma legada posicional\n"
        << "  ./build/fem_tl [N_elem_init] [num_meshes] [f] [L] [C] [V_in]\n";
}
} // namespace

int main(int argc, char **argv)
{
    const fs::path out_dir = PROJECT_OUT_DIR;
    fs::create_directories(out_dir);

    // Parametros padrao
    int N_elem_init = 2;
    int num_meshes = 1;
    double f = 1.0; // Hz
    double L = 1.0;
    double C = 1.0;
    double V_in = 1.0;

    try
    {
        const bool legacy_positional =
            argc > 1 && std::string(argv[1]).rfind("--", 0) != 0;

        if (legacy_positional)
        {
            if (argc >= 2)
                N_elem_init = std::max(1, std::atoi(argv[1]));
            if (argc >= 3)
                num_meshes = std::max(1, std::atoi(argv[2]));
            if (argc >= 4)
                f = std::atof(argv[3]);
            if (argc >= 5)
                L = std::atof(argv[4]);
            if (argc >= 6)
                C = std::atof(argv[5]);
            if (argc >= 7)
                V_in = std::atof(argv[6]);
        }
        else
        {
            for (int i = 1; i < argc; ++i)
            {
                const std::string arg = argv[i];
                if (arg == "--help" || arg == "-h")
                {
                    print_help();
                    return 0;
                }
                if (arg == "--n-elem-init")
                    N_elem_init = std::max(1, std::stoi(require_value(i, argc, argv)));
                else if (arg == "--num-meshes")
                    num_meshes = std::max(1, std::stoi(require_value(i, argc, argv)));
                else if (arg == "--freq")
                    f = std::stod(require_value(i, argc, argv));
                else if (arg == "--inductance")
                    L = std::stod(require_value(i, argc, argv));
                else if (arg == "--capacitance")
                    C = std::stod(require_value(i, argc, argv));
                else if (arg == "--vin")
                    V_in = std::stod(require_value(i, argc, argv));
                else
                    throw std::runtime_error("Opcao desconhecida: " + arg);
            }
        }
    }
    catch (const std::exception &ex)
    {
        std::cerr << "Erro: " << ex.what() << "\n";
        print_help();
        return 1;
    }

    if (f <= 0.0 || L <= 0.0 || C <= 0.0)
    {
        std::cerr << "Erro: freq, L e C devem ser positivos.\n";
        return 1;
    }

    // Constantes e grandezas derivadas
    double c = 1.0 / std::sqrt(L * C);
    double omega = 2.0 * std::numbers::pi * f;
    double lambda = c / f;
    double ell = lambda / 2.0; // comprimento da linha = lambda/2
    double beta = 2.0 * std::numbers::pi / lambda;

    std::cout << "Parametros:\n"
              << "  f=" << f << " Hz, L=" << L << ", C=" << C
              << ", c=" << c << " m/s, lambda=" << lambda
              << ", ell=lambda/2=" << ell << ", V_in=" << V_in << "\n"
              << "  N_elem_init=" << N_elem_init
              << ", num_meshes=" << num_meshes << "\n"
              << "  out_dir=" << out_dir << "\n";

    std::vector<double> conv_h_over_lambda;
    std::vector<double> conv_rms_err;

    int N_elem = N_elem_init;
    for (int stage = 1; stage <= num_meshes; ++stage)
    {
        double h = ell / static_cast<double>(N_elem);

        // Solver FEM no grid atual
        std::vector<double> V = FEM_1D_solver(N_elem, h, omega, L, C, V_in);

        // Pós-processamento em N_int pontos
        int N_int = N_elem * 10;
        std::vector<double> z;
        std::vector<double> V_fem = FEM_postprocess(V, N_elem, ell, h, N_int, z);

        // "Analitico"
        std::vector<double> V_anl(N_int, 0.0);
        for (int i = 0; i < N_int; ++i)
            V_anl[i] = -std::cos(2.0 * beta * z[i] * ell);

        // Erro RMS
        double se = 0.0;
        for (int i = 0; i < N_int; ++i)
        {
            double e = V_anl[i] - V_fem[i];
            se += e * e;
        }
        double rms_err = std::sqrt(se / static_cast<double>(N_int));

        conv_h_over_lambda.push_back(h / lambda);
        conv_rms_err.push_back(rms_err);

        // Salva perfil FEM vs Analítico para este estágio
        {
            std::ostringstream name;
            name << "fem_profile_stage_" << stage << ".csv";
            save_csv((out_dir / name.str()).string(),
                     {z, V_fem, V_anl},
                     {"z", "V_fem", "V_analitico"});
        }

        // Salva valores nos nos
        {
            std::ostringstream name;
            name << "fem_nodes_stage_" << stage << ".csv";
            std::vector<double> zn(V.size());
            for (size_t k = 0; k < V.size(); ++k)
                zn[k] = h * k;
            save_csv((out_dir / name.str()).string(),
                     {zn, V},
                     {"z_node", "V_node"});
        }

        std::cout << "Stage " << stage << ": N_elem=" << N_elem
                  << ", h/lambda=" << (h / lambda)
                  << ", RMS err=" << rms_err << "\n";

        // Refinamento
        N_elem *= 2;
    }

    save_csv((out_dir / "fem_convergence.csv").string(),
             {conv_h_over_lambda, conv_rms_err},
             {"h_over_lambda", "rms_err"});

    {
        std::ofstream meta(out_dir / "fem_run_metadata.csv");
        meta << std::setprecision(10) << std::fixed;
        meta << "key,value\n";
        meta << "N_elem_init," << N_elem_init << "\n";
        meta << "num_meshes," << num_meshes << "\n";
        meta << "freq_hz," << f << "\n";
        meta << "L," << L << "\n";
        meta << "C," << C << "\n";
        meta << "c," << c << "\n";
        meta << "lambda," << lambda << "\n";
        meta << "ell," << ell << "\n";
        meta << "V_in," << V_in << "\n";
        if (conv_h_over_lambda.size() > 1)
        {
            std::vector<double> log_h(conv_h_over_lambda.size());
            std::vector<double> log_err(conv_rms_err.size());
            for (size_t i = 0; i < conv_h_over_lambda.size(); ++i)
            {
                log_h[i] = std::log(conv_h_over_lambda[i]);
                log_err[i] = std::log(conv_rms_err[i]);
            }

            double sum_x = 0.0;
            double sum_y = 0.0;
            double sum_xx = 0.0;
            double sum_xy = 0.0;
            for (size_t i = 0; i < log_h.size(); ++i)
            {
                sum_x += log_h[i];
                sum_y += log_err[i];
                sum_xx += log_h[i] * log_h[i];
                sum_xy += log_h[i] * log_err[i];
            }
            const double n = static_cast<double>(log_h.size());
            const double denom = n * sum_xx - sum_x * sum_x;
            if (std::abs(denom) > 1e-14)
            {
                const double slope = (n * sum_xy - sum_x * sum_y) / denom;
                meta << "loglog_slope," << slope << "\n";
            }
        }
    }

    std::cout << "Arquivos gerados:\n"
              << "  " << (out_dir / "fem_profile_stage_#.csv") << " (z, V_fem, V_analitico)\n"
              << "  " << (out_dir / "fem_nodes_stage_#.csv") << "   (z_node, V_node)\n"
              << "  " << (out_dir / "fem_convergence.csv") << "     (h_over_lambda, rms_err)\n"
              << "  " << (out_dir / "fem_run_metadata.csv") << "    (metadados da rodada)\n";
    return 0;
}
