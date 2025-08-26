#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include <numbers> // aqui está o std::numbers::pi
#include "FEM_1D_solver.h"
#include "FEM_pp.h"

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

int main(int argc, char **argv)
{
    // Parâmetros padrão
    int N_elem_init = 2;
    int num_meshes = 1;
    double f = 1.0; // Hz
    double L = 1.0;
    double C = 1.0;
    double V_in = 1.0; // Voltage at source end of tx line.

    // Leitura de argumentos opcionais:
    // ./fem_tl [N_elem_init] [num_meshes] [f] [L] [C] [V_in]
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

    // Constantes e grandezas derivadas
    // const double pi = 3.14159265358979323846;
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
              << ", num_meshes=" << num_meshes << "\n";

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
            save_csv(name.str(),
                     {z, V_fem, V_anl},
                     {"z", "V_fem", "V_analitico"});
        }

        // (Opcional) salvar valores nos nós
        {
            std::ostringstream name;
            name << "fem_nodes_stage_" << stage << ".csv";
            std::vector<double> zn(V.size());
            for (size_t k = 0; k < V.size(); ++k)
                zn[k] = h * k;
            save_csv(name.str(),
                     {zn, V},
                     {"z_node", "V_node"});
        }

        std::cout << "Stage " << stage << ": N_elem=" << N_elem
                  << ", h/lambda=" << (h / lambda)
                  << ", RMS err=" << rms_err << "\n";

        // Refinamento
        N_elem *= 2;
    }

    // Salva curva de convergência
    save_csv("fem_convergence.csv",
             {conv_h_over_lambda, conv_rms_err},
             {"h_over_lambda", "rms_err"});

    std::cout << "Arquivos gerados:\n"
              << "  fem_profile_stage_#.csv (z, V_fem, V_analitico)\n"
              << "  fem_nodes_stage_#.csv   (z_node, V_node)\n"
              << "  fem_convergence.csv     (h_over_lambda, rms_err)\n";
    return 0;
}