#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <tuple>
#include <filesystem>
#include <iomanip>
#include <string>
#include <stdexcept>

#include "trimesh.h"
#include "edgemake.h"
#include "free_dof.h"
#include "renumber_dof_LTQN.h"
#include "sandt_LTQN.h"
#include "free_nodes.h"
#include "TEeig_err.h"

using namespace Eigen;

#include "globals.h"

const double eps_0 = 8.854e-12;
const double mu_0 = 4 * M_PI * 1e-7;
const double eps_r = 1.0;
const double mu_r = 1.0;

namespace
{
struct Eigen2DLtqnConfig
{
    double a = 2.286e-2;
    double b = 1.016e-2;
    int x_mesh = 4;
    int y_mesh = 2;
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
        << "Uso: ./build/Eigen2D_LTQN [opcoes]\n"
        << "  --a VALOR        largura do guia (padrao 2.286e-2)\n"
        << "  --b VALOR        altura do guia (padrao 1.016e-2)\n"
        << "  --x-mesh N       divisao da malha em x (padrao 4)\n"
        << "  --y-mesh N       divisao da malha em y (padrao 2)\n"
        << "  --help           mostra esta ajuda\n";
}

void write_csv(const std::string &path, const Eigen::MatrixXd &M)
{
    std::ofstream f(path);
    f.setf(std::ios::scientific);
    f << std::setprecision(17);
    for (int i = 0; i < M.rows(); ++i)
    {
        for (int j = 0; j < M.cols(); ++j)
        {
            if (j)
                f << ',';
            f << M(i, j);
        }
        f << '\n';
    }
}

void write_csv(const std::string &path, const Eigen::VectorXd &V)
{
    std::ofstream f(path);
    f.setf(std::ios::scientific);
    f << std::setprecision(17);
    for (int i = 0; i < V.size(); ++i)
        f << V(i) << '\n';
}
} // namespace

void Eigen2D_LTQN(const Eigen2DLtqnConfig &config)
{
    const std::filesystem::path out_dir = PROJECT_OUT_DIR;
    std::filesystem::create_directories(out_dir);

    // Dimensões do guia
    const double a = config.a;
    const double b = config.b;

    std::cout << "Parametros do Eigen2D_LTQN:\n"
              << "  a=" << a << ", b=" << b
              << ", x_mesh=" << config.x_mesh << ", y_mesh=" << config.y_mesh << "\n";

    // Gerar a malha triangular
    std::vector<std::vector<double>> x_nodes, y_nodes;
    trimesh(a, b, config.x_mesh, config.y_mesh, x_nodes, y_nodes);
    // A função também deve preencher ELEMENTS e NODE_COORD

    // Criar as arestas e associar com os elementos
    edgemake(); // Função para preencher EDGES, ELEMENT_EDGES e NUM_EDGES

    // Determina os graus de liberdade livres nas arestas (PEC boundaries)
    std::vector<int> edge_free_flag = free_dof(a, b);

    // Renumera os graus de liberdade: CTLN, LTLN, faces LTQN
    std::vector<int> dof_e1, dof_e2;
    std::vector<int> dof_f1(NUM_ELEMENTS), dof_f2(NUM_ELEMENTS);
    renumber_dof_LTQN(edge_free_flag, dof_e1, dof_e2, dof_f1, dof_f2); // Atualiza NUM_DOFS também

    // Inicializa matrizes globais S e T
    std::vector<std::vector<double>> S(NUM_DOFS, std::vector<double>(NUM_DOFS, 0.0));
    std::vector<std::vector<double>> T(NUM_DOFS, std::vector<double>(NUM_DOFS, 0.0));

    // Loop sobre todos os elementos para montar as matrizes globais
    for (int ielem = 0; ielem < NUM_ELEMENTS; ++ielem)
    {
        std::vector<int> trinodes = ELEMENTS[ielem];

        double x1 = NODE_COORD[trinodes[0]][0], y1 = NODE_COORD[trinodes[0]][1];
        double x2 = NODE_COORD[trinodes[1]][0], y2 = NODE_COORD[trinodes[1]][1];
        double x3 = NODE_COORD[trinodes[2]][0], y3 = NODE_COORD[trinodes[2]][1];

        auto [S_elem, T_elem] = sandt_LTQN(x1, y1, x2, y2, x3, y3); // retorna matrizes 8x8

        int ll_f1 = dof_f1[ielem];
        int ll_f2 = dof_f2[ielem];

        for (int jedge = 0; jedge < 3; ++jedge)
        {
            int jj_e1 = dof_e1[ELEMENT_EDGES[ielem][jedge]];
            int jj_e2 = dof_e2[ELEMENT_EDGES[ielem][jedge]];

            for (int kedge = 0; kedge < 3; ++kedge)
            {
                int kk_e1 = dof_e1[ELEMENT_EDGES[ielem][kedge]];
                int kk_e2 = dof_e2[ELEMENT_EDGES[ielem][kedge]];

                if (jj_e1 >= 0 && jj_e2 >= 0 && kk_e1 >= 0 && kk_e2 >= 0)
                {
                    S[jj_e1][kk_e1] += S_elem[jedge][kedge];
                    S[jj_e1][kk_e2] += S_elem[jedge][kedge + 3];
                    S[jj_e2][kk_e2] += S_elem[jedge + 3][kedge + 3];

                    T[jj_e1][kk_e1] += T_elem[jedge][kedge];
                    T[jj_e1][kk_e2] += T_elem[jedge][kedge + 3];
                    T[jj_e2][kk_e2] += T_elem[jedge + 3][kedge + 3];
                }
            }

            if (jj_e1 >= 0 && jj_e2 >= 0)
            {
                S[jj_e1][ll_f1] += S_elem[jedge][6];
                S[jj_e1][ll_f2] += S_elem[jedge][7];
                S[jj_e2][ll_f1] += S_elem[jedge + 3][6];
                S[jj_e2][ll_f2] += S_elem[jedge + 3][7];

                T[jj_e1][ll_f1] += T_elem[jedge][6];
                T[jj_e1][ll_f2] += T_elem[jedge][7];
                T[jj_e2][ll_f1] += T_elem[jedge + 3][6];
                T[jj_e2][ll_f2] += T_elem[jedge + 3][7];
            }
        }

        S[ll_f1][ll_f1] += S_elem[6][6];
        S[ll_f1][ll_f2] += S_elem[6][7];
        S[ll_f2][ll_f2] += S_elem[7][7];

        T[ll_f1][ll_f1] += T_elem[6][6];
        T[ll_f1][ll_f2] += T_elem[6][7];
        T[ll_f2][ll_f2] += T_elem[7][7];
    }

    // Preencher as partes inferiores simétricas
    for (int j = 1; j < NUM_DOFS; ++j)
    {
        for (int k = 0; k < j; ++k)
        {
            S[j][k] = S[k][j];
            T[j][k] = T[k][j];
        }
    }

    // Converter std::vector para Eigen::MatrixXd
    MatrixXd S_mat(NUM_DOFS, NUM_DOFS), T_mat(NUM_DOFS, NUM_DOFS);
    for (int i = 0; i < NUM_DOFS; ++i)
    {
        for (int j = 0; j < NUM_DOFS; ++j)
        {
            S_mat(i, j) = S[i][j] / mu_r;
            T_mat(i, j) = T[i][j] * eps_r;
        }
    }

    write_csv((out_dir / "cpp_S_ltqn_mat.csv").string(), S_mat);
    write_csv((out_dir / "cpp_T_ltqn_mat.csv").string(), T_mat);

    // Resolver o problema generalizado de autovalores
    GeneralizedSelfAdjointEigenSolver<MatrixXd> solver(S_mat, T_mat);
    if (solver.info() != Eigen::Success)
    {
        std::cerr << "Falha ao resolver o problema generalizado LTQN.\n";
        return;
    }

    VectorXd eigvals = solver.eigenvalues();
    MatrixXd eigvecs = solver.eigenvectors();

    write_csv((out_dir / "cpp_eigvals_ltqn.csv").string(), eigvals);
    write_csv((out_dir / "cpp_eigvecs_ltqn.csv").string(), eigvecs);

    constexpr double zero_tol = 1e-6;

    // Extração dos autovalores e ordenação
    std::vector<std::pair<double, int>> indexed_kc;
    for (int i = 0; i < eigvals.size(); ++i)
    {
        double lambda = eigvals[i];
        if (std::abs(lambda) <= zero_tol)
            lambda = 0.0;
        if (lambda >= 0.0)
            indexed_kc.emplace_back(std::sqrt(lambda), i);
    }
    std::sort(indexed_kc.begin(), indexed_kc.end());

    // Avaliação de erro relativa com base nos autovalores esperados
    std::vector<int> node_flag;
    int num_free_nodes = 0;
    free_nodes(a, b, node_flag, num_free_nodes);
    int num_free_edges = std::count_if(edge_free_flag.begin(), edge_free_flag.end(), [](int f)
                                       { return f == 1; });
    int num_zero_eigvals = num_free_nodes + num_free_edges;

    std::vector<double> TEeigvalues;
    for (const auto& [kc, _] : indexed_kc) {
        TEeigvalues.push_back(kc);
    }
    // Avaliação de erro dos primeiros modos úteis
    std::vector<double> rel_err = TEeig_err(a, b, TEeigvalues, std::min(NUM_DOFS - num_zero_eigvals, 8), num_zero_eigvals);

    std::ofstream summary((out_dir / "eigdata_LTQN_modes.csv").string());
    summary << "rank,kc,internal_index,relative_error\n";
    summary << std::scientific << std::setprecision(10);
    for (std::size_t i = 0; i < indexed_kc.size(); ++i)
    {
        double err = -1.0;
        const int physical_rank = static_cast<int>(i) - num_zero_eigvals;
        if (physical_rank >= 0 && physical_rank < static_cast<int>(rel_err.size()))
            err = rel_err[physical_rank];
        summary << (i + 1) << ","
                << indexed_kc[i].first << ","
                << indexed_kc[i].second << ","
                << err << "\n";
    }
    summary.close();

    std::ofstream ferr((out_dir / "eigdata_LTQN_relerr.csv").string());
    ferr << "mode,rel_err\n";
    ferr << std::scientific << std::setprecision(10);
    for (std::size_t i = 0; i < rel_err.size(); ++i)
    {
        ferr << (i + 1) << "," << rel_err[i] << "\n";
    }
    ferr.close();

    std::ofstream meta((out_dir / "eigen2d_ltqn_summary.csv").string());
    meta << "key,value\n";
    meta << "a," << a << "\n";
    meta << "b," << b << "\n";
    meta << "x_mesh," << config.x_mesh << "\n";
    meta << "y_mesh," << config.y_mesh << "\n";
    meta << "num_nodes," << NUM_NODES << "\n";
    meta << "num_elements," << NUM_ELEMENTS << "\n";
    meta << "num_edges," << NUM_EDGES << "\n";
    meta << "num_dofs," << NUM_DOFS << "\n";
    meta << "num_zero_eigvals," << num_zero_eigvals << "\n";
    meta.close();

    std::cout << "Modos nulos/espurios ignorados: " << num_zero_eigvals << "\n";
    for (std::size_t i = 0; i < rel_err.size(); ++i)
    {
        std::cout << "Modo TE " << (i + 1)
                  << " | kc = " << TEeigvalues[num_zero_eigvals + static_cast<int>(i)]
                  << " | erro relativo = " << rel_err[i] << "\n";
    }
}

int main(int argc, char **argv)
{
    Eigen2DLtqnConfig config;
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
            if (arg == "--a")
                config.a = std::stod(require_value(i, argc, argv));
            else if (arg == "--b")
                config.b = std::stod(require_value(i, argc, argv));
            else if (arg == "--x-mesh")
                config.x_mesh = std::stoi(require_value(i, argc, argv));
            else if (arg == "--y-mesh")
                config.y_mesh = std::stoi(require_value(i, argc, argv));
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

    if (config.a <= 0.0 || config.b <= 0.0 || config.x_mesh <= 0 || config.y_mesh <= 0)
    {
        std::cerr << "Erro: parametros geometricos e malha devem ser positivos.\n";
        return 1;
    }

    Eigen2D_LTQN(config);
    return 0;
}
