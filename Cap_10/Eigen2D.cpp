#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <fstream>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <string>
#include <filesystem>
#include <stdexcept>

#include "trimesh.h"
#include "edgemake.h"
#include "free_dof.h"
#include "renumber_dof.h"
#include "sandt.h"
#include "free_nodes.h"
#include "plot_field.h"
#include "simplex2D.h"
#include "whitney.h"
#include "TEeig_err.h"

using namespace Eigen;

#include "globals.h"

const double eps_0 = 8.854e-12;
const double mu_0 = 4.0 * M_PI * 1e-7;
const double eps_r = 1.0;
const double mu_r = 1.0;

namespace
{
struct Eigen2DConfig
{
    double a = 2.286e-2;
    double b = 1.016e-2;
    int x_mesh = 8;
    int y_mesh = 4;
    int field_grid_x = 21;
    int field_grid_y = 11;
    int num_physical_modes = 6;
    int num_spurious_modes = 6;
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
        << "Uso: ./build/Eigen2D [opcoes]\n"
        << "  --a VALOR                largura do guia (padrao 2.286e-2)\n"
        << "  --b VALOR                altura do guia (padrao 1.016e-2)\n"
        << "  --x-mesh N               divisao da malha em x (padrao 8)\n"
        << "  --y-mesh N               divisao da malha em y (padrao 4)\n"
        << "  --field-grid-x N         grade de avaliacao dos campos em x (padrao 21)\n"
        << "  --field-grid-y N         grade de avaliacao dos campos em y (padrao 11)\n"
        << "  --num-physical-modes N   numero de modos fisicos para exportar (padrao 6)\n"
        << "  --num-spurious-modes N   numero de modos espurios para exportar (padrao 6)\n"
        << "  --help                   mostra esta ajuda\n";
}
} // namespace

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
    {
            f << V(i) << '\n';
    }
}

void Eigen2D(const Eigen2DConfig &config)
{
    const std::filesystem::path out_dir = PROJECT_OUT_DIR;
    std::filesystem::create_directories(out_dir);

    const double a = config.a;
    const double b = config.b;

    std::cout << "Parametros do Eigen2D:\n"
              << "  a=" << a << ", b=" << b
              << ", x_mesh=" << config.x_mesh << ", y_mesh=" << config.y_mesh << "\n"
              << "  field_grid_x=" << config.field_grid_x << ", field_grid_y=" << config.field_grid_y
              << ", num_physical_modes=" << config.num_physical_modes
              << ", num_spurious_modes=" << config.num_spurious_modes << "\n";

    std::vector<std::vector<double>> x_nodes, y_nodes;
    trimesh(a, b, config.x_mesh, config.y_mesh, x_nodes, y_nodes);

    // for (int i = 0; i < NUM_NODES; ++i)
    //     std::cout << "Node " << i << ": (" << NODE_COORD[i][0] << ", " << NODE_COORD[i][1] << ")\n";

    for (int e = 0; e < NUM_ELEMENTS; ++e)
    {
        double xc = 0.0, yc = 0.0;
        for (int j = 0; j < 3; ++j)
        {
            xc += NODE_COORD[ELEMENTS[e][j]][0];
            yc += NODE_COORD[ELEMENTS[e][j]][1];
        }
        xc /= 3.0;
        yc /= 3.0;
        // std::cout << "Element " << e << " center: (" << xc << ", " << yc << ")\n";
    }

    // std::cout << "Press ENTER to continue...\n";
    // std::cin.get();

    // Geração das arestas a partir dos elementos
    edgemake(); // preenche EDGES, ELEMENT_EDGES, atualiza NUM_EDGES

    // Determina quais arestas estão em contorno PEC e devem ser prescritas
    std::vector<int> dof_e1_free_flag = free_dof(a, b);

    // Renumera apenas os DOFs livres
    std::vector<int> dof_e1 = renumber_dof(dof_e1_free_flag); // também atualiza NUM_DOFS

    // for (int j = 0; j < dof_e1.size(); ++j)
    //     std::cout << "dof_e1[" << j << "] = " << dof_e1[j] << "\n";

    // Inicializa matrizes globais
    std::vector<std::vector<double>> S(NUM_DOFS, std::vector<double>(NUM_DOFS, 0.0));
    std::vector<std::vector<double>> T(NUM_DOFS, std::vector<double>(NUM_DOFS, 0.0));

    for (int ielem = 0; ielem < NUM_ELEMENTS; ++ielem)
    {
        const auto &trinodes = ELEMENTS[ielem];
        // std::cout << "Tri[" << ielem << "] = (" << trinodes[0] << ", " << trinodes[1] << ", " << trinodes[2] << ")\n";

        for (int j = 0; j < 3; ++j)
        {
            if (trinodes[j] < 0 || trinodes[j] >= NODE_COORD.size())
            {
                std::cerr << "Erro: índice inválido no elemento " << ielem << ": trinodes[" << j << "] = "
                          << trinodes[j] << " (NODE_COORD tem tamanho " << NODE_COORD.size() << ")\n";
                exit(EXIT_FAILURE);
            }
        }

        double x1 = NODE_COORD[trinodes[0]][0], y1 = NODE_COORD[trinodes[0]][1];
        double x2 = NODE_COORD[trinodes[1]][0], y2 = NODE_COORD[trinodes[1]][1];
        double x3 = NODE_COORD[trinodes[2]][0], y3 = NODE_COORD[trinodes[2]][1];

        // std::cout << "Calling sandt with nodes: " << "(" << x1 << ", " << y1 << "), " << "(" << x2 << ", " << y2 << "), " << "(" << x3 << ", " << y3 << ")\n";

        auto [S_elem, T_elem] = sandt(x1, y1, x2, y2, x3, y3); // retorna 3x3
                                                               // std::cout << "Saiu de sandt para o elemento [" << ielem << "]\n";

        // for (int jedge = 0; jedge < 3; ++jedge)
        // {
        //    for (int kedge = 0; kedge < 3; ++kedge)
        //    {
        // std::cout << S_elem[jedge][kedge] << "\t";
        // std::cout << T_elem[jedge][kedge] << "\t";
        //    }
        // std::cout << std::endl;
        // }

        for (int jedge = 0; jedge < 3; ++jedge)
        {
            int edge_j = ELEMENT_EDGES[ielem][jedge];
            int jj = dof_e1[edge_j];
            for (int kedge = 0; kedge < 3; ++kedge)
            {
                int edge_k = ELEMENT_EDGES[ielem][kedge];
                int kk = dof_e1[edge_k];
                if (jj >= 0 && kk >= 0)
                {
                    S[jj][kk] += S_elem[jedge][kedge];
                    T[jj][kk] += T_elem[jedge][kedge];
                }
            }
        }
    }

    MatrixXd S_mat(NUM_DOFS, NUM_DOFS), T_mat(NUM_DOFS, NUM_DOFS);
    for (int i = 0; i < NUM_DOFS; ++i)
    {
        for (int j = 0; j < NUM_DOFS; ++j)
        {
            S_mat(i, j) = S[i][j] / mu_r;
            T_mat(i, j) = T[i][j] * eps_r;
        }
    }
    // std::cout << "S_mat dimensions: " << S_mat.rows() << " x " << S_mat.cols() << std::endl;
    /*
        // Resolver o problema generalizado: (S/mu_r) x = lambda (T*eps_r) x
        GeneralizedEigenSolver<MatrixXd> solver(S_mat, T_mat);

        // Autovalores (podem ser complexos)
        VectorXcd eigvals = solver.eigenvalues();

        // Tirar raiz quadrada de cada autovalor
        std::vector<std::complex<double>> sqrt_vals;
        sqrt_vals.reserve(eigvals.size());
        for (int i = 0; i < eigvals.size(); ++i) {
            sqrt_vals.push_back(std::sqrt(eigvals[i]));
        }

        // Ordenar (por parte real, como no MATLAB sort)
        std::sort(sqrt_vals.begin(), sqrt_vals.end(),
                  [](const std::complex<double>& a, const std::complex<double>& b) {
                      return a.real() < b.real();
                  });

        // Imprimir resultado
        std::cout << "Autovalores (sqrt ordenados):\n";
        for (auto& val : sqrt_vals) {
            std::cout << val << "\n";
        }


        */

    GeneralizedSelfAdjointEigenSolver<MatrixXd> solver(S_mat, T_mat);
    if (solver.info() != Eigen::Success)
    {
        std::cerr << "Falha ao resolver o problema generalizado de autovalores.\n";
        return;
    }

    write_csv((out_dir / "cpp_S_mat.csv").string(), S_mat);
    write_csv((out_dir / "cpp_T_mat.csv").string(), T_mat);

    VectorXd eigvals = solver.eigenvalues();
    MatrixXd eigvecs = solver.eigenvectors();

    write_csv((out_dir / "cpp_eigvals.csv").string(), eigvals);
    write_csv((out_dir / "cpp_eigvecs.csv").string(), eigvecs);

    MatrixXd V = eigvecs; // cópia
                          // for (int i = 0; i < V.cols(); ++i) {
    //     double nB = std::sqrt(V.col(i).transpose() * T_mat * V.col(i));
    //     if (nB > 0) V.col(i) /= nB;
    // }
    //  for (int i = 0; i < eigvals.size(); ++i)
    //      std::cout << "Autovalor[" << i << "]: " << eigvals[i] << std::endl;
    //  for (int i = 0; i < eigvals.size(); ++i)
    //     std::cout << "eig[" << i << "]: " << eigvals[i] << std::endl;

    constexpr double zero_tol = 1e-6;

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

    std::vector<int> node_flag;
    int num_free_nodes = 0;
    free_nodes(a, b, node_flag, num_free_nodes);

    int start_idx = num_free_nodes;
    int available_modes = indexed_kc.size() > static_cast<std::size_t>(start_idx)
                              ? static_cast<int>(indexed_kc.size()) - start_idx
                              : 0;
    int num_plot_modes = std::min(config.num_physical_modes, available_modes);
    if (num_plot_modes <= 0)
    {
        std::cerr << "Nao ha modos fisicos suficientes para plotar.\n";
        return;
    }

    std::cout << "Modos nulos/espurios ignorados: " << start_idx << "\n";
    // Geração da grade de avaliação
    std::vector<double> XX, YY;
    int NX = config.field_grid_x;
    int NY = config.field_grid_y;
    for (int i = 0; i <= NX; ++i)
        XX.push_back((double)i / NX * a);
    for (int j = 0; j <= NY; ++j)
        YY.push_back((double)j / NY * b);

    std::vector<double> TEeigvalues;
    TEeigvalues.reserve(indexed_kc.size());
    for (const auto &[kc, _] : indexed_kc)
        TEeigvalues.push_back(kc);

    std::vector<double> rel_err =
        TEeig_err(a, b, TEeigvalues, std::min(available_modes, 8), num_free_nodes);

    std::ofstream summary((out_dir / "eigdata_whitney_modes.csv").string());
    summary << "rank,kind,kc,internal_index,relative_error\n";
    summary << std::scientific << std::setprecision(10);
    for (std::size_t i = 0; i < indexed_kc.size(); ++i)
    {
        const bool is_physical = static_cast<int>(i) >= start_idx;
        const int physical_rank = static_cast<int>(i) - start_idx;
        double err = -1.0;
        if (is_physical && physical_rank >= 0 && physical_rank < static_cast<int>(rel_err.size()))
            err = rel_err[physical_rank];
        summary << (i + 1) << ","
                << (is_physical ? "physical" : "spurious") << ","
                << indexed_kc[i].first << ","
                << indexed_kc[i].second << ","
                << err << "\n";
    }
    summary.close();

    std::ofstream meta((out_dir / "eigen2d_summary.csv").string());
    meta << "key,value\n";
    meta << "a," << a << "\n";
    meta << "b," << b << "\n";
    meta << "x_mesh," << config.x_mesh << "\n";
    meta << "y_mesh," << config.y_mesh << "\n";
    meta << "num_nodes," << NUM_NODES << "\n";
    meta << "num_elements," << NUM_ELEMENTS << "\n";
    meta << "num_edges," << NUM_EDGES << "\n";
    meta << "num_dofs," << NUM_DOFS << "\n";
    meta << "num_free_nodes," << num_free_nodes << "\n";
    meta << "num_physical_modes_exported," << num_plot_modes << "\n";
    meta.close();

    for (int ii = 0; ii < num_plot_modes; ++ii)
    {
        int idx = indexed_kc[start_idx + ii].second;

        VectorXd eigmode = eigvecs.col(idx);

        std::cout << "Modo TE " << (ii + 1)
                  << " | indice interno " << idx
                  << " | kc = " << indexed_kc[start_idx + ii].first << "\n";

        std::vector<double> dofs(eigmode.data(), eigmode.data() + eigmode.size());

        plot_field(dofs, dof_e1, XX, YY, 3, 2, ii + 1, indexed_kc[start_idx + ii].first, "field_kc_");
    }

    const int num_spurious_modes = std::min(config.num_spurious_modes, start_idx);
    for (int ii = 0; ii < num_spurious_modes; ++ii)
    {
        int idx = indexed_kc[ii].second;
        VectorXd eigmode = eigvecs.col(idx);
        std::vector<double> dofs(eigmode.data(), eigmode.data() + eigmode.size());
        std::ostringstream prefix;
        prefix << "spur_field_mode_" << std::setw(2) << std::setfill('0') << (ii + 1) << "_kc_";
        plot_field(dofs, dof_e1, XX, YY, 3, 2, ii + 1, indexed_kc[ii].first, prefix.str());
    }
}

int main(int argc, char **argv)
{
    Eigen2DConfig config;
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
            else if (arg == "--field-grid-x")
                config.field_grid_x = std::stoi(require_value(i, argc, argv));
            else if (arg == "--field-grid-y")
                config.field_grid_y = std::stoi(require_value(i, argc, argv));
            else if (arg == "--num-physical-modes")
                config.num_physical_modes = std::stoi(require_value(i, argc, argv));
            else if (arg == "--num-spurious-modes")
                config.num_spurious_modes = std::stoi(require_value(i, argc, argv));
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

    if (config.a <= 0.0 || config.b <= 0.0 || config.x_mesh <= 0 || config.y_mesh <= 0 ||
        config.field_grid_x <= 0 || config.field_grid_y <= 0 ||
        config.num_physical_modes < 0 || config.num_spurious_modes < 0)
    {
        std::cerr << "Erro: parametros geometricos, malha e contagens devem ser validos.\n";
        return 1;
    }

    Eigen2D(config);
    return 0;
}
