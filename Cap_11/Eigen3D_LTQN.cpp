#include "include/globals3d.hpp"
#include "include/brick_mesh.hpp"
#include "include/internal_tet_mesh.hpp"
#include "include/read_gmsh2.hpp"
#include "include/edgemake3D.hpp"
#include "include/facemake3D.hpp"
#include "include/free_dof3D_LTQN.hpp"
#include "include/renumber_dof3D.hpp"
#include "include/renumber_dof3D_LTQN.hpp"
#include "include/sandt3D.hpp"
#include "include/sandt3D_LTQN.hpp"
#include "include/free_nodes3D.hpp"
#include "include/eig_err3D.hpp"
#include "include/avg_mesh_length.hpp"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Eigen;
namespace fs = std::filesystem;

namespace
{
struct Config
{
    int element_order = 1;
    bool internal_mesh = true;
    std::string mesh_file;
    bool symmetry_flag = false;
    bool eigenvalues_only = true;
    double L_x = 1.0;
    double L_y = 0.5;
    double L_z = 0.75;
    int N_x = 4;
    int N_y = 2;
    int N_z = 3;
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
        << "Uso: ./build/Eigen3D_LTQN [opcoes]\n"
        << "  --element-order 1|2   1 para CTLN, 2 para LTQN (padrao 1)\n"
        << "  --internal-mesh 0|1   gera malha internamente (padrao 1)\n"
        << "  --mesh-file CAMINHO   arquivo .msh quando --internal-mesh 0\n"
        << "  --symmetry 0|1        usa planos de simetria magnetica (padrao 0)\n"
        << "  --eigenvalues-only 0|1 calcula apenas autovalores (padrao 1)\n"
        << "  --lx VALOR            dimensao x da cavidade (padrao 1.0)\n"
        << "  --ly VALOR            dimensao y da cavidade (padrao 0.5)\n"
        << "  --lz VALOR            dimensao z da cavidade (padrao 0.75)\n"
        << "  --nx N                bricks em x para malha interna (padrao 4)\n"
        << "  --ny N                bricks em y para malha interna (padrao 2)\n"
        << "  --nz N                bricks em z para malha interna (padrao 3)\n"
        << "  --help                mostra esta ajuda\n";
}

std::string default_mesh_file()
{
    return "../original_matlab/Chapter 11/Gmsh files/box_20.msh";
}

std::vector<double> sorted_kc_from_values(const VectorXd &values)
{
    constexpr double zero_tol = 1e-9;
    std::vector<double> eigvalues;
    eigvalues.reserve(values.size());
    for (int index = 0; index < values.size(); ++index)
    {
        double lambda = values[index];
        if (std::abs(lambda) <= zero_tol)
            lambda = 0.0;
        if (lambda < 0.0 && std::abs(lambda) <= 1e-7)
            lambda = 0.0;
        if (lambda >= 0.0)
            eigvalues.push_back(std::sqrt(lambda));
    }
    std::sort(eigvalues.begin(), eigvalues.end());
    return eigvalues;
}
}

int main(int argc, char **argv)
{
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
            if (arg == "--element-order")
                config.element_order = std::stoi(require_value(i, argc, argv));
            else if (arg == "--internal-mesh")
                config.internal_mesh = (std::stoi(require_value(i, argc, argv)) != 0);
            else if (arg == "--mesh-file")
                config.mesh_file = require_value(i, argc, argv);
            else if (arg == "--symmetry")
                config.symmetry_flag = (std::stoi(require_value(i, argc, argv)) != 0);
            else if (arg == "--eigenvalues-only")
                config.eigenvalues_only = (std::stoi(require_value(i, argc, argv)) != 0);
            else if (arg == "--lx")
                config.L_x = std::stod(require_value(i, argc, argv));
            else if (arg == "--ly")
                config.L_y = std::stod(require_value(i, argc, argv));
            else if (arg == "--lz")
                config.L_z = std::stod(require_value(i, argc, argv));
            else if (arg == "--nx")
                config.N_x = std::stoi(require_value(i, argc, argv));
            else if (arg == "--ny")
                config.N_y = std::stoi(require_value(i, argc, argv));
            else if (arg == "--nz")
                config.N_z = std::stoi(require_value(i, argc, argv));
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

    if (config.element_order != 1 && config.element_order != 2)
    {
        std::cerr << "Erro: --element-order deve ser 1 ou 2.\n";
        return 1;
    }

    initialize_local_topology3d();
    reset_global_mesh3d();

    const fs::path out_dir = PROJECT_OUT_DIR;
    fs::create_directories(out_dir);

    if (config.internal_mesh)
    {
        NODE_COORD = brick_mesh(config.symmetry_flag ? config.L_x / 2.0 : config.L_x,
                                config.symmetry_flag ? config.L_y / 2.0 : config.L_y,
                                config.L_z,
                                config.N_x, config.N_y, config.N_z);
        ELEMENTS = split_brick_mesh_to_tets(config.N_x, config.N_y, config.N_z);
    }
    else
    {
        if (config.mesh_file.empty())
            config.mesh_file = default_mesh_file();
        const GmshReadResult mesh = read_gmsh2(config.mesh_file);
        NODE_COORD = mesh.nodes;
        ELEMENTS = mesh.tets;
    }
    for (auto &tet : ELEMENTS)
        std::sort(tet.begin(), tet.end());
    NUM_NODES = static_cast<int>(NODE_COORD.size());
    NUM_ELEMENTS = static_cast<int>(ELEMENTS.size());

    edgemake3D();
    facemake3D();
    const double toler = std::min({config.L_x, config.L_y, config.L_z}) * 1e-5;
    const auto [edge_free_flag, face_free_flag] =
        free_dof3D_LTQN(config.L_x, config.L_y, config.L_z, NODE_COORD,
                        config.symmetry_flag, toler);

    std::vector<double> eigvalues;
    int num_free_edges = 0;
    int num_free_faces = 0;

    if (config.element_order == 1)
    {
        const auto [dof_e1, free_edges] = renumber_dof3D(edge_free_flag);
        num_free_edges = free_edges;
        NUM_DOFS = num_free_edges;

        MatrixXd S = MatrixXd::Zero(NUM_DOFS, NUM_DOFS);
        MatrixXd T = MatrixXd::Zero(NUM_DOFS, NUM_DOFS);
        for (int ielem = 0; ielem < NUM_ELEMENTS; ++ielem)
        {
            const auto [S_elem, T_elem] =
                sandt3D(NODE_COORD[ELEMENTS[ielem][0]],
                        NODE_COORD[ELEMENTS[ielem][1]],
                        NODE_COORD[ELEMENTS[ielem][2]],
                        NODE_COORD[ELEMENTS[ielem][3]]);

            for (int jedge = 0; jedge < 6; ++jedge)
            {
                const int jj = dof_e1[ELEMENT_EDGES[ielem][jedge]];
                for (int kedge = 0; kedge < 6; ++kedge)
                {
                    const int kk = dof_e1[ELEMENT_EDGES[ielem][kedge]];
                    if (jj >= 0 && kk >= 0)
                    {
                        S(jj, kk) += S_elem(jedge, kedge);
                        T(jj, kk) += T_elem(jedge, kedge);
                    }
                }
            }
        }

        if (config.eigenvalues_only)
        {
            GeneralizedSelfAdjointEigenSolver<MatrixXd> solver;
            solver.compute(S, T, EigenvaluesOnly);
            if (solver.info() != Success)
            {
                std::cerr << "Falha ao resolver o autoproblema CTLN em Eigen3D_LTQN.\n";
                return 1;
            }
            eigvalues = sorted_kc_from_values(solver.eigenvalues());
        }
        else
        {
            GeneralizedSelfAdjointEigenSolver<MatrixXd> solver(S, T);
            if (solver.info() != Success)
            {
                std::cerr << "Falha ao resolver o autoproblema CTLN em Eigen3D_LTQN.\n";
                return 1;
            }
            eigvalues = sorted_kc_from_values(solver.eigenvalues());
        }
    }
    else
    {
        const RenumberLtqn3DResult dofs =
            renumber_dof3D_LTQN(edge_free_flag, face_free_flag);
        num_free_edges = dofs.num_free_edges;
        num_free_faces = dofs.num_free_faces;
        NUM_DOFS = dofs.num_dofs;

        MatrixXd S = MatrixXd::Zero(NUM_DOFS, NUM_DOFS);
        MatrixXd T = MatrixXd::Zero(NUM_DOFS, NUM_DOFS);
        for (int ielem = 0; ielem < NUM_ELEMENTS; ++ielem)
        {
            const auto [S_elem, T_elem] =
                sandt3D_LTQN(NODE_COORD[ELEMENTS[ielem][0]],
                             NODE_COORD[ELEMENTS[ielem][1]],
                             NODE_COORD[ELEMENTS[ielem][2]],
                             NODE_COORD[ELEMENTS[ielem][3]]);

            for (int jedge = 0; jedge < 6; ++jedge)
            {
                const int jj_e1 = dofs.dof_e1[ELEMENT_EDGES[ielem][jedge]];
                const int jj_e2 = dofs.dof_e2[ELEMENT_EDGES[ielem][jedge]];
                for (int kedge = 0; kedge < 6; ++kedge)
                {
                    const int kk_e1 = dofs.dof_e1[ELEMENT_EDGES[ielem][kedge]];
                    const int kk_e2 = dofs.dof_e2[ELEMENT_EDGES[ielem][kedge]];
                    if (jj_e1 >= 0 && kk_e1 >= 0)
                    {
                        S(jj_e1, kk_e1) += S_elem(jedge, kedge);
                        T(jj_e1, kk_e1) += T_elem(jedge, kedge);
                        S(jj_e2, kk_e1) += S_elem(jedge + 6, kedge);
                        T(jj_e2, kk_e1) += T_elem(jedge + 6, kedge);
                        S(jj_e2, kk_e2) += S_elem(jedge + 6, kedge + 6);
                        T(jj_e2, kk_e2) += T_elem(jedge + 6, kedge + 6);
                        S(kk_e1, jj_e2) = S(jj_e2, kk_e1);
                        T(kk_e1, jj_e2) = T(jj_e2, kk_e1);
                    }
                }
            }

            for (int jface = 0; jface < 4; ++jface)
            {
                const int jj_f1 = dofs.dof_f1[ELEMENT_FACES[ielem][jface]];
                const int jj_f2 = dofs.dof_f2[ELEMENT_FACES[ielem][jface]];
                for (int kedge = 0; kedge < 6; ++kedge)
                {
                    const int kk_e1 = dofs.dof_e1[ELEMENT_EDGES[ielem][kedge]];
                    const int kk_e2 = dofs.dof_e2[ELEMENT_EDGES[ielem][kedge]];
                    if (jj_f1 >= 0 && kk_e1 >= 0)
                    {
                        S(jj_f1, kk_e1) += S_elem(jface + 12, kedge);
                        T(jj_f1, kk_e1) += T_elem(jface + 12, kedge);
                        S(jj_f1, kk_e2) += S_elem(jface + 12, kedge + 6);
                        T(jj_f1, kk_e2) += T_elem(jface + 12, kedge + 6);
                        S(jj_f2, kk_e1) += S_elem(jface + 16, kedge);
                        T(jj_f2, kk_e1) += T_elem(jface + 16, kedge);
                        S(jj_f2, kk_e2) += S_elem(jface + 16, kedge + 6);
                        T(jj_f2, kk_e2) += T_elem(jface + 16, kedge + 6);

                        S(kk_e1, jj_f1) = S(jj_f1, kk_e1);
                        T(kk_e1, jj_f1) = T(jj_f1, kk_e1);
                        S(kk_e2, jj_f1) = S(jj_f1, kk_e2);
                        T(kk_e2, jj_f1) = T(jj_f1, kk_e2);
                        S(kk_e1, jj_f2) = S(jj_f2, kk_e1);
                        T(kk_e1, jj_f2) = T(jj_f2, kk_e1);
                        S(kk_e2, jj_f2) = S(jj_f2, kk_e2);
                        T(kk_e2, jj_f2) = T(jj_f2, kk_e2);
                    }
                }

                for (int kface = 0; kface < 4; ++kface)
                {
                    const int kk_f1 = dofs.dof_f1[ELEMENT_FACES[ielem][kface]];
                    const int kk_f2 = dofs.dof_f2[ELEMENT_FACES[ielem][kface]];
                    if (jj_f1 >= 0 && kk_f1 >= 0)
                    {
                        S(jj_f1, kk_f1) += S_elem(jface + 12, kface + 12);
                        T(jj_f1, kk_f1) += T_elem(jface + 12, kface + 12);
                        S(jj_f2, kk_f1) += S_elem(jface + 16, kface + 12);
                        T(jj_f2, kk_f1) += T_elem(jface + 16, kface + 12);
                        S(jj_f2, kk_f2) += S_elem(jface + 16, kface + 16);
                        T(jj_f2, kk_f2) += T_elem(jface + 16, kface + 16);
                        S(kk_f1, jj_f2) = S(jj_f2, kk_f1);
                        T(kk_f1, jj_f2) = T(jj_f2, kk_f1);
                    }
                }
            }
        }

        if (config.eigenvalues_only)
        {
            GeneralizedSelfAdjointEigenSolver<MatrixXd> solver;
            solver.compute(S, T, EigenvaluesOnly);
            if (solver.info() != Success)
            {
                std::cerr << "Falha ao resolver o autoproblema LTQN.\n";
                return 1;
            }
            eigvalues = sorted_kc_from_values(solver.eigenvalues());
        }
        else
        {
            GeneralizedSelfAdjointEigenSolver<MatrixXd> solver(S, T);
            if (solver.info() != Success)
            {
                std::cerr << "Falha ao resolver o autoproblema LTQN.\n";
                return 1;
            }
            eigvalues = sorted_kc_from_values(solver.eigenvalues());
        }
    }

    const auto [node_flag, num_free_nodes] =
        free_nodes3D(config.L_x, config.L_y, config.L_z, NODE_COORD,
                     config.symmetry_flag, toler);
    (void)node_flag;
    const int num_zero_eig =
        (config.element_order == 1) ? num_free_nodes : (num_free_nodes + num_free_edges);
    const int usable_modes = std::min(NUM_DOFS - num_zero_eig, 8);
    const std::vector<double> rel_err =
        eig_err3D(config.L_x, config.L_y, config.L_z, eigvalues, usable_modes, num_zero_eig);
    const double h = avg_mesh_length(NODE_COORD);

    std::ofstream modes(out_dir / "eigdata3D_ltqn_modes.csv");
    modes << "rank,kc,relative_error\n";
    modes << std::scientific << std::setprecision(10);
    for (std::size_t index = 0; index < eigvalues.size(); ++index)
    {
        double err = -1.0;
        const int physical_index = static_cast<int>(index) - num_zero_eig;
        if (physical_index >= 0 && physical_index < static_cast<int>(rel_err.size()))
            err = rel_err[physical_index];
        modes << (index + 1) << "," << eigvalues[index] << "," << err << "\n";
    }

    std::ofstream summary(out_dir / "eigen3d_ltqn_summary.csv");
    summary << "key,value\n";
    summary << "element_order," << config.element_order << "\n";
    summary << "internal_mesh," << (config.internal_mesh ? 1 : 0) << "\n";
    summary << "mesh_file," << (config.mesh_file.empty() ? "internal" : config.mesh_file) << "\n";
    summary << "symmetry_flag," << (config.symmetry_flag ? 1 : 0) << "\n";
    summary << "eigenvalues_only," << (config.eigenvalues_only ? 1 : 0) << "\n";
    summary << "L_x," << config.L_x << "\n";
    summary << "L_y," << config.L_y << "\n";
    summary << "L_z," << config.L_z << "\n";
    summary << "N_x," << config.N_x << "\n";
    summary << "N_y," << config.N_y << "\n";
    summary << "N_z," << config.N_z << "\n";
    summary << "num_nodes," << NUM_NODES << "\n";
    summary << "num_elements," << NUM_ELEMENTS << "\n";
    summary << "num_edges," << NUM_EDGES << "\n";
    summary << "num_faces," << NUM_FACES << "\n";
    summary << "num_dofs," << NUM_DOFS << "\n";
    summary << "num_free_edges," << num_free_edges << "\n";
    summary << "num_free_faces," << num_free_faces << "\n";
    summary << "num_zero_eig," << num_zero_eig << "\n";
    summary << "avg_mesh_length," << h << "\n";

    std::cout << "Arquivos gerados: eigdata3D_ltqn_modes.csv e eigen3d_ltqn_summary.csv\n";
    return 0;
}
