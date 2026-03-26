#include "ComputeRho_c.hpp"
#include "FillVVector.hpp"
#include "FillZMatrixByEdge.hpp"
#include "FillZMatrixByFace.hpp"
#include "PostProcMoM.hpp"
#include "edge_conx_elem.hpp"
#include "edgemake_MoM.hpp"
#include "find_local_dofs.hpp"
#include "mom3d_globals.hpp"
#include "outside_edge.hpp"
#include "renumber_RWG.hpp"
#include "trimesh3D.hpp"

#include <Eigen/Dense>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

int require_int(const std::string& flag, const char* value)
{
    if (value == nullptr) {
        throw std::runtime_error("Missing value for " + flag);
    }
    return std::stoi(value);
}

void print_help()
{
    std::cout
        << "Uso: ./build/MoM3D_demo [opcoes]\n"
        << "  --prob-type <5|6>  Seleciona o caso de teste de RWG82 (default: 6)\n"
        << "  --sing <0|1>       Usa o esquema de integral singular (default: 0)\n"
        << "  --quad-pts <n>     Numero de pontos da quadratura triangular (default: 6)\n"
        << "  --xmesh <n>        Sobrescreve a malha em x\n"
        << "  --ymesh <n>        Sobrescreve a malha em y\n"
        << "  --help             Mostra esta ajuda\n";
}

}  // namespace

int main(int argc, char** argv)
{
    using clock = std::chrono::steady_clock;

    cap06::initialize_mom3d_globals();
    std::filesystem::create_directories(cap06::project_out_dir());

    int ProbType = 6;
    bool sing = false;
    int quad_pts = 6;
    int xmesh_override = -1;
    int ymesh_override = -1;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--prob-type") {
            ProbType = require_int(arg, (i + 1 < argc) ? argv[++i] : nullptr);
        } else if (arg == "--sing") {
            sing = require_int(arg, (i + 1 < argc) ? argv[++i] : nullptr) != 0;
        } else if (arg == "--quad-pts") {
            quad_pts = require_int(arg, (i + 1 < argc) ? argv[++i] : nullptr);
        } else if (arg == "--xmesh") {
            xmesh_override = require_int(arg, (i + 1 < argc) ? argv[++i] : nullptr);
        } else if (arg == "--ymesh") {
            ymesh_override = require_int(arg, (i + 1 < argc) ? argv[++i] : nullptr);
        } else if (arg == "--help") {
            print_help();
            return 0;
        } else {
            throw std::runtime_error("Unknown argument: " + arg);
        }
    }

    const double eps_0 = 8.854e-12;
    const double mu_0 = 4.0 * cap06::PI * 1.0e-7;
    const double eta_0 = std::sqrt(mu_0 / eps_0);
    const double c = 1.0 / std::sqrt(eps_0 * mu_0);
    const double freq = c;
    const double omega = 2.0 * cap06::PI * freq;
    const double lambda = c / freq;
    const double k = 2.0 * cap06::PI / lambda;
    const double EMag = 10.0;
    const double theta_0 = 0.0;
    const double phi_0 = 0.0;

    double L = 0.0;
    double W = 0.0;
    int Xmesh = 0;
    int Ymesh = 0;

    switch (ProbType) {
    case 5:
        L = 0.15 * lambda;
        W = L;
        Xmesh = 6;
        Ymesh = 5;
        break;
    case 6:
        L = 1.0 * lambda;
        W = L;
        Xmesh = 6;
        Ymesh = 7;
        break;
    default:
        throw std::runtime_error("Unknown problem.");
    }

    if (xmesh_override > 0) {
        Xmesh = xmesh_override;
    }
    if (ymesh_override > 0) {
        Ymesh = ymesh_override;
    }

    const auto overall_start = clock::now();
    const auto prep_start = clock::now();

    cap06::trimesh3D(L, W, Xmesh, Ymesh);

    Eigen::MatrixXd r_c = Eigen::MatrixXd::Zero(cap06::NUM_ELEMENTS, 3);
    for (int ielem = 0; ielem < cap06::NUM_ELEMENTS; ++ielem) {
        r_c.row(ielem) = (
            cap06::NODE_COORD.row(cap06::ELEMENTS(ielem, 0)) +
            cap06::NODE_COORD.row(cap06::ELEMENTS(ielem, 1)) +
            cap06::NODE_COORD.row(cap06::ELEMENTS(ielem, 2))) /
            3.0;
    }

    cap06::edgemake_MoM();
    const Eigen::VectorXi dof_int = cap06::outside_edge(L, W);
    const auto [dof_RWG, dof2edge] = cap06::renumber_RWG(dof_int);
    cap06::edge_conx_elem(dof_RWG);
    cap06::find_local_dofs(dof_RWG);
    const cap06::RhoCenters rho = cap06::ComputeRho_c(r_c);

    const double time_preprocessing =
        std::chrono::duration<double>(clock::now() - prep_start).count();

    const auto fill_edge_start = clock::now();
    const Eigen::MatrixXcd Z = cap06::FillZMatrixByEdge(
        omega, eps_0, mu_0, k, r_c, rho.rho_c_pls, rho.rho_c_mns, quad_pts, sing, dof2edge);
    const double time_matrix_fill_by_edge =
        std::chrono::duration<double>(clock::now() - fill_edge_start).count();

    const auto fill_face_start = clock::now();
    const Eigen::MatrixXcd Z_1 = cap06::FillZMatrixByFace(
        omega, eps_0, mu_0, k, r_c, rho.rho_c_pls, rho.rho_c_mns, quad_pts, sing, dof2edge, dof_RWG);
    const double time_matrix_fill_by_face =
        std::chrono::duration<double>(clock::now() - fill_face_start).count();
    const double matrix_diff = (Z - Z_1).cwiseAbs().maxCoeff();

    const Eigen::VectorXcd V = cap06::FillVVector(rho.rho_c_pls, rho.rho_c_mns, EMag, theta_0, phi_0, dof2edge);

    const auto solve_start = clock::now();
    const Eigen::VectorXcd I = Z.partialPivLu().solve(V);
    const double time_matrix_solve = std::chrono::duration<double>(clock::now() - solve_start).count();

    const Eigen::VectorXcd I_1 = Z_1.partialPivLu().solve(V);
    const double current_diff = (I - I_1).cwiseAbs().maxCoeff();

    cap06::PostProcMoM(I, EMag, dof2edge, eta_0, L, W, Xmesh, Ymesh, ProbType, quad_pts, sing);
    const double time_overall = std::chrono::duration<double>(clock::now() - overall_start).count();

    {
        std::ofstream node_file(cap06::project_out_dir() / "mom3d_mesh_nodes.csv");
        node_file << "node,x_m,y_m,z_m\n";
        for (int inode = 0; inode < cap06::NUM_NODES; ++inode) {
            node_file << inode + 1 << ',' << cap06::NODE_COORD(inode, 0) << ',' << cap06::NODE_COORD(inode, 1)
                      << ',' << cap06::NODE_COORD(inode, 2) << '\n';
        }
    }

    {
        std::ofstream elem_file(cap06::project_out_dir() / "mom3d_mesh_elements.csv");
        elem_file << "element,node1,node2,node3,xc_m,yc_m,zc_m\n";
        for (int ielem = 0; ielem < cap06::NUM_ELEMENTS; ++ielem) {
            elem_file << ielem + 1 << ',' << cap06::ELEMENTS(ielem, 0) + 1 << ','
                      << cap06::ELEMENTS(ielem, 1) + 1 << ',' << cap06::ELEMENTS(ielem, 2) + 1 << ','
                      << r_c(ielem, 0) << ',' << r_c(ielem, 1) << ',' << r_c(ielem, 2) << '\n';
        }
    }

    {
        const Eigen::VectorXcd I_norm = eta_0 * I / EMag;
        std::ofstream current_file(cap06::project_out_dir() / "mom3d_currents.csv");
        current_file << "dof,edge,node1,node2,edge_center_x_m,edge_center_y_m,edge_center_z_m,abs_i,real_i,imag_i,abs_i_norm\n";
        for (int idof = 0; idof < cap06::NUM_DOFS; ++idof) {
            const int edge = dof2edge(idof);
            const int node1 = cap06::EDGES(edge, 0);
            const int node2 = cap06::EDGES(edge, 1);
            const Eigen::RowVector3d edge_center =
                0.5 * (cap06::NODE_COORD.row(node1) + cap06::NODE_COORD.row(node2));
            current_file << idof + 1 << ',' << edge + 1 << ',' << node1 + 1 << ',' << node2 + 1 << ','
                         << edge_center(0) << ',' << edge_center(1) << ',' << edge_center(2) << ','
                         << std::abs(I(idof)) << ',' << I(idof).real() << ',' << I(idof).imag() << ','
                         << std::abs(I_norm(idof)) << '\n';
        }
    }

    {
        std::ofstream summary_file(cap06::project_out_dir() / "mom3d_summary.csv");
        summary_file << "key,value\n";
        summary_file << "prob_type," << ProbType << '\n';
        summary_file << "sing," << static_cast<int>(sing) << '\n';
        summary_file << "quad_pts," << quad_pts << '\n';
        summary_file << "L_m," << L << '\n';
        summary_file << "W_m," << W << '\n';
        summary_file << "Xmesh," << Xmesh << '\n';
        summary_file << "Ymesh," << Ymesh << '\n';
        summary_file << "num_nodes," << cap06::NUM_NODES << '\n';
        summary_file << "num_elements," << cap06::NUM_ELEMENTS << '\n';
        summary_file << "num_edges," << cap06::NUM_EDGES << '\n';
        summary_file << "num_dofs," << cap06::NUM_DOFS << '\n';
        summary_file << "matrix_diff_max_abs," << matrix_diff << '\n';
        summary_file << "current_diff_max_abs," << current_diff << '\n';
        summary_file << "time_preprocessing_s," << time_preprocessing << '\n';
        summary_file << "time_matrix_fill_by_edge_s," << time_matrix_fill_by_edge << '\n';
        summary_file << "time_matrix_fill_by_face_s," << time_matrix_fill_by_face << '\n';
        summary_file << "time_matrix_solve_s," << time_matrix_solve << '\n';
        summary_file << "time_overall_s," << time_overall << '\n';
    }

    std::cout << "MoM3D_demo finalizado.\n";
    std::cout << "Problema: " << ProbType << ", DOFs: " << cap06::NUM_DOFS
              << ", diferenca maxima entre matrizes: " << matrix_diff << '\n';
    std::cout << "Diferenca maxima entre correntes: " << current_diff << '\n';
    return 0;
}
