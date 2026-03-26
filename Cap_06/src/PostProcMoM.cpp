#include "PostProcMoM.hpp"

#include "mom3d_globals.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <tuple>
#include <vector>

namespace cap06 {

void PostProcMoM(
    const Eigen::VectorXcd& I,
    const double EMag,
    const Eigen::VectorXi& dof2edge,
    const double eta_0,
    const double L,
    const double W,
    const int Xmesh,
    const int Ymesh,
    const int ProbType,
    const int quad_pts,
    const bool sing)
{
    const Eigen::VectorXcd I_norm = eta_0 * I / EMag;

    Eigen::MatrixXd ell_vec = Eigen::MatrixXd::Zero(NUM_EDGES, 3);
    for (int iedge = 0; iedge < NUM_EDGES; ++iedge) {
        ell_vec.row(iedge) = NODE_COORD.row(EDGES(iedge, 1)) - NODE_COORD.row(EDGES(iedge, 0));
        ell_vec.row(iedge) /= ELL(iedge);
    }

    const double tol = ELL.maxCoeff() / 1.0e6;
    std::vector<std::pair<double, double>> cut_AA;
    std::vector<std::pair<double, double>> cut_BB;

    if ((Xmesh % 2) == 0) {
        for (int idof = 0; idof < NUM_DOFS; ++idof) {
            const int edge = dof2edge(idof);
            if (std::abs(ell_vec(edge, 1) - 1.0) < tol &&
                std::abs(NODE_COORD(EDGES(edge, 0), 0) - L / 2.0) < tol) {
                const double y_vert = 0.5 * (NODE_COORD(EDGES(edge, 0), 1) + NODE_COORD(EDGES(edge, 1), 1));
                cut_AA.emplace_back(y_vert, std::abs(I_norm(idof)));
            }
        }
    } else {
        std::cerr << "Aviso: Xmesh deve ser par para o corte AA.\n";
    }

    if ((Ymesh % 2) != 0) {
        for (int idof = 0; idof < NUM_DOFS; ++idof) {
            const int edge = dof2edge(idof);
            const double y_c = 0.5 * (NODE_COORD(EDGES(edge, 0), 1) + NODE_COORD(EDGES(edge, 1), 1));
            if (std::abs(ell_vec(edge, 1) - 1.0) < tol && std::abs(y_c - L / 2.0) < tol) {
                const double x_hor = 0.5 * (NODE_COORD(EDGES(edge, 0), 0) + NODE_COORD(EDGES(edge, 1), 0));
                cut_BB.emplace_back(x_hor, std::abs(I_norm(idof)));
            }
        }
    } else {
        std::cerr << "Aviso: Ymesh deve ser impar para o corte BB.\n";
    }

    auto sorter = [](const auto& lhs, const auto& rhs) { return lhs.first < rhs.first; };
    std::sort(cut_AA.begin(), cut_AA.end(), sorter);
    std::sort(cut_BB.begin(), cut_BB.end(), sorter);

    std::vector<std::pair<double, double>> ref_AA;
    std::vector<std::pair<double, double>> ref_BB;
    switch (ProbType) {
    case 5:
        ref_AA = {
            {L * 1.0 / 10.0, 1.8},
            {L * 3.0 / 10.0, 1.1},
            {L * 5.0 / 10.0, 1.1},
            {L * 7.0 / 10.0, 1.1},
            {L * 9.0 / 10.0, 1.8},
        };
        ref_BB = {
            {W * 1.0 / 6.0, 0.8},
            {W * 2.0 / 6.0, 1.0},
            {W * 3.0 / 6.0, 1.1},
            {W * 4.0 / 6.0, 1.0},
            {W * 5.0 / 6.0, 0.8},
        };
        break;
    case 6:
        ref_AA = {
            {L * 1.0 / 14.0, 4.6},
            {L * 3.0 / 14.0, 2.6},
            {L * 5.0 / 14.0, 2.8},
            {L * 7.0 / 14.0, 2.9},
            {L * 9.0 / 14.0, 2.8},
            {L * 11.0 / 14.0, 2.6},
            {L * 13.0 / 14.0, 4.6},
        };
        ref_BB = {
            {W * 1.0 / 6.0, 1.7},
            {W * 2.0 / 6.0, 2.5},
            {W * 3.0 / 6.0, 2.9},
            {W * 4.0 / 6.0, 2.5},
            {W * 5.0 / 6.0, 1.7},
        };
        break;
    default:
        break;
    }

    std::ofstream cut_file(project_out_dir() / "mom3d_cuts.csv");
    cut_file << "series,coordinate_m,value\n";
    for (const auto& [coord, value] : cut_AA) {
        cut_file << "AA_computed," << coord << ',' << value << '\n';
    }
    for (const auto& [coord, value] : ref_AA) {
        cut_file << "AA_reference," << coord << ',' << value << '\n';
    }
    for (const auto& [coord, value] : cut_BB) {
        cut_file << "BB_computed," << coord << ',' << value << '\n';
    }
    for (const auto& [coord, value] : ref_BB) {
        cut_file << "BB_reference," << coord << ',' << value << '\n';
    }

    std::ofstream metadata_file(project_out_dir() / "mom3d_postproc_summary.csv");
    metadata_file << "key,value\n";
    metadata_file << "prob_type," << ProbType << '\n';
    metadata_file << "quad_pts," << quad_pts << '\n';
    metadata_file << "sing," << static_cast<int>(sing) << '\n';
    metadata_file << "num_cut_aa," << cut_AA.size() << '\n';
    metadata_file << "num_cut_bb," << cut_BB.size() << '\n';
}

}  // namespace cap06
