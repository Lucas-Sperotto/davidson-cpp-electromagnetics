#include "renumber_RWG.hpp"

#include "mom3d_globals.hpp"

#include <vector>

namespace cap06 {

std::pair<Eigen::VectorXi, Eigen::VectorXi> renumber_RWG(const Eigen::VectorXi& dof_free_flag)
{
    Eigen::VectorXi dof_num = Eigen::VectorXi::Zero(NUM_EDGES);
    std::vector<int> dof2edge_vec;

    int counter = 0;
    for (int i_edge = 0; i_edge < NUM_EDGES; ++i_edge) {
        if (dof_free_flag(i_edge) != 0) {
            ++counter;
            dof_num(i_edge) = counter;
            dof2edge_vec.push_back(i_edge);
        }
    }

    NUM_DOFS = counter;
    Eigen::VectorXi dof2edge(NUM_DOFS);
    for (int idof = 0; idof < NUM_DOFS; ++idof) {
        dof2edge(idof) = dof2edge_vec[idof];
    }

    return {dof_num, dof2edge};
}

}  // namespace cap06
