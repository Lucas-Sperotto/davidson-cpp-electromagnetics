#include "outside_edge.hpp"

#include "mom3d_globals.hpp"

#include <limits>

namespace cap06 {

Eigen::VectorXi outside_edge(const double a, const double b)
{
    Eigen::VectorXi dof_flag = Eigen::VectorXi::Ones(NUM_EDGES);
    const double tol = std::numeric_limits<double>::epsilon();

    for (int i_edge = 0; i_edge < NUM_EDGES; ++i_edge) {
        const int node1 = EDGES(i_edge, 0);
        const int node2 = EDGES(i_edge, 1);

        if (std::abs(NODE_COORD(node1, 1)) < tol && std::abs(NODE_COORD(node2, 1)) < tol) {
            dof_flag(i_edge) = 0;
        }
        if (std::abs(NODE_COORD(node1, 1) - b) < tol && std::abs(NODE_COORD(node2, 1) - b) < tol) {
            dof_flag(i_edge) = 0;
        }
        if (std::abs(NODE_COORD(node1, 0)) < tol && std::abs(NODE_COORD(node2, 0)) < tol) {
            dof_flag(i_edge) = 0;
        }
        if (std::abs(NODE_COORD(node1, 0) - a) < tol && std::abs(NODE_COORD(node2, 0) - a) < tol) {
            dof_flag(i_edge) = 0;
        }
    }

    return dof_flag;
}

}  // namespace cap06
