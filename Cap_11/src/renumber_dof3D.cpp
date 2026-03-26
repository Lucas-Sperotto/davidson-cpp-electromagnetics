#include "../include/renumber_dof3D.hpp"
#include "../include/globals3d.hpp"

#include <utility>
#include <vector>

// Translation of renumber_dof.m:
// renumber free edge DOFs for CTLN/Whitney elements.
std::pair<std::vector<int>, int> renumber_dof3D(const std::vector<int> &edge_free_flag)
{
    std::vector<int> dof_e1(NUM_EDGES, -1);
    int counter = 0;
    for (int i_edge = 0; i_edge < NUM_EDGES; ++i_edge)
    {
        if (edge_free_flag[i_edge])
            dof_e1[i_edge] = counter++;
    }
    return {dof_e1, counter};
}
