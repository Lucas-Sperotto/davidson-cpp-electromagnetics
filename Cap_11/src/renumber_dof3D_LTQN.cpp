#include "../include/renumber_dof3D_LTQN.hpp"
#include "../include/globals3d.hpp"

#include <vector>

// Translation of renumber_dof_LTQN.m:
// renumber the mixed 2nd-order edge and face DOFs.
RenumberLtqn3DResult renumber_dof3D_LTQN(const std::vector<int> &edge_free_flag,
                                         const std::vector<int> &face_free_flag)
{
    RenumberLtqn3DResult result;
    result.dof_e1.assign(NUM_EDGES, -1);
    result.dof_e2.assign(NUM_EDGES, -1);
    result.dof_f1.assign(NUM_FACES, -1);
    result.dof_f2.assign(NUM_FACES, -1);

    int counter = 0;
    for (int i_edge = 0; i_edge < NUM_EDGES; ++i_edge)
    {
        if (edge_free_flag[i_edge])
        {
            ++result.num_free_edges;
            result.dof_e1[i_edge] = counter++;
            result.dof_e2[i_edge] = counter++;
        }
    }

    for (int i_face = 0; i_face < NUM_FACES; ++i_face)
    {
        if (face_free_flag[i_face])
        {
            ++result.num_free_faces;
            result.dof_f1[i_face] = counter++;
            result.dof_f2[i_face] = counter++;
        }
    }

    result.num_dofs = counter;
    return result;
}
