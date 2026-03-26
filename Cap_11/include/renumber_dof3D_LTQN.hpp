#ifndef CAP_11_RENUMBER_DOF3D_LTQN_HPP
#define CAP_11_RENUMBER_DOF3D_LTQN_HPP

#include <vector>

struct RenumberLtqn3DResult
{
    std::vector<int> dof_e1;
    std::vector<int> dof_e2;
    std::vector<int> dof_f1;
    std::vector<int> dof_f2;
    int num_free_edges = 0;
    int num_free_faces = 0;
    int num_dofs = 0;
};

RenumberLtqn3DResult renumber_dof3D_LTQN(const std::vector<int> &edge_free_flag,
                                         const std::vector<int> &face_free_flag);

#endif
