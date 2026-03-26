#ifndef CAP_11_FREE_DOF3D_LTQN_HPP
#define CAP_11_FREE_DOF3D_LTQN_HPP

#include <utility>
#include <vector>

std::pair<std::vector<int>, std::vector<int>> free_dof3D_LTQN(
    double a,
    double b,
    double c,
    const std::vector<std::vector<double>> &node_coord,
    bool symmetry_flag,
    double toler);

#endif
