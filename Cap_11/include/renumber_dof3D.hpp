#ifndef CAP_11_RENUMBER_DOF3D_HPP
#define CAP_11_RENUMBER_DOF3D_HPP

#include <utility>
#include <vector>

std::pair<std::vector<int>, int> renumber_dof3D(const std::vector<int> &edge_free_flag);

#endif
