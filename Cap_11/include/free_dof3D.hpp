#ifndef CAP_11_FREE_DOF3D_HPP
#define CAP_11_FREE_DOF3D_HPP

#include <vector>

std::vector<int> free_dof3D(double a,
                            double b,
                            double c,
                            const std::vector<std::vector<double>> &node_coord,
                            bool symmetry_flag,
                            double toler);

#endif
