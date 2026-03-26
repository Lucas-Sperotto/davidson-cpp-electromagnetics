#ifndef CAP_11_FREE_NODES3D_HPP
#define CAP_11_FREE_NODES3D_HPP

#include <utility>
#include <vector>

std::pair<std::vector<int>, int> free_nodes3D(double a,
                                              double b,
                                              double c,
                                              const std::vector<std::vector<double>> &node_coord,
                                              bool symmetry_flag,
                                              double toler);

#endif
