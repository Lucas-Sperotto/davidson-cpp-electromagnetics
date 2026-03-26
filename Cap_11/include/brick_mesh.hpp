#ifndef CAP_11_BRICK_MESH_HPP
#define CAP_11_BRICK_MESH_HPP

#include <vector>

std::vector<std::vector<double>> brick_mesh(double x_len,
                                            double y_len,
                                            double z_len,
                                            int N_x,
                                            int N_y,
                                            int N_z);

#endif
