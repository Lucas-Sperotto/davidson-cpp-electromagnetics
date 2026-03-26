#include "../include/brick_mesh.hpp"

#include <vector>

// Translation of brick_mesh.m:
// generate the regular brick-grid vertices in column-major order.
std::vector<std::vector<double>> brick_mesh(double x_len,
                                            double y_len,
                                            double z_len,
                                            int N_x,
                                            int N_y,
                                            int N_z)
{
    const double delta_x = x_len / static_cast<double>(N_x);
    const double delta_y = y_len / static_cast<double>(N_y);
    const double delta_z = z_len / static_cast<double>(N_z);

    std::vector<std::vector<double>> vertices((N_x + 1) * (N_y + 1) * (N_z + 1),
                                              std::vector<double>(3, 0.0));

    for (int kk = 0; kk <= N_z; ++kk)
    {
        for (int jj = 0; jj <= N_y; ++jj)
        {
            for (int ii = 0; ii <= N_x; ++ii)
            {
                const int node_counter =
                    ii + jj * (N_x + 1) + kk * (N_y + 1) * (N_x + 1);
                vertices[node_counter][0] = ii * delta_x;
                vertices[node_counter][1] = jj * delta_y;
                vertices[node_counter][2] = kk * delta_z;
            }
        }
    }

    return vertices;
}
