#include "../include/internal_tet_mesh.hpp"

#include <vector>

namespace
{
int node_index(int ii, int jj, int kk, int N_x, int N_y)
{
    return ii + jj * (N_x + 1) + kk * (N_y + 1) * (N_x + 1);
}
}

// Helper for the internal mesh path used by Eigen3D_*:
// split each brick into six tetrahedra along one body diagonal.
std::vector<std::vector<int>> split_brick_mesh_to_tets(int N_x, int N_y, int N_z)
{
    std::vector<std::vector<int>> tets;
    tets.reserve(6 * N_x * N_y * N_z);

    for (int kk = 0; kk < N_z; ++kk)
    {
        for (int jj = 0; jj < N_y; ++jj)
        {
            for (int ii = 0; ii < N_x; ++ii)
            {
                const int n000 = node_index(ii, jj, kk, N_x, N_y);
                const int n100 = node_index(ii + 1, jj, kk, N_x, N_y);
                const int n010 = node_index(ii, jj + 1, kk, N_x, N_y);
                const int n110 = node_index(ii + 1, jj + 1, kk, N_x, N_y);
                const int n001 = node_index(ii, jj, kk + 1, N_x, N_y);
                const int n101 = node_index(ii + 1, jj, kk + 1, N_x, N_y);
                const int n011 = node_index(ii, jj + 1, kk + 1, N_x, N_y);
                const int n111 = node_index(ii + 1, jj + 1, kk + 1, N_x, N_y);

                tets.push_back({n000, n100, n110, n111});
                tets.push_back({n000, n110, n010, n111});
                tets.push_back({n000, n010, n011, n111});
                tets.push_back({n000, n011, n001, n111});
                tets.push_back({n000, n001, n101, n111});
                tets.push_back({n000, n101, n100, n111});
            }
        }
    }

    return tets;
}
