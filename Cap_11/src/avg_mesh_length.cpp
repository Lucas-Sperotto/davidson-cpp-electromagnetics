#include "../include/avg_mesh_length.hpp"
#include "../include/globals3d.hpp"

#include <cmath>
#include <numeric>
#include <vector>

namespace
{
double edge_length(const std::vector<double> &a, const std::vector<double> &b)
{
    return std::sqrt(std::pow(b[0] - a[0], 2) +
                     std::pow(b[1] - a[1], 2) +
                     std::pow(b[2] - a[2], 2));
}
}

// Translation of avg_mesh_length.m:
// average mesh length of the tetrahedral mesh.
double avg_mesh_length(const std::vector<std::vector<double>> &vertices)
{
    std::vector<double> h_elem(NUM_ELEMENTS, 0.0);
    for (int ielem = 0; ielem < NUM_ELEMENTS; ++ielem)
    {
        std::vector<double> length(6, 0.0);
        for (int jedge = 0; jedge < 6; ++jedge)
        {
            const auto &vertex_1 =
                vertices[ELEMENTS[ielem][LOCALEDGENODES[jedge][0]]];
            const auto &vertex_2 =
                vertices[ELEMENTS[ielem][LOCALEDGENODES[jedge][1]]];
            length[jedge] = edge_length(vertex_1, vertex_2);
        }
        h_elem[ielem] =
            std::accumulate(length.begin(), length.end(), 0.0) / 6.0;
    }

    return std::accumulate(h_elem.begin(), h_elem.end(), 0.0) /
           static_cast<double>(h_elem.size());
}
