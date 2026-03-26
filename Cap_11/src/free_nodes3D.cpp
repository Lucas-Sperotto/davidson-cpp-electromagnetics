#include "../include/free_nodes3D.hpp"
#include "../include/globals3d.hpp"

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

// Translation of free_nodes3D.m:
// flag nodes as free or prescribed for the PEC cavity.
std::pair<std::vector<int>, int> free_nodes3D(double a,
                                              double b,
                                              double c,
                                              const std::vector<std::vector<double>> &node_coord,
                                              bool symmetry_flag,
                                              double toler)
{
    std::vector<int> node_flag(NUM_NODES, 1);

    for (int inode = 0; inode < NUM_NODES; ++inode)
    {
        if (std::abs(node_coord[inode][0]) < toler)
            node_flag[inode] = 0;
        if (!symmetry_flag && std::abs(node_coord[inode][0] - a) < toler)
            node_flag[inode] = 0;
        if (std::abs(node_coord[inode][1]) < toler)
            node_flag[inode] = 0;
        if (!symmetry_flag && std::abs(node_coord[inode][1] - b) < toler)
            node_flag[inode] = 0;
        if (std::abs(node_coord[inode][2]) < toler)
            node_flag[inode] = 0;
        if (std::abs(node_coord[inode][2] - c) < toler)
            node_flag[inode] = 0;
    }

    const int num_free_nodes =
        static_cast<int>(std::count(node_flag.begin(), node_flag.end(), 1));
    return {node_flag, num_free_nodes};
}
