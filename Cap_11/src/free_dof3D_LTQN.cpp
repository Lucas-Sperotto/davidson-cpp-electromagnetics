#include "../include/free_dof3D_LTQN.hpp"
#include "../include/globals3d.hpp"

#include <cmath>
#include <utility>
#include <vector>

// Translation of free_dof3D_LTQN.m:
// flag edge and face DOFs as free or prescribed.
std::pair<std::vector<int>, std::vector<int>> free_dof3D_LTQN(
    double a,
    double b,
    double c,
    const std::vector<std::vector<double>> &node_coord,
    bool symmetry_flag,
    double toler)
{
    std::vector<int> edge_free_flag(NUM_EDGES, 1);
    std::vector<int> face_free_flag(NUM_FACES, 1);

    for (int i_edge = 0; i_edge < NUM_EDGES; ++i_edge)
    {
        const int node1 = EDGES[i_edge][0];
        const int node2 = EDGES[i_edge][1];

        if (std::abs(node_coord[node1][0]) < toler &&
            std::abs(node_coord[node2][0]) < toler)
            edge_free_flag[i_edge] = 0;
        if (!symmetry_flag &&
            std::abs(node_coord[node1][0] - a) < toler &&
            std::abs(node_coord[node2][0] - a) < toler)
            edge_free_flag[i_edge] = 0;
        if (std::abs(node_coord[node1][1]) < toler &&
            std::abs(node_coord[node2][1]) < toler)
            edge_free_flag[i_edge] = 0;
        if (!symmetry_flag &&
            std::abs(node_coord[node1][1] - b) < toler &&
            std::abs(node_coord[node2][1] - b) < toler)
            edge_free_flag[i_edge] = 0;
        if (std::abs(node_coord[node1][2]) < toler &&
            std::abs(node_coord[node2][2]) < toler)
            edge_free_flag[i_edge] = 0;
        if (std::abs(node_coord[node1][2] - c) < toler &&
            std::abs(node_coord[node2][2] - c) < toler)
            edge_free_flag[i_edge] = 0;
    }

    for (int i_face = 0; i_face < NUM_FACES; ++i_face)
    {
        const int node1 = FACES[i_face][0];
        const int node2 = FACES[i_face][1];
        const int node3 = FACES[i_face][2];

        if (std::abs(node_coord[node1][0]) < toler &&
            std::abs(node_coord[node2][0]) < toler &&
            std::abs(node_coord[node3][0]) < toler)
            face_free_flag[i_face] = 0;
        if (!symmetry_flag &&
            std::abs(node_coord[node1][0] - a) < toler &&
            std::abs(node_coord[node2][0] - a) < toler &&
            std::abs(node_coord[node3][0] - a) < toler)
            face_free_flag[i_face] = 0;
        if (std::abs(node_coord[node1][1]) < toler &&
            std::abs(node_coord[node2][1]) < toler &&
            std::abs(node_coord[node3][1]) < toler)
            face_free_flag[i_face] = 0;
        if (!symmetry_flag &&
            std::abs(node_coord[node1][1] - b) < toler &&
            std::abs(node_coord[node2][1] - b) < toler &&
            std::abs(node_coord[node3][1] - b) < toler)
            face_free_flag[i_face] = 0;
        if (std::abs(node_coord[node1][2]) < toler &&
            std::abs(node_coord[node2][2]) < toler &&
            std::abs(node_coord[node3][2]) < toler)
            face_free_flag[i_face] = 0;
        if (std::abs(node_coord[node1][2] - c) < toler &&
            std::abs(node_coord[node2][2] - c) < toler &&
            std::abs(node_coord[node3][2] - c) < toler)
            face_free_flag[i_face] = 0;
    }

    return {edge_free_flag, face_free_flag};
}
