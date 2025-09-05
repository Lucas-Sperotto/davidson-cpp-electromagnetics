#include <vector>
#include <cmath>
#include <cassert>

#include "globals.h"

// Gera malha triangular uniforme
void trimesh(double a, double b, int Nx, int Ny,
             std::vector<std::vector<double>> &x_nodes,
             std::vector<std::vector<double>> &y_nodes)
{
    NUM_ELEMENTS = 2 * Nx * Ny;
    ELEMENTS.resize(NUM_ELEMENTS, std::vector<int>(3));

    // Define os elementos
    for (int jj = 0; jj < Ny; ++jj)
    {
        for (int ii = 0; ii < Nx; ++ii)
        {
            int elem_offset = jj * 2 * Nx;
            int node_offset = jj * (Nx + 1);

            ELEMENTS[ii + elem_offset] = {ii + node_offset,
                                          ii + 1 + node_offset,
                                          ii + Nx + 1 + node_offset};

            ELEMENTS[ii + elem_offset + Nx] = {ii + 1 + node_offset,
                                               ii + Nx + 1 + node_offset,
                                               ii + Nx + 2 + node_offset};
        }
    }

    NUM_NODES = (Nx + 1) * (Ny + 1);
    NODE_COORD.resize(NUM_NODES, std::vector<double>(2));

    double del_x = a / (double)Nx;
    double del_y = b / (double)Ny;
    int counter = 0;

    for (int jj = 0; jj <= Ny; ++jj)
    {
        for (int ii = 0; ii <= Nx; ++ii)
        {
            NODE_COORD[counter][0] = ii * del_x;
            NODE_COORD[counter][1] = jj * del_y;
            ++counter;
        }
    }

    x_nodes = {{0.0, a / (double)Nx, a}};
    y_nodes = {{0.0, b / (double)Ny, b}};
}
