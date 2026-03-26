#include "trimesh3D.hpp"

#include "mom3d_globals.hpp"

namespace cap06 {

void trimesh3D(const double a, const double b, const int Nx, const int Ny)
{
    NUM_ELEMENTS = 2 * Nx * Ny;
    ELEMENTS.resize(NUM_ELEMENTS, 3);

    for (int jj = 0; jj < Ny; ++jj) {
        for (int ii = 0; ii < Nx; ++ii) {
            const int elem_offset = jj * 2 * Nx;
            const int node_offset = jj * (Nx + 1);

            ELEMENTS.row(ii + elem_offset) << ii + node_offset,
                                              ii + 1 + node_offset,
                                              ii + Nx + 1 + node_offset;
            ELEMENTS.row(ii + elem_offset + Nx) << ii + 1 + node_offset,
                                                   ii + Nx + 1 + node_offset,
                                                   ii + Nx + 2 + node_offset;
        }
    }

    NUM_NODES = (Nx + 1) * (Ny + 1);
    NODE_COORD.resize(NUM_NODES, 3);

    const double del_x = a / static_cast<double>(Nx);
    const double del_y = b / static_cast<double>(Ny);

    int counter = 0;
    for (int jj = 0; jj <= Ny; ++jj) {
        for (int ii = 0; ii <= Nx; ++ii) {
            NODE_COORD.row(counter) << ii * del_x, jj * del_y, 0.0;
            ++counter;
        }
    }

    EDGES.resize(0, 2);
    ELEMENT_EDGES.resize(0, 3);
    ELL.resize(0);
    NUM_EDGES = 0;
}

}  // namespace cap06
