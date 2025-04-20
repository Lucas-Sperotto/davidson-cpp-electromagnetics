#include <vector>
#include <cmath>
#include <limits>
#include "globals.h"

std::vector<int> free_dof(double a, double b)
{
    std::vector<int> dof_flag(NUM_EDGES, 1); // Inicializa todos como livres
    double eps = std::numeric_limits<double>::epsilon();

    for (int i_edge = 0; i_edge < NUM_EDGES; ++i_edge)
    {
        int node1 = EDGES[i_edge][0];
        int node2 = EDGES[i_edge][1];

        double x1 = NODE_COORD[node1][0], y1 = NODE_COORD[node1][1];
        double x2 = NODE_COORD[node2][0], y2 = NODE_COORD[node2][1];

        // Condições de contorno PEC (x=0, x=a, y=0, y=b)
        bool on_y0 = std::abs(y1) < eps && std::abs(y2) < eps;
        bool on_yb = std::abs(y1 - b) < eps && std::abs(y2 - b) < eps;
        bool on_x0 = std::abs(x1) < eps && std::abs(x2) < eps;
        bool on_xa = std::abs(x1 - a) < eps && std::abs(x2 - a) < eps;

        if (on_y0 || on_yb || on_x0 || on_xa)
        {
            dof_flag[i_edge] = 0; // Prescrito
        }
    }

    return dof_flag;
}
