#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>

#include "globals.h"

// Marca os nós com condição de contorno inhomogênea (prescrita no condutor central)
void prescr_nodes_mstrip(double a, double b, double h, double w,
                         std::vector<int>& node_flag, int& num_pre_nodes) {
    node_flag.assign(NUM_NODES, 0); // Inicialmente todos não prescritos
    double eps = std::numeric_limits<double>::epsilon();

    for (int inode = 0; inode < NUM_NODES; ++inode) {
        double x = NODE_COORD[inode][0];
        double y = NODE_COORD[inode][1];

        if (std::abs(y - h) < eps && std::abs(x) <= (w / 2.0))
            node_flag[inode] = 1;
    }

    num_pre_nodes = std::count(node_flag.begin(), node_flag.end(), 1);
}
