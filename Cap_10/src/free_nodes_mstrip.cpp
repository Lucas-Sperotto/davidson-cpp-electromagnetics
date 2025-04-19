#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>

extern int NUM_NODES;
extern std::vector<std::vector<double>> NODE_COORD;

void free_nodes_mstrip(double a, double b, double h, double w,
                       std::vector<int>& node_flag, int& num_free_nodes) {
    node_flag.resize(NUM_NODES, 1); // Inicializa todos como livres
    double eps = std::numeric_limits<double>::epsilon();

    for (int inode = 0; inode < NUM_NODES; ++inode) {
        double x = NODE_COORD[inode][0];
        double y = NODE_COORD[inode][1];

        // PEC: y = 0
        if (std::abs(y) < eps)
            node_flag[inode] = 0;

        // PEC: x = a/2
        if (std::abs(x - a / 2.0) < eps)
            node_flag[inode] = 0;

        // PEC: y = b
        if (std::abs(y - b) < eps)
            node_flag[inode] = 0;

        // Centro condutor (inhomogeneamente prescrito): y = h, |x| <= w/2
        if (std::abs(y - h) < eps && std::abs(x) <= w / 2.0)
            node_flag[inode] = 0;
    }

    num_free_nodes = std::count(node_flag.begin(), node_flag.end(), 1);
}
