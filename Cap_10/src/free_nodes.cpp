#include <vector>
#include <cmath>
#include <algorithm>

extern int NUM_NODES;
extern std::vector<std::vector<double>> NODE_COORD;

void free_nodes(double a, double b, std::vector<int>& node_flag, int& num_free_nodes) {
    node_flag.resize(NUM_NODES, 1);  // Começa marcando todos como livres

    for (int inode = 0; inode < NUM_NODES; ++inode) {
        double x = NODE_COORD[inode][0];
        double y = NODE_COORD[inode][1];
        double eps = std::numeric_limits<double>::epsilon();

        if (std::abs(x) < eps || std::abs(x - a) < eps ||
            std::abs(y) < eps || std::abs(y - b) < eps) {
            node_flag[inode] = 0;
        }
    }

    num_free_nodes = std::count(node_flag.begin(), node_flag.end(), 1);
}
