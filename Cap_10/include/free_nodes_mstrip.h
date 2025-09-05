#ifndef FREE_NODES_MSTRIP_H
#define FREE_NODES_MSTRIP_H

#include <vector>

void free_nodes_mstrip(double a, double b, double h, double w,
                       std::vector<int> &node_flag, int &num_free_nodes);

#endif