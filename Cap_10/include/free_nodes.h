#ifndef FREE_NODES_H
#define FREE_NODES_H

#include <vector>

void free_nodes(double a, double b, std::vector<int>& node_flag, int& num_free_nodes);
void free_nodes_mstrip(double a, double b, double h, double w, std::vector<int>& node_flag, int& num_free_nodes);

#endif