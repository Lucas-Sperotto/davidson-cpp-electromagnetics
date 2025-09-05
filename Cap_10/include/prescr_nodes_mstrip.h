#ifndef PRE_SCR_NODES_MSTRIP_H
#define PRE_SCR_NODES_MSTRIP_H

#include <vector>

// Marca os nós com condição de contorno inhomogênea (prescrita no condutor central)
void prescr_nodes_mstrip(double a, double b, double h, double w,
                         std::vector<int> &node_flag, int &num_pre_nodes);

#endif