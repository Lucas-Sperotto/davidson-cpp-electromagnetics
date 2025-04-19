#ifndef GLOBALS_H
#define GLOBALS_H

#include <vector>

// Índices
extern int NUM_NODES;
extern int NUM_ELEMENTS;
extern int NUM_EDGES;
extern int NUM_DOFS;

// Malha e topologia
extern std::vector<std::vector<int>> ELEMENTS;
extern std::vector<std::vector<int>> ELEMENT_EDGES;
extern std::vector<std::vector<int>> EDGES;
extern std::vector<std::vector<int>> EDGECONXELEM;
extern std::vector<std::vector<int>> EDGECONXELEMEDGES;
extern std::vector<std::vector<double>> NODE_COORD;
extern std::vector<std::vector<int>> LOCALEDGENODES;

// Para funções como find_local_dofs
extern std::vector<int> dof_RWG;
extern std::vector<std::vector<int>> DOFLOCALNUM;

#endif // GLOBALS_H
