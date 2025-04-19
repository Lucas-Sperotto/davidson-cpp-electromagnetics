#include "globals.h"

// Índices
int NUM_NODES = 0;
int NUM_ELEMENTS = 0;
int NUM_EDGES = 0;
int NUM_DOFS = 0;

// Malha e topologia
std::vector<std::vector<int>> ELEMENTS;
std::vector<std::vector<int>> ELEMENT_EDGES;
std::vector<std::vector<int>> EDGES;
std::vector<std::vector<int>> EDGECONXELEM;
std::vector<std::vector<int>> EDGECONXELEMEDGES;
std::vector<std::vector<double>> NODE_COORD;
std::vector<std::vector<int>> LOCALEDGENODES = {
    {0, 1}, {0, 2}, {1, 2}
};

// Para funções como find_local_dofs
std::vector<int> dof_RWG;
std::vector<std::vector<int>> DOFLOCALNUM;
