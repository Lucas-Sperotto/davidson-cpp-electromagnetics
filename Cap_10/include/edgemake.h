#ifndef EDGEMAKE_H
#define EDGEMAKE_H

#include <vector>

// Declaração da função
void edgemake();

// Declarações externas (devem estar definidas em outro lugar no código)
extern int NUM_ELEMENTS, NUM_EDGES;
extern std::vector<std::vector<int>> ELEMENTS;
extern std::vector<std::vector<int>> ELEMENT_EDGES;
extern std::vector<std::vector<int>> EDGES;
extern std::vector<std::vector<int>> LOCALEDGENODES;

#endif // EDGEMAKE_H
