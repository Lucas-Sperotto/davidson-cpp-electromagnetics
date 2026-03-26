#ifndef CAP_11_GLOBALS3D_HPP
#define CAP_11_GLOBALS3D_HPP

#include <vector>

extern int NUM_NODES;
extern int NUM_ELEMENTS;
extern int NUM_EDGES;
extern int NUM_FACES;
extern int NUM_DOFS;

extern std::vector<std::vector<int>> ELEMENTS;
extern std::vector<std::vector<int>> ELEMENT_EDGES;
extern std::vector<std::vector<int>> ELEMENT_FACES;
extern std::vector<std::vector<int>> EDGES;
extern std::vector<std::vector<int>> FACES;
extern std::vector<std::vector<double>> NODE_COORD;
extern std::vector<std::vector<int>> LOCALEDGENODES;
extern std::vector<std::vector<int>> LOCALFACENODES;

void reset_global_mesh3d();
void initialize_local_topology3d();

#endif
