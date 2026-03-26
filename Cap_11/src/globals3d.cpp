#include "../include/globals3d.hpp"

int NUM_NODES = 0;
int NUM_ELEMENTS = 0;
int NUM_EDGES = 0;
int NUM_FACES = 0;
int NUM_DOFS = 0;

std::vector<std::vector<int>> ELEMENTS;
std::vector<std::vector<int>> ELEMENT_EDGES;
std::vector<std::vector<int>> ELEMENT_FACES;
std::vector<std::vector<int>> EDGES;
std::vector<std::vector<int>> FACES;
std::vector<std::vector<double>> NODE_COORD;
std::vector<std::vector<int>> LOCALEDGENODES;
std::vector<std::vector<int>> LOCALFACENODES;

void reset_global_mesh3d()
{
    NUM_NODES = 0;
    NUM_ELEMENTS = 0;
    NUM_EDGES = 0;
    NUM_FACES = 0;
    NUM_DOFS = 0;

    ELEMENTS.clear();
    ELEMENT_EDGES.clear();
    ELEMENT_FACES.clear();
    EDGES.clear();
    FACES.clear();
    NODE_COORD.clear();
}

void initialize_local_topology3d()
{
    LOCALEDGENODES = {
        {0, 1},
        {0, 2},
        {0, 3},
        {1, 2},
        {1, 3},
        {2, 3}};

    LOCALFACENODES = {
        {0, 1, 2},
        {0, 1, 3},
        {0, 2, 3},
        {1, 2, 3}};
}
