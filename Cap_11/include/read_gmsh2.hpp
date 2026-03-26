#ifndef CAP_11_READ_GMSH2_HPP
#define CAP_11_READ_GMSH2_HPP

#include <string>
#include <vector>

struct GmshReadResult
{
    int num_nodes = 0;
    int line_counter = 0;
    int tri_counter = 0;
    int tet_counter = 0;
    std::vector<std::vector<double>> nodes;
    std::vector<std::vector<int>> lines;
    std::vector<std::vector<int>> triangles;
    std::vector<std::vector<int>> tets;
};

GmshReadResult read_gmsh2(const std::string &filename);

#endif
