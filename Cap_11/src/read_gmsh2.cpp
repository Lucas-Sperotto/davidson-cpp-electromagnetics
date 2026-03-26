#include "../include/read_gmsh2.hpp"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
std::string trim_cr(std::string line)
{
    if (!line.empty() && line.back() == '\r')
        line.pop_back();
    return line;
}

void expect_line(std::ifstream &stream, const std::string &expected)
{
    std::string line;
    if (!std::getline(stream, line) || trim_cr(line) != expected)
        throw std::runtime_error("Formato inesperado no arquivo Gmsh.");
}
}

// Translation of read_gmsh2.m:
// read a Gmsh v2 ASCII mesh with lines, triangles and tetrahedra.
GmshReadResult read_gmsh2(const std::string &filename)
{
    std::ifstream input(filename);
    if (!input)
        throw std::runtime_error("Nao foi possivel abrir o arquivo de malha: " + filename);

    GmshReadResult result;
    std::string line;

    expect_line(input, "$MeshFormat");
    std::getline(input, line); // e.g. 2 0 8
    line = trim_cr(line);
    expect_line(input, "$EndMeshFormat");
    expect_line(input, "$Nodes");

    std::getline(input, line);
    line = trim_cr(line);
    result.num_nodes = std::stoi(line);
    result.nodes.assign(result.num_nodes, std::vector<double>(3, 0.0));
    for (int inode = 0; inode < result.num_nodes; ++inode)
    {
        std::getline(input, line);
        line = trim_cr(line);
        std::istringstream iss(line);
        int node_id = 0;
        double x = 0.0, y = 0.0, z = 0.0;
        iss >> node_id >> x >> y >> z;
        result.nodes[node_id - 1] = {x, y, z};
    }

    expect_line(input, "$EndNodes");
    expect_line(input, "$Elements");
    std::getline(input, line);
    line = trim_cr(line);
    const int num_elements = std::stoi(line);

    for (int element_index = 0; element_index < num_elements; ++element_index)
    {
        std::getline(input, line);
        line = trim_cr(line);
        std::istringstream iss(line);
        int elem_num = 0;
        int elem_type = 0;
        int num_tags = 0;
        iss >> elem_num >> elem_type >> num_tags;
        for (int tag_index = 0; tag_index < num_tags; ++tag_index)
        {
            int dummy = 0;
            iss >> dummy;
        }

        if (elem_type == 1)
        {
            int n1 = 0, n2 = 0;
            iss >> n1 >> n2;
            result.lines.push_back({n1 - 1, n2 - 1});
        }
        else if (elem_type == 2)
        {
            int n1 = 0, n2 = 0, n3 = 0;
            iss >> n1 >> n2 >> n3;
            result.triangles.push_back({n1 - 1, n2 - 1, n3 - 1});
        }
        else if (elem_type == 4)
        {
            int n1 = 0, n2 = 0, n3 = 0, n4 = 0;
            iss >> n1 >> n2 >> n3 >> n4;
            std::vector<int> tet = {n1 - 1, n2 - 1, n3 - 1, n4 - 1};
            std::sort(tet.begin(), tet.end());
            result.tets.push_back(tet);
        }
        else
        {
            // The MATLAB reader ignores unsupported element types too.
        }
    }

    expect_line(input, "$EndElements");
    result.line_counter = static_cast<int>(result.lines.size());
    result.tri_counter = static_cast<int>(result.triangles.size());
    result.tet_counter = static_cast<int>(result.tets.size());
    return result;
}
