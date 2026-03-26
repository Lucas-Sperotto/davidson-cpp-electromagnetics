#include "../include/edgemake3D.hpp"
#include "../include/globals3d.hpp"

#include <map>
#include <utility>

// Translation of edgemake3D.m:
// assign unique global edges to each tetrahedron.
void edgemake3D()
{
    EDGES.clear();
    ELEMENT_EDGES.assign(NUM_ELEMENTS, std::vector<int>(6, -1));

    std::map<std::pair<int, int>, int> edge_map;
    for (int ielem = 0; ielem < NUM_ELEMENTS; ++ielem)
    {
        for (int jedge = 0; jedge < 6; ++jedge)
        {
            const int n1 = ELEMENTS[ielem][LOCALEDGENODES[jedge][0]];
            const int n2 = ELEMENTS[ielem][LOCALEDGENODES[jedge][1]];
            const std::pair<int, int> key = {n1, n2};

            auto [it, inserted] =
                edge_map.emplace(key, static_cast<int>(EDGES.size()));
            if (inserted)
                EDGES.push_back({n1, n2});

            ELEMENT_EDGES[ielem][jedge] = it->second;
        }
    }

    NUM_EDGES = static_cast<int>(EDGES.size());
}
