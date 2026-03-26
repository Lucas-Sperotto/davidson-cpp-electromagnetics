#include "../include/facemake3D.hpp"
#include "../include/globals3d.hpp"

#include <map>
#include <tuple>

// Translation of facemake3D.m:
// assign unique global faces to each tetrahedron.
void facemake3D()
{
    FACES.clear();
    ELEMENT_FACES.assign(NUM_ELEMENTS, std::vector<int>(4, -1));

    std::map<std::tuple<int, int, int>, int> face_map;
    for (int ielem = 0; ielem < NUM_ELEMENTS; ++ielem)
    {
        for (int jface = 0; jface < 4; ++jface)
        {
            const int n1 = ELEMENTS[ielem][LOCALFACENODES[jface][0]];
            const int n2 = ELEMENTS[ielem][LOCALFACENODES[jface][1]];
            const int n3 = ELEMENTS[ielem][LOCALFACENODES[jface][2]];
            const auto key = std::make_tuple(n1, n2, n3);

            auto [it, inserted] =
                face_map.emplace(key, static_cast<int>(FACES.size()));
            if (inserted)
                FACES.push_back({n1, n2, n3});

            ELEMENT_FACES[ielem][jface] = it->second;
        }
    }

    NUM_FACES = static_cast<int>(FACES.size());
}
