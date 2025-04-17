#include <vector>
#include <algorithm>

extern int NUM_ELEMENTS, NUM_EDGES;
extern std::vector<std::vector<int>> ELEMENTS;         // ELEMENTS[e][0..2] são os vértices
extern std::vector<std::vector<int>> ELEMENT_EDGES;    // ELEMENT_EDGES[e][0..2]
extern std::vector<std::vector<int>> EDGES;            // EDGES[i] = {n1, n2}
extern std::vector<std::vector<int>> LOCALEDGENODES;   // LOCALEDGENODES[3] = {{0,1}, {0,2}, {1,2}}

void edgemake() {
    EDGES.clear();
    ELEMENT_EDGES.resize(NUM_ELEMENTS, std::vector<int>(3, -1));
    EDGES.push_back({ELEMENTS[0][0], ELEMENTS[0][1]});
    EDGES.push_back({ELEMENTS[0][0], ELEMENTS[0][2]});
    EDGES.push_back({ELEMENTS[0][1], ELEMENTS[0][2]});
    ELEMENT_EDGES[0] = {0, 1, 2};

    int edge_counter = 3;
    for (int ielem = 1; ielem < NUM_ELEMENTS; ++ielem) {
        for (int jedge = 0; jedge < 3; ++jedge) {
            int ni = ELEMENTS[ielem][LOCALEDGENODES[jedge][0]];
            int nj = ELEMENTS[ielem][LOCALEDGENODES[jedge][1]];

            // Ordenar os índices para garantir consistência
            int n1 = std::min(ni, nj);
            int n2 = std::max(ni, nj);
            std::vector<int> current_edge = {n1, n2};

            bool found = false;
            for (int k = 0; k < edge_counter; ++k) {
                if (EDGES[k] == current_edge) {
                    ELEMENT_EDGES[ielem][jedge] = k;
                    found = true;
                    break;
                }
            }

            if (!found) {
                EDGES.push_back(current_edge);
                ELEMENT_EDGES[ielem][jedge] = edge_counter;
                ++edge_counter;
            }
        }
    }

    NUM_EDGES = edge_counter;
}
