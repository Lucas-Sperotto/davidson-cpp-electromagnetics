#include "edgemake_MoM.hpp"

#include "mom3d_globals.hpp"

namespace cap06 {

void edgemake_MoM()
{
    EDGES = Eigen::MatrixXi::Zero(3 * NUM_ELEMENTS, 2);
    ELEMENT_EDGES = Eigen::MatrixXi::Zero(NUM_ELEMENTS, 3);

    EDGES.row(0) << ELEMENTS(0, LOCALEDGENODES(0, 0)), ELEMENTS(0, LOCALEDGENODES(0, 1));
    EDGES.row(1) << ELEMENTS(0, LOCALEDGENODES(1, 0)), ELEMENTS(0, LOCALEDGENODES(1, 1));
    EDGES.row(2) << ELEMENTS(0, LOCALEDGENODES(2, 0)), ELEMENTS(0, LOCALEDGENODES(2, 1));
    ELEMENT_EDGES.row(0) << 0, 1, 2;

    int edge_counter = 3;
    for (int ielem = 1; ielem < NUM_ELEMENTS; ++ielem) {
        for (int jedge = 0; jedge < 3; ++jedge) {
            const int node_a = ELEMENTS(ielem, LOCALEDGENODES(jedge, 0));
            const int node_b = ELEMENTS(ielem, LOCALEDGENODES(jedge, 1));

            bool new_edge = true;
            for (int kedge = 0; kedge < edge_counter; ++kedge) {
                if (EDGES(kedge, 0) == node_a && EDGES(kedge, 1) == node_b) {
                    new_edge = false;
                    ELEMENT_EDGES(ielem, jedge) = kedge;
                    break;
                }
            }

            if (new_edge) {
                EDGES.row(edge_counter) << node_a, node_b;
                ELEMENT_EDGES(ielem, jedge) = edge_counter;
                ++edge_counter;
            }
        }
    }

    NUM_EDGES = edge_counter;
    EDGES.conservativeResize(NUM_EDGES, 2);
    ELL.resize(NUM_EDGES);

    for (int iedge = 0; iedge < NUM_EDGES; ++iedge) {
        ELL(iedge) = (NODE_COORD.row(EDGES(iedge, 1)) - NODE_COORD.row(EDGES(iedge, 0))).norm();
    }
}

}  // namespace cap06
