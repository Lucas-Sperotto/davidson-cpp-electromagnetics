#include <vector>

#include "globals.h"

void find_local_dofs() {
    DOFLOCALNUM.resize(dof_RWG.size(), std::vector<int>(2, -1)); // Inicializa com -1

    for (int iedge = 0; iedge < NUM_EDGES; ++iedge) {
        bool edge_found = false;
        bool terminate = false;

        for (int jelem = 0; jelem < NUM_ELEMENTS && !terminate; ++jelem) {
            for (int kedge = 0; kedge < 3; ++kedge) {
                if (ELEMENT_EDGES[jelem][kedge] == iedge) {
                    if (!edge_found) {
                        DOFLOCALNUM[dof_RWG[iedge]][0] = kedge;
                        edge_found = true;
                    } else {
                        DOFLOCALNUM[dof_RWG[iedge]][1] = kedge;
                        terminate = true;
                        break;
                    }
                }
            }
        }
    }
}
