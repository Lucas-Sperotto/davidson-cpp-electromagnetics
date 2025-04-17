#include <vector>

extern int NUM_ELEMENTS;
extern int NUM_EDGES;
extern std::vector<std::vector<int>> ELEMENT_EDGES;  // ELEMENT_EDGES[elem][local_edge]
extern std::vector<int> dof_RWG;                     // Mapeamento edge -> dof index
extern std::vector<std::vector<int>> DOFLOCALNUM;    // DOFLOCALNUM[dof][2]

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
