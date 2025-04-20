#include <vector>

#include "globals.h"

// Renumera os graus de liberdade LTQN: dois por aresta (e1, e2) e dois por elemento (f1, f2)
void renumber_dof_LTQN(const std::vector<int>& dof_free_flag,
                       std::vector<int>& dof_e1,
                       std::vector<int>& dof_e2,
                       std::vector<int>& dof_f1,
                       std::vector<int>& dof_f2) {
    dof_e1.assign(NUM_EDGES, 0); // 0 indica grau prescrito
    dof_e2.assign(NUM_EDGES, 0);
    dof_f1.resize(NUM_ELEMENTS);
    dof_f2.resize(NUM_ELEMENTS);

    int counter = 0;

    // Primeira família (e1)
    for (int i_edge = 0; i_edge < NUM_EDGES; ++i_edge) {
        if (dof_free_flag[i_edge]) {
            ++counter;
            dof_e1[i_edge] = counter;
        }
    }

    // Segunda família (e2)
    for (int i_edge = 0; i_edge < NUM_EDGES; ++i_edge) {
        if (dof_free_flag[i_edge]) {
            ++counter;
            dof_e2[i_edge] = counter;
        }
    }

    // DOFs de face
    for (int i_elem = 0; i_elem < NUM_ELEMENTS; ++i_elem) {
        ++counter;
        dof_f1[i_elem] = counter;
    }
    for (int i_elem = 0; i_elem < NUM_ELEMENTS; ++i_elem) {
        ++counter;
        dof_f2[i_elem] = counter;
    }

    NUM_DOFS = counter;
}
