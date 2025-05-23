#include <vector>
#include "globals.h"

// Renumera graus de liberdade nas arestas (elementos de Whitney)
// Retorna vetor dof_e1 onde dof_e1[i] = índice global se livre, 0 se prescrito
std::vector<int> renumber_dof(const std::vector<int> &dof_free_flag)
{
    std::vector<int> dof_e1(NUM_EDGES, 0); // Inicializa todos como prescritos (0)
    int counter = 0;

    for (int i_edge = 0; i_edge < NUM_EDGES; ++i_edge)
    {
        if (dof_free_flag[i_edge])
        {
            dof_e1[i_edge] = counter;
            ++counter;
        }
    }

    NUM_DOFS = counter; // Atualiza número total de DOFs livres
    return dof_e1;
}
