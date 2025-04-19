#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <tuple>

#include "trimesh.h"
#include "edgemake.h"
#include "free_dof.h"
#include "renumber_dof_LTQN.h"
#include "sandt_LTQN.h"
#include "free_nodes.h"
#include "TEeig_err.h"

using namespace Eigen;

#include "globals.h"

const double eps_0 = 8.854e-12;
const double mu_0 = 4 * M_PI * 1e-7;
const double eps_r = 1.0;
const double mu_r = 1.0;

void Eigen2D_LTQN()
{
    // Dimensões do guia
    double a = 2.286e-2; // X-band guide
    double b = 1.016e-2;

    // Gerar a malha triangular
    std::vector<std::vector<double>> x_nodes, y_nodes;
    trimesh(a, b, 4, 2, x_nodes, y_nodes); // função a ser implementada
    // A função também deve preencher ELEMENTS e NODE_COORD

    // (Opcional) Imprimir coordenadas e IDs de nós/elementos — substitui triplot + text do MATLAB
    for (int i = 0; i < NUM_NODES; ++i)
    {
        std::cout << "Node " << i << ": (" << NODE_COORD[i][0] << ", " << NODE_COORD[i][1] << ")\n";
    }
    for (int e = 0; e < NUM_ELEMENTS; ++e)
    {
        double xc = 0.0, yc = 0.0;
        for (int j = 0; j < 3; ++j)
        {
            xc += NODE_COORD[ELEMENTS[e][j]][0];
            yc += NODE_COORD[ELEMENTS[e][j]][1];
        }
        xc /= 3.0;
        yc /= 3.0;
        std::cout << "Element " << e << " center: (" << xc << ", " << yc << ")\n";
    }

    // Pausa manual: aguarda o usuário
    std::cout << "Press ENTER to continue...\n";
    std::cin.get();

    // Criar as arestas e associar com os elementos
    edgemake(); // Função para preencher EDGES, ELEMENT_EDGES e NUM_EDGES

    // Determina os graus de liberdade livres nas arestas (PEC boundaries)
    std::vector<int> edge_free_flag = free_dof(a, b);

    // Renumera os graus de liberdade: CTLN, LTLN, faces LTQN
    std::vector<int> dof_e1, dof_e2;
    std::vector<int> dof_f1(NUM_ELEMENTS), dof_f2(NUM_ELEMENTS);
    renumber_dof_LTQN(edge_free_flag, dof_e1, dof_e2, dof_f1, dof_f2); // Atualiza NUM_DOFS também

    // Inicializa matrizes globais S e T
    std::vector<std::vector<double>> S(NUM_DOFS, std::vector<double>(NUM_DOFS, 0.0));
    std::vector<std::vector<double>> T(NUM_DOFS, std::vector<double>(NUM_DOFS, 0.0));

    // Loop sobre todos os elementos para montar as matrizes globais
    for (int ielem = 0; ielem < NUM_ELEMENTS; ++ielem)
    {
        std::vector<int> trinodes = ELEMENTS[ielem];

        double x1 = NODE_COORD[trinodes[0]][0], y1 = NODE_COORD[trinodes[0]][1];
        double x2 = NODE_COORD[trinodes[1]][0], y2 = NODE_COORD[trinodes[1]][1];
        double x3 = NODE_COORD[trinodes[2]][0], y3 = NODE_COORD[trinodes[2]][1];

        auto [S_elem, T_elem] = sandt_LTQN(x1, y1, x2, y2, x3, y3); // retorna matrizes 8x8

        int ll_f1 = dof_f1[ielem];
        int ll_f2 = dof_f2[ielem];

        for (int jedge = 0; jedge < 3; ++jedge)
        {
            int jj_e1 = dof_e1[ELEMENT_EDGES[ielem][jedge]];
            int jj_e2 = dof_e2[ELEMENT_EDGES[ielem][jedge]];

            for (int kedge = 0; kedge < 3; ++kedge)
            {
                int kk_e1 = dof_e1[ELEMENT_EDGES[ielem][kedge]];
                int kk_e2 = dof_e2[ELEMENT_EDGES[ielem][kedge]];

                if (jj_e1 && kk_e1)
                {
                    S[jj_e1][kk_e1] += S_elem[jedge][kedge];
                    S[jj_e1][kk_e2] += S_elem[jedge][kedge + 3];
                    S[jj_e2][kk_e2] += S_elem[jedge + 3][kedge + 3];

                    T[jj_e1][kk_e1] += T_elem[jedge][kedge];
                    T[jj_e1][kk_e2] += T_elem[jedge][kedge + 3];
                    T[jj_e2][kk_e2] += T_elem[jedge + 3][kedge + 3];
                }
            }

            if (jj_e1)
            {
                S[jj_e1][ll_f1] += S_elem[jedge][6];
                S[jj_e1][ll_f2] += S_elem[jedge][7];
                S[jj_e2][ll_f1] += S_elem[jedge + 3][6];
                S[jj_e2][ll_f2] += S_elem[jedge + 3][7];

                T[jj_e1][ll_f1] += T_elem[jedge][6];
                T[jj_e1][ll_f2] += T_elem[jedge][7];
                T[jj_e2][ll_f1] += T_elem[jedge + 3][6];
                T[jj_e2][ll_f2] += T_elem[jedge + 3][7];
            }
        }

        S[ll_f1][ll_f1] += S_elem[6][6];
        S[ll_f1][ll_f2] += S_elem[6][7];
        S[ll_f2][ll_f2] += S_elem[7][7];

        T[ll_f1][ll_f1] += T_elem[6][6];
        T[ll_f1][ll_f2] += T_elem[6][7];
        T[ll_f2][ll_f2] += T_elem[7][7];
    }

    // Preencher as partes inferiores simétricas
    for (int j = 1; j < NUM_DOFS; ++j)
    {
        for (int k = 0; k < j; ++k)
        {
            S[j][k] = S[k][j];
            T[j][k] = T[k][j];
        }
    }

    // Converter std::vector para Eigen::MatrixXd
    MatrixXd S_mat(NUM_DOFS, NUM_DOFS), T_mat(NUM_DOFS, NUM_DOFS);
    for (int i = 0; i < NUM_DOFS; ++i)
    {
        for (int j = 0; j < NUM_DOFS; ++j)
        {
            S_mat(i, j) = S[i][j] / mu_r;
            T_mat(i, j) = T[i][j] * eps_r;
        }
    }

    // Resolver o problema generalizado de autovalores
    GeneralizedSelfAdjointEigenSolver<MatrixXd> solver(S_mat, T_mat);
    VectorXd eigvals = solver.eigenvalues();
    MatrixXd eigvecs = solver.eigenvectors();

    // Extração dos autovalores positivos e ordenação
    std::vector<std::pair<double, int>> indexed_kc;
    for (int i = 0; i < eigvals.size(); ++i)
    {
        double val = eigvals[i];
        if (val > 0)
            indexed_kc.emplace_back(std::sqrt(val), i);
    }
    std::sort(indexed_kc.begin(), indexed_kc.end());

    // Avaliação de erro relativa com base nos autovalores esperados
    std::vector<int> node_flag;
    int num_free_nodes = 0;
    free_nodes(a, b, node_flag, num_free_nodes);
    int num_free_edges = std::count_if(edge_free_flag.begin(), edge_free_flag.end(), [](int f)
                                       { return f == 1; });
    int num_zero_eigvals = num_free_nodes + num_free_edges;

    std::vector<double> TEeigvalues;
    for (const auto& [kc, _] : indexed_kc) {
        TEeigvalues.push_back(kc);
    }
    // Avaliação de erro dos primeiros modos úteis
    std::vector<double> rel_err = TEeig_err(a, b, TEeigvalues, std::min(NUM_DOFS - num_zero_eigvals, 8), num_zero_eigvals);
    


    // Opcional: salvar dados
    std::string filename = "eigdata_LTQN_" + std::to_string(NUM_ELEMENTS) + ".txt";
    std::ofstream fout(filename);
    for (const auto &[val, idx] : indexed_kc)
    {
        fout << val << "\n";
    }
    fout.close();
}

int main()
{
    Eigen2D_LTQN();
    return 0;
}
