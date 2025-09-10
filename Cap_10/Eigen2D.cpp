#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <fstream>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <string>

#include "trimesh.h"
#include "edgemake.h"
#include "free_dof.h"
#include "renumber_dof.h"
#include "sandt.h"
#include "free_nodes.h"
#include "plot_field.h"
#include "simplex2D.h"
#include "whitney.h"

using namespace Eigen;

#include "globals.h"

const double eps_0 = 8.854e-12;
const double mu_0 = 4.0 * M_PI * 1e-7;
const double eps_r = 1.0;
const double mu_r = 1.0;

void write_csv(const std::string &path, const Eigen::MatrixXd &M)
{
    std::ofstream f(path);
    f.setf(std::ios::scientific);
    f << std::setprecision(17);
    for (int i = 0; i < M.rows(); ++i)
    {
        for (int j = 0; j < M.cols(); ++j)
        {
            if (j)
                f << ',';
            f << M(i, j);
        }
        f << '\n';
    }
}

void write_csv(const std::string &path, const Eigen::VectorXd &V)
{
    std::ofstream f(path);
    f.setf(std::ios::scientific);
    f << std::setprecision(17);
    for (int i = 0; i < V.size(); ++i)
    {
            f << V(i) << '\n';
    }
}

void Eigen2D()
{
    double a = 2.286e-2; // X band guide in cm
    double b = 1.016e-2;

    std::vector<std::vector<double>> x_nodes, y_nodes;
    trimesh(a, b, 8, 4, x_nodes, y_nodes); // também preenche ELEMENTS e NODE_COORD

    // for (int i = 0; i < NUM_NODES; ++i)
    //     std::cout << "Node " << i << ": (" << NODE_COORD[i][0] << ", " << NODE_COORD[i][1] << ")\n";

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
        // std::cout << "Element " << e << " center: (" << xc << ", " << yc << ")\n";
    }

    // std::cout << "Press ENTER to continue...\n";
    // std::cin.get();

    // Geração das arestas a partir dos elementos
    edgemake(); // preenche EDGES, ELEMENT_EDGES, atualiza NUM_EDGES

    // Determina quais arestas estão em contorno PEC e devem ser prescritas
    std::vector<int> dof_e1_free_flag = free_dof(a, b);

    // Renumera apenas os DOFs livres
    std::vector<int> dof_e1 = renumber_dof(dof_e1_free_flag); // também atualiza NUM_DOFS

    // for (int j = 0; j < dof_e1.size(); ++j)
    //     std::cout << "dof_e1[" << j << "] = " << dof_e1[j] << "\n";

    // Inicializa matrizes globais
    std::vector<std::vector<double>> S(NUM_DOFS, std::vector<double>(NUM_DOFS, 0.0));
    std::vector<std::vector<double>> T(NUM_DOFS, std::vector<double>(NUM_DOFS, 0.0));

    for (int ielem = 0; ielem < NUM_ELEMENTS; ++ielem)
    {
        const auto &trinodes = ELEMENTS[ielem];
        // std::cout << "Tri[" << ielem << "] = (" << trinodes[0] << ", " << trinodes[1] << ", " << trinodes[2] << ")\n";

        for (int j = 0; j < 3; ++j)
        {
            if (trinodes[j] < 0 || trinodes[j] >= NODE_COORD.size())
            {
                std::cerr << "Erro: índice inválido no elemento " << ielem << ": trinodes[" << j << "] = "
                          << trinodes[j] << " (NODE_COORD tem tamanho " << NODE_COORD.size() << ")\n";
                exit(EXIT_FAILURE);
            }
        }

        double x1 = NODE_COORD[trinodes[0]][0], y1 = NODE_COORD[trinodes[0]][1];
        double x2 = NODE_COORD[trinodes[1]][0], y2 = NODE_COORD[trinodes[1]][1];
        double x3 = NODE_COORD[trinodes[2]][0], y3 = NODE_COORD[trinodes[2]][1];

        // std::cout << "Calling sandt with nodes: " << "(" << x1 << ", " << y1 << "), " << "(" << x2 << ", " << y2 << "), " << "(" << x3 << ", " << y3 << ")\n";

        auto [S_elem, T_elem] = sandt(x1, y1, x2, y2, x3, y3); // retorna 3x3
                                                               // std::cout << "Saiu de sandt para o elemento [" << ielem << "]\n";

        // for (int jedge = 0; jedge < 3; ++jedge)
        // {
        //    for (int kedge = 0; kedge < 3; ++kedge)
        //    {
        // std::cout << S_elem[jedge][kedge] << "\t";
        // std::cout << T_elem[jedge][kedge] << "\t";
        //    }
        // std::cout << std::endl;
        // }

        for (int jedge = 0; jedge < 3; ++jedge)
        {
            int edge_j = ELEMENT_EDGES[ielem][jedge];
            int jj = dof_e1[edge_j];
            for (int kedge = 0; kedge < 3; ++kedge)
            {
                int edge_k = ELEMENT_EDGES[ielem][kedge];
                int kk = dof_e1[edge_k];
                if (jj >= 0 && jj >= 0)
                {
                    S[jj][kk] += S_elem[jedge][kedge];
                    T[jj][kk] += T_elem[jedge][kedge];
                }
            }
        }
    }

    MatrixXd S_mat(NUM_DOFS, NUM_DOFS), T_mat(NUM_DOFS, NUM_DOFS);
    for (int i = 0; i < NUM_DOFS; ++i)
    {
        for (int j = 0; j < NUM_DOFS; ++j)
        {
            std::cout << T[i][j] << "\t";
            S_mat(i, j) = S[i][j] / mu_r;
            T_mat(i, j) = T[i][j] * eps_r;
        }
        std::cout << std::endl;
    }
    // std::cout << "S_mat dimensions: " << S_mat.rows() << " x " << S_mat.cols() << std::endl;
    /*
        // Resolver o problema generalizado: (S/mu_r) x = lambda (T*eps_r) x
        GeneralizedEigenSolver<MatrixXd> solver(S_mat, T_mat);

        // Autovalores (podem ser complexos)
        VectorXcd eigvals = solver.eigenvalues();

        // Tirar raiz quadrada de cada autovalor
        std::vector<std::complex<double>> sqrt_vals;
        sqrt_vals.reserve(eigvals.size());
        for (int i = 0; i < eigvals.size(); ++i) {
            sqrt_vals.push_back(std::sqrt(eigvals[i]));
        }

        // Ordenar (por parte real, como no MATLAB sort)
        std::sort(sqrt_vals.begin(), sqrt_vals.end(),
                  [](const std::complex<double>& a, const std::complex<double>& b) {
                      return a.real() < b.real();
                  });

        // Imprimir resultado
        std::cout << "Autovalores (sqrt ordenados):\n";
        for (auto& val : sqrt_vals) {
            std::cout << val << "\n";
        }


        */

    GeneralizedSelfAdjointEigenSolver<MatrixXd> solver(S_mat, T_mat);
    write_csv("cpp_S_mat.csv", S_mat);
    write_csv("cpp_T_mat.csv", T_mat);

    VectorXd eigvals = solver.eigenvalues();
    MatrixXd eigvecs = solver.eigenvectors();

    write_csv("cpp_eigvals.csv", eigvals);
    write_csv("cpp_eigvecs.csv", eigvecs);

    MatrixXd V = eigvecs; // cópia
                          // for (int i = 0; i < V.cols(); ++i) {
    //     double nB = std::sqrt(V.col(i).transpose() * T_mat * V.col(i));
    //     if (nB > 0) V.col(i) /= nB;
    // }
    //  for (int i = 0; i < eigvals.size(); ++i)
    //      std::cout << "Autovalor[" << i << "]: " << eigvals[i] << std::endl;
    //  for (int i = 0; i < eigvals.size(); ++i)
    //     std::cout << "eig[" << i << "]: " << eigvals[i] << std::endl;

    for (int i = 0; i < eigvals.size(); ++i)
        std::cout << "kc[" << i << "]: " << std::sqrt(eigvals[i]) << std::endl;

    std::vector<std::pair<double, int>> indexed_kc;
    for (int i = 0; i < eigvals.size(); ++i)
        // if (eigvals[i] > 0)
        indexed_kc.emplace_back(std::sqrt(eigvals[i]), i);
    // std::sort(indexed_kc.begin(), indexed_kc.end());
    //  std::cout << "Indexed kc size: " << indexed_kc.size() << std::endl;

    // for (int i = 0; i < indexed_kc.size(); ++i)
    //    std::cout << "kc[" << i << "]: " << indexed_kc[i].first << std::endl;

    std::vector<int> node_flag;
    int num_free_nodes = 0;
    free_nodes(a, b, node_flag, num_free_nodes);

    int start_idx = num_free_nodes;
    int num_plot_modes = 6;
    std::cout << "num_free_nodes: " << num_free_nodes << std::endl;
    // Geração da grade de avaliação
    std::vector<double> XX, YY;
    int NX = 21, NY = 11;
    for (int i = 0; i <= NX; ++i)
        XX.push_back((double)i / NX * a);
    for (int j = 0; j <= NY; ++j)
        YY.push_back((double)j / NY * b);

    for (int ii = 0; ii < num_plot_modes; ++ii)
    {
        int idx = indexed_kc[start_idx + ii].second;

        std::cout << "idx: " << idx << std::endl;
        VectorXd eigmode = eigvecs.col(idx);
       // double maximo = 0.0;
       // for (int j = 0; j < eigmode.size(); ++j)
        //if(maximo < eigmode[j])
        //maximo = eigmode[j];

//for (int j = 0; j < eigmode.size(); ++j)
       // eigmode[j] /= maximo;


        for (int j = 0; j < eigmode.size(); ++j)
        std::cout << "eigmode["<< j <<"]: " << eigmode[j] << std::endl;
        
        std::vector<double> dofs(eigmode.data(), eigmode.data() + eigmode.size());

        plot_field(dofs, dof_e1, XX, YY, 3, 2, ii + 1, indexed_kc[start_idx + ii].first);
    }

    // for (int i = 0; i < V.col(start_idx).size(); ++i)
    //     std::cout << "eigvecs: " << V.col(start_idx)[i] << std::endl;
    //  Também plota modos espúrios (primeiros modos antes de free_nodes)
    //  for (int ii = 0; ii < num_plot_modes; ++ii)
    //  {
    //     int idx = indexed_kc[ii].second;
    //    VectorXd eigmode = eigvecs.col(idx);
    //     std::vector<double> dofs(eigmode.data(), eigmode.data() + eigmode.size());

    //    plot_field(dofs, dof_e1, XX, YY, 3, 2, ii + 1, indexed_kc[ii].first);
    // }
}

int main()
{
    Eigen2D(); // Chama sua função principal que faz toda a lógica

    /*
        // Definir malha
        NODE_COORD = {{0.0,0.0}, {1.0,0.0}, {0.0,1.0}};
        ELEMENTS   = {{0,1,2}}; // em C++ começa do 0

        double xp = 0.2, yp = 0.3;

        auto [lambda, area] = simplex2D(0, xp, yp);
        std::cout << "Lambda em C++:\n";
        for (double l : lambda) std::cout << l << " ";
        std::cout << "\nÁrea = " << area << "\n";

        auto W = whitney(0, xp, yp);
        std::cout << "Whitney em C++:\n";
        for (auto &row : W) {
            std::cout << row[0] << "\t" << row[1] << "\n";
        }

    */
    return 0;
}