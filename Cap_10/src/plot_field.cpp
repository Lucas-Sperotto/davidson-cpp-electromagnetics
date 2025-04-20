#include <vector>
#include <iostream>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "globals.h"
#include "simplex2D.h"
#include "whitney.h"

// Avalia e armazena o campo vetorial (para visualização posterior)
void plot_field(const std::vector<double> &dofs,
                const std::vector<int> &dof_e1,
                const std::vector<double> &XX,
                const std::vector<double> &YY,
                int rows, int cols, int plotnum,
                double eigvalue)
{

    int Nx = XX.size();
    int Ny = YY.size();

    std::vector<std::vector<std::vector<double>>> vec_field(Nx,
                                                            std::vector<std::vector<double>>(Ny, std::vector<double>(2, 0.0)));

    for (int ii = 0; ii < Nx; ++ii)
    {
        for (int jj = 0; jj < Ny; ++jj)
        {
            double xc = XX[ii];
            double yc = YY[jj];

            bool point_found = false;

            for (int i_elem = 0; i_elem < NUM_ELEMENTS && !point_found; ++i_elem)
            {
                const auto &nodes = ELEMENTS[i_elem];
                double x1 = NODE_COORD[nodes[0]][0], y1 = NODE_COORD[nodes[0]][1];
                double x2 = NODE_COORD[nodes[1]][0], y2 = NODE_COORD[nodes[1]][1];
                double x3 = NODE_COORD[nodes[2]][0], y3 = NODE_COORD[nodes[2]][1];

                double x_min = std::min(std::min(x1, x2), x3);
                double x_max = std::max(std::max(x1, x2), x3);
                double y_min = std::min(std::min(y1, y2), y3);
                double y_max = std::max(std::max(y1, y2), y3);

                if (xc >= x_min && xc <= x_max && yc >= y_min && yc <= y_max)
                {
                    std::vector<double> lambda = simplex2D(i_elem, xc, yc).first; // coordenadas baricêntricas
                    double tol = 1e-4;

                    if (lambda[0] >= -tol && lambda[0] <= 1 + tol &&
                        lambda[1] >= -tol && lambda[1] <= 1 + tol &&
                        lambda[2] >= -tol && lambda[2] <= 1 + tol)
                    {

                        std::vector<int> edges = ELEMENT_EDGES[i_elem];
                        std::vector<double> dofs_tri(3, 0.0);

                        for (int k = 0; k < 3; ++k)
                        {
                            int global_dof = dof_e1[edges[k]];
                            if (global_dof > 0)
                                dofs_tri[k] = dofs[global_dof];
                        }

                        std::vector<std::vector<double>> w = whitney(i_elem, xc, yc);
                        for (int k = 0; k < 3; ++k)
                        {
                            vec_field[ii][jj][0] += dofs_tri[k] * w[k][0];
                            vec_field[ii][jj][1] += dofs_tri[k] * w[k][1];
                        }

                        point_found = true;
                    }
                }
            }

            if (!point_found)
            {
                std::cerr << "Warning: point (" << xc << "," << yc << ") not found in any triangle\n";
            }
        }
    }

    // Cria pasta "out" se não existir
    std::filesystem::create_directory("../out");

    // Formata o nome do arquivo
    std::ostringstream filename;
    filename << "../out/field_kc_" << std::fixed << std::setprecision(6) << eigvalue << ".csv";

    // Abre o arquivo
    std::ofstream fout(filename.str());
    if (!fout.is_open())
    {
        std::cerr << "Erro ao abrir arquivo: " << filename.str() << "\n";
        return;
    }

    // Escreve cabeçalho
    fout << "x,y,Ex,Ey\n";

    // Salva os dados
    for (int ii = 0; ii < Nx; ++ii)
    {
        for (int jj = 0; jj < Ny; ++jj)
        {
            fout << XX[ii] << "," << YY[jj] << ","
                 << vec_field[ii][jj][0] << "," << vec_field[ii][jj][1] << "\n";
        }
    }

    fout.close();
    std::cout << "Campo avaliado e salvo em: " << filename.str() << "\n";

    std::cout << "Campo avaliado para k_c = " << eigvalue << "\n";
}
