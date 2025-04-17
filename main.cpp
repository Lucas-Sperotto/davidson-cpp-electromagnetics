#include <iostream>
#include <vector>
#include "include/free_nodes.h"
#include "include/trimesh.h"

int main() {
    std::cout << "Teste do Capítulo 10: Geração de malha e nós livres.\n";

    // Parâmetros da malha
    double a = 0.02286;
    double b = 0.01016;
    int Nx = 4, Ny = 2;

    // Geração da malha
    std::vector<std::vector<double>> x_nodes, y_nodes;
    trimesh(a, b, Nx, Ny, x_nodes, y_nodes);
    std::cout << "Malha gerada com sucesso: " << x_nodes[0].size() << " pontos x, "
              << y_nodes[0].size() << " pontos y.\n";

    // Identificação de nós livres
    std::vector<int> node_flag;
    int num_free_nodes = 0;
    free_nodes(a, b, node_flag, num_free_nodes);
    std::cout << "Número de nós livres detectados: " << num_free_nodes << "\n";

    return 0;
}
