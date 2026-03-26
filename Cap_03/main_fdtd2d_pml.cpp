#include "include/fdtd2d_pml.hpp"

#include <iostream>
#include <stdexcept>
#include <string>

namespace
{
std::string require_value(int &index, int argc, char **argv)
{
    if (index + 1 >= argc)
        throw std::runtime_error(std::string("Falta valor para ") + argv[index]);
    return argv[++index];
}

void print_help()
{
    std::cout
        << "Uso: ./build/fdtd2d_pml [opcoes]\n"
        << "  --cyl-present           inclui cilindro PEC\n"
        << "  --line-source           usa fonte de linha no lugar da interface scat/tot\n"
        << "  --refine N              fator de refinamento\n"
        << "  --pulse-compress X      compressao do pulso gaussiano\n"
        << "  --d-cell N              espessura da PML em celulas\n"
        << "  --poly-m N              ordem do perfil polinomial da PML\n"
        << "  --output ARQ            nome do CSV de saida em out/\n"
        << "  --meta ARQ              nome do CSV de metadados em out/\n";
}
}

int main(int argc, char **argv) {
    Fdtd2DPmlConfig config;

    try
    {
        for (int i = 1; i < argc; ++i)
        {
            std::string arg = argv[i];
            if (arg == "--help" || arg == "-h")
            {
                print_help();
                return 0;
            }
            if (arg == "--cyl-present")
                config.cyl_present = true;
            else if (arg == "--line-source")
                config.line_source = true;
            else if (arg == "--refine")
                config.refine = std::stoi(require_value(i, argc, argv));
            else if (arg == "--pulse-compress")
                config.pulse_compress = std::stod(require_value(i, argc, argv));
            else if (arg == "--d-cell")
                config.d_cell = std::stoi(require_value(i, argc, argv));
            else if (arg == "--poly-m")
                config.poly_m = std::stoi(require_value(i, argc, argv));
            else if (arg == "--output")
                config.output_filename = require_value(i, argc, argv);
            else if (arg == "--meta")
                config.metadata_filename = require_value(i, argc, argv);
            else
                throw std::runtime_error("Opcao desconhecida: " + arg);
        }

        run_fdtd2d_pml_simulation(config);
    }
    catch (const std::exception &ex)
    {
        std::cerr << "Erro: " << ex.what() << "\n";
        print_help();
        return 1;
    }

    return 0;
}
