#include "fdtd3d.hpp"

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
        << "Uso: ./build/fdtd3d [opcoes]\n"
        << "  --refine N              fator de refinamento\n"
        << "  --seed N                semente do ruido inicial\n"
        << "  --fft-output ARQ        nome do CSV de espectro em out/\n"
        << "  --time-output ARQ       nome do CSV do sinal temporal em out/\n"
        << "  --meta ARQ              nome do CSV de metadados em out/\n";
}
}

int main(int argc, char **argv) {
    Fdtd3DConfig config;

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
            if (arg == "--refine")
                config.refine = std::stoi(require_value(i, argc, argv));
            else if (arg == "--seed")
                config.random_seed = static_cast<unsigned int>(std::stoul(require_value(i, argc, argv)));
            else if (arg == "--fft-output")
                config.fft_output_filename = require_value(i, argc, argv);
            else if (arg == "--time-output")
                config.time_output_filename = require_value(i, argc, argv);
            else if (arg == "--meta")
                config.metadata_filename = require_value(i, argc, argv);
            else
                throw std::runtime_error("Opcao desconhecida: " + arg);
        }

        run_fdtd3d_simulation(config);
    }
    catch (const std::exception &ex)
    {
        std::cerr << "Erro: " << ex.what() << "\n";
        print_help();
        return 1;
    }

    return 0;
}
