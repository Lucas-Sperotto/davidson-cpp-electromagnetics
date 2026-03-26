#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace fs = std::filesystem;

int main()
{
    // Translation of FETD_FDTD.m:
    // analyze CTLN functions for the FETD-FDTD equivalence note.
    const fs::path out_dir = PROJECT_OUT_DIR;
    fs::create_directories(out_dir);

    std::ofstream out(out_dir / "fetd_fdtd_ctlN_funcs.csv");
    out << "y,we1,we2,we1_dot_we2,is_quad_point,quad_value\n";
    out << std::scientific << std::setprecision(10);

    for (int index = 0; index <= 100; ++index)
    {
        const double y = index * 0.01;
        const double we1 = 1.0 - y;
        const double we2 = y;
        const double we1_dot_we2 = we1 * we2;
        const int is_quad_point = (index == 0 || index == 100) ? 1 : 0;
        const double quad_value = 0.0;
        out << y << "," << we1 << "," << we2 << ","
            << we1_dot_we2 << "," << is_quad_point << "," << quad_value << "\n";
    }

    std::ofstream summary(out_dir / "fetd_fdtd_summary.csv");
    summary << "key,value\n";
    summary << "quad_point_1,0\n";
    summary << "quad_point_2,1\n";
    summary << "quad_value_1,0\n";
    summary << "quad_value_2,0\n";

    std::cout << "Arquivos gerados: fetd_fdtd_ctlN_funcs.csv e fetd_fdtd_summary.csv\n";
    return 0;
}
