// main.cpp
#include "include/fdtd2d.hpp"
#include "include/fdtd2d_pml.hpp"

int main() {
//#ifdef RUN_PML
    run_fdtd2d_pml_simulation();
//#else
    run_fdtd2d_simulation();
//#endif
    return 0;
}