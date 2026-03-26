#ifndef FDTD2D_PML_HPP
#define FDTD2D_PML_HPP

#include <string>

struct Fdtd2DPmlConfig
{
    bool cyl_present = false;
    bool line_source = false;
    int refine = 1;
    double pulse_compress = 1.0;
    int d_cell = 10;
    int poly_m = 3;
    int snapshot_step = 0;
    std::string output_filename = "ey_point1_pml.csv";
    std::string metadata_filename = "ey_point1_pml_meta.csv";
    std::string snapshot_prefix;
};

void run_fdtd2d_pml_simulation(const Fdtd2DPmlConfig &config = {});

#endif // FDTD2D_PML_HPP
