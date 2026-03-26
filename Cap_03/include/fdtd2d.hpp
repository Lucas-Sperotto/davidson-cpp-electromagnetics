// include/fdtd2d.hpp
#ifndef FDTD2D_HPP
#define FDTD2D_HPP

#include <string>

struct Fdtd2DConfig
{
    bool cyl_present = false;
    int refine = 1;
    double pulse_compress = 1.0;
    int snapshot_step = 0;
    std::string output_filename = "ey_point1.csv";
    std::string metadata_filename = "ey_point1_meta.csv";
    std::string snapshot_prefix;
};

void run_fdtd2d_simulation(const Fdtd2DConfig &config = {});

#endif // FDTD2D_HPP
