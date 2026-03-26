// include/fdtd3d.hpp
#ifndef FDTD3D_HPP
#define FDTD3D_HPP

#include <string>

struct Fdtd3DConfig
{
    int refine = 2;
    unsigned int random_seed = 12345;
    std::string fft_output_filename = "hz_center_fft.csv";
    std::string time_output_filename = "hz_center_time.csv";
    std::string metadata_filename = "hz_center_meta.csv";
};

void run_fdtd3d_simulation(const Fdtd3DConfig &config = {});

#endif // FDTD3D_HPP
