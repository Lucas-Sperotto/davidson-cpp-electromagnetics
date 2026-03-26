#pragma once

#include "mom3d_globals.hpp"

namespace cap06 {

struct SingularIntegral {
    Complex value;
    Complex xi1;
    Complex xi2;
    Complex xi3;
};

SingularIntegral intg_sing_SGF(
    double k,
    const Vector3d& r,
    const Vector3d& r1,
    const Vector3d& r2,
    const Vector3d& r3,
    int num_pts_rad,
    int num_pts_trans);

}  // namespace cap06
