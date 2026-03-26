#pragma once

#include "mom3d_globals.hpp"

namespace cap06 {

struct IntegralPQ {
    Complex Ipq;
    Complex Ipq_xi;
    Complex Ipq_eta;
    Complex Ipq_zeta;
};

IntegralPQ Int_pq(int p, int q, const Vector3d& r_cp, double k, int quad_pts, bool sing);

}  // namespace cap06
