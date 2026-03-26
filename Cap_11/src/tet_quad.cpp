#include "../include/tet_quad.hpp"

#include <Eigen/Dense>

#include <stdexcept>

// Translation of tet_quad.m:
// Keast tetrahedral quadrature rules.
std::pair<Eigen::VectorXd, Eigen::MatrixXd> tet_quad(int n)
{
    if (n != 11)
        throw std::runtime_error("tet_quad: regra nao implementada.");

    const double zz = 1.0 / 4.0;
    const double aa = 0.714285714285714285e-01;
    const double bb = 0.785714285714285714;
    const double cc = 0.399403576166799219;
    const double dd = 0.100596423833200785;
    const double w1 = -0.1315555555555555550e-01;
    const double w2 = 0.7622222222222222222e-02;
    const double w3 = 0.2488888888888888880e-01;

    Eigen::MatrixXd lambda(11, 4);
    lambda << zz, zz, zz, zz,
        aa, aa, aa, bb,
        aa, aa, bb, aa,
        aa, bb, aa, aa,
        bb, aa, aa, aa,
        cc, cc, dd, dd,
        cc, dd, cc, dd,
        cc, dd, dd, cc,
        dd, dd, cc, cc,
        dd, cc, dd, cc,
        dd, cc, cc, dd;

    Eigen::VectorXd w(11);
    w << w1, w2, w2, w2, w2, w3, w3, w3, w3, w3, w3;
    w *= 6.0;

    return {w, lambda};
}
