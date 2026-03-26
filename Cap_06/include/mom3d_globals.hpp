#pragma once

#include <Eigen/Dense>

#include <complex>
#include <filesystem>

namespace cap06 {

using Complex = std::complex<double>;
using Vector3d = Eigen::Vector3d;
using Vector3cd = Eigen::Matrix<Complex, 3, 1>;

extern Eigen::MatrixXi ELEMENTS;
extern Eigen::MatrixXd NODE_COORD;
extern Eigen::MatrixXi EDGES;
extern Eigen::MatrixXi ELEMENT_EDGES;
extern Eigen::MatrixXi EDGECONXELEMS;
extern Eigen::MatrixXi DOFLOCALNUM;
extern Eigen::MatrixXi ELEMENT_PLS_MNS;
extern Eigen::VectorXd ELL;
extern Eigen::Matrix<int, 3, 2> LOCALEDGENODES;
extern Eigen::Vector3i LOCALVERTEX;
extern int NUM_NODES;
extern int NUM_ELEMENTS;
extern int NUM_EDGES;
extern int NUM_DOFS;

constexpr double PI = 3.141592653589793238462643383279502884;

inline Complex j_complex()
{
    return Complex(0.0, 1.0);
}

std::filesystem::path project_out_dir();
Complex plain_dot(const Vector3cd& lhs, const Vector3d& rhs);
void initialize_mom3d_globals();

}  // namespace cap06
