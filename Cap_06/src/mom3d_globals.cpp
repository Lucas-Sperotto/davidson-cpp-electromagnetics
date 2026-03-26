#include "mom3d_globals.hpp"

namespace cap06 {

Eigen::MatrixXi ELEMENTS;
Eigen::MatrixXd NODE_COORD;
Eigen::MatrixXi EDGES;
Eigen::MatrixXi ELEMENT_EDGES;
Eigen::MatrixXi EDGECONXELEMS;
Eigen::MatrixXi DOFLOCALNUM;
Eigen::MatrixXi ELEMENT_PLS_MNS;
Eigen::VectorXd ELL;
Eigen::Matrix<int, 3, 2> LOCALEDGENODES;
Eigen::Vector3i LOCALVERTEX;
int NUM_NODES = 0;
int NUM_ELEMENTS = 0;
int NUM_EDGES = 0;
int NUM_DOFS = 0;

std::filesystem::path project_out_dir()
{
    return std::filesystem::path(PROJECT_OUT_DIR);
}

Complex plain_dot(const Vector3cd& lhs, const Vector3d& rhs)
{
    return lhs(0) * rhs(0) + lhs(1) * rhs(1) + lhs(2) * rhs(2);
}

void initialize_mom3d_globals()
{
    LOCALEDGENODES << 1, 2,
                      0, 2,
                      0, 1;
    LOCALVERTEX << 0, 1, 2;

    ELEMENTS.resize(0, 3);
    NODE_COORD.resize(0, 3);
    EDGES.resize(0, 2);
    ELEMENT_EDGES.resize(0, 3);
    EDGECONXELEMS.resize(0, 2);
    DOFLOCALNUM.resize(0, 2);
    ELEMENT_PLS_MNS.resize(0, 3);
    ELL.resize(0);
    NUM_NODES = 0;
    NUM_ELEMENTS = 0;
    NUM_EDGES = 0;
    NUM_DOFS = 0;
}

}  // namespace cap06
