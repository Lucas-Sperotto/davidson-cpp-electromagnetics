#include "edge_conx_elem.hpp"

#include "mom3d_globals.hpp"

namespace cap06 {

void edge_conx_elem(const Eigen::VectorXi& dof_RWG)
{
    EDGECONXELEMS = Eigen::MatrixXi::Constant(NUM_DOFS, 2, -1);
    ELEMENT_PLS_MNS = Eigen::MatrixXi::Zero(NUM_ELEMENTS, 3);

    for (int iedge = 0; iedge < NUM_EDGES; ++iedge) {
        const int dof_id = dof_RWG(iedge);
        if (dof_id > 0) {
            int counter = 0;
            for (int jelem = 0; jelem < NUM_ELEMENTS; ++jelem) {
                for (int kedge = 0; kedge < 3; ++kedge) {
                    if (ELEMENT_EDGES(jelem, kedge) == iedge) {
                        EDGECONXELEMS(dof_id - 1, counter) = jelem;
                        ELEMENT_PLS_MNS(jelem, kedge) = (counter == 0) ? +1 : -1;
                        ++counter;
                    }
                }
            }
        }
    }
}

}  // namespace cap06
