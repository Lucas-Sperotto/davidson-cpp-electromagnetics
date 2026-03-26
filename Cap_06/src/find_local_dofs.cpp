#include "find_local_dofs.hpp"

#include "mom3d_globals.hpp"

namespace cap06 {

void find_local_dofs(const Eigen::VectorXi& dof_RWG)
{
    DOFLOCALNUM = Eigen::MatrixXi::Constant(NUM_DOFS, 2, -1);

    for (int iedge = 0; iedge < NUM_EDGES; ++iedge) {
        bool edge_found = false;
        bool terminate = false;

        for (int jelem = 0; jelem < NUM_ELEMENTS; ++jelem) {
            for (int kedge = 0; kedge < 3; ++kedge) {
                if (ELEMENT_EDGES(jelem, kedge) == iedge) {
                    const int dof_id = dof_RWG(iedge);
                    if (!edge_found) {
                        if (dof_id > 0) {
                            DOFLOCALNUM(dof_id - 1, 0) = kedge;
                        }
                        edge_found = true;
                    } else {
                        if (dof_id > 0) {
                            DOFLOCALNUM(dof_id - 1, 1) = kedge;
                        }
                        terminate = true;
                        break;
                    }
                }
            }

            if (terminate) {
                break;
            }
        }
    }
}

}  // namespace cap06
