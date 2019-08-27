#pragma once

#include <floattetwild/Parameters.h>
#include <floattetwild/Logger.hpp>

namespace floatTetWild {

int tetrahedralization(GEO::Mesh&       sf_mesh,
                       Parameters       params,
                       Eigen::MatrixXd& VO,
                       Eigen::MatrixXi& TO,
                       int              boolean_op    = -1,
                       bool             skip_simplify = false);

}
