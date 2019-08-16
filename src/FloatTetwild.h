#pragma once

#include <floattetwild/Parameters.h>
#include <igl/Timer.h>
#include <floattetwild/Logger.hpp>
#include <floattetwild/MeshIO.hpp>
#include <ifstream>

namespace floatTetWild {

int tetrahedralization(GEO::Mesh&       sf_mesh,
                       Parameters       params,
                       Eigen::MatrixXd& VO,
                       Eigen::MatrixXi& TO,
                       int              boolean_op    = -1,
                       bool             skip_simplify = false);

}
