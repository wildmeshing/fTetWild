// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

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
