// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#ifndef FLOATTETWILD_EDGESPLITTING_H
#define FLOATTETWILD_EDGESPLITTING_H

#include <floattetwild/Mesh.hpp>
#include <floattetwild/AABBWrapper.h>

namespace floatTetWild {
    void edge_splitting(Mesh &mesh, const AABBWrapper& tree);

    bool split_an_edge(Mesh &mesh, int v1_id, int v2_id, bool is_repush, std::vector<std::array<int, 2>> &new_edges,
                       std::vector<bool> &is_splittable, const AABBWrapper& tree);
}

#endif //FLOATTETWILD_EDGESPLITTING_H
