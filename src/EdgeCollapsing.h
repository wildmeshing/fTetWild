// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#ifndef FLOATTETWILD_EDGECOLLAPSING_H
#define FLOATTETWILD_EDGECOLLAPSING_H

#include <floattetwild/Mesh.hpp>
#include <floattetwild/AABBWrapper.h>

namespace floatTetWild {
    void edge_collapsing(Mesh& mesh, const AABBWrapper& tree);
    int collapse_an_edge(Mesh& mesh, int v1_id, int v2_id, const AABBWrapper& tree,
            std::vector<std::array<int, 2>>& new_edges, int ts, std::vector<int>& tet_tss,
            bool is_check_quality = true, bool is_update_tss = true);

    bool is_edge_freezed(Mesh& mesh, int v1_id, int v2_id);
    bool is_collapsable_bbox(Mesh& mesh, int v1_id, int v2_id);
    bool is_collapsable_length(Mesh& mesh, int v1_id, int v2_id, Scalar l_2);
    bool is_collapsable_boundary(Mesh& mesh, int v1_id, int v2_id, const AABBWrapper& tree);
}

#endif //FLOATTETWILD_EDGECOLLAPSING_H
