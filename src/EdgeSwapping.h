// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#ifndef FLOATTETWILD_EDGESWAPPING_H
#define FLOATTETWILD_EDGESWAPPING_H

#include <floattetwild/Mesh.hpp>

namespace floatTetWild {
    extern bool is_es_check;

    void edge_swapping(Mesh& mesh);
    bool remove_an_edge_32(Mesh& mesh, int v1_id, int v2_id, const std::vector<int>& n12_t_ids, std::vector<std::array<int, 2>>& new_edges);
    bool remove_an_edge_44(Mesh& mesh, int v1_id, int v2_id, const std::vector<int>& n12_t_ids, std::vector<std::array<int, 2>>& new_edges);
    bool remove_an_edge_56(Mesh& mesh, int v1_id, int v2_id, const std::vector<int>& n12_t_ids, std::vector<std::array<int, 2>>& new_edges);
}

#endif //FLOATTETWILD_EDGESWAPPING_H
