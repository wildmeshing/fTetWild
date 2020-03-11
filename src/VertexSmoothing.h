// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#ifndef FLOATTETWILD_VERTEXSMOOTHING_H
#define FLOATTETWILD_VERTEXSMOOTHING_H

#include <floattetwild/Mesh.hpp>
#include <floattetwild/AABBWrapper.h>


namespace floatTetWild {
    void vertex_smoothing(Mesh& mesh, const AABBWrapper& tree);
    bool project_and_check(Mesh& mesh, int v_id, Vector3& p, const AABBWrapper& tree, bool is_sf, std::vector<Scalar>& new_qs);
    bool find_new_pos(Mesh& mesh, const int v_id, Vector3& p);
}

#endif //FLOATTETWILD_VERTEXSMOOTHING_H
