// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#ifndef FLOATTETWILD_SIMPLIFICATION_H
#define FLOATTETWILD_SIMPLIFICATION_H

#include <floattetwild/Mesh.hpp>
#include <floattetwild/AABBWrapper.h>

namespace floatTetWild {
    void simplify(std::vector<Vector3>& input_vertices, std::vector<Vector3i>& input_faces, std::vector<int>& input_tags,
            const AABBWrapper& tree, const Parameters& params, bool skip_simplify = false);
    bool remove_duplicates(std::vector<Vector3>& input_vertices, std::vector<Vector3i>& input_faces, std::vector<int>& input_tags, const Parameters& params);
    void collapsing(std::vector<Vector3>& input_vertices, std::vector<Vector3i>& input_faces, const AABBWrapper& sf_tree, const Parameters& params,
                    std::vector<bool>& is_v_removed, std::vector<bool>& is_f_removed, std::vector<std::unordered_set<int>>& conn_fs);
    void swapping(std::vector<Vector3>& input_vertices, std::vector<Vector3i>& input_faces, const AABBWrapper& sf_tree, const Parameters& params,
                  std::vector<bool>& is_v_removed, std::vector<bool>& is_f_removed, std::vector<std::unordered_set<int>>& conn_fs);
    void flattening(std::vector<Vector3>& input_vertices, std::vector<Vector3i>& input_faces, const AABBWrapper& sf_tree, const Parameters& params);

                    bool is_out_envelope(const std::array<Vector3, 3>& vs, const AABBWrapper& tree, const Parameters& params);
    Scalar get_angle_cos(const Vector3& p, const Vector3& p1, const Vector3& p2);

    void check_surface(std::vector<Vector3>& input_vertices, std::vector<Vector3i>& input_faces, const std::vector<bool>& is_f_removed,
                       const AABBWrapper& tree, const Parameters& params);

    void output_component(const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces, const std::vector<int>& input_tags);
}

#endif //FLOATTETWILD_SIMPLIFICATION_H
