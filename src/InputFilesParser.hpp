// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Teseo Schneider <teseo.schneider@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#pragma once

#include <floattetwild/Types.hpp>

#include <geogram/mesh/mesh.h>

#include <vector>
#include <map>
#include <string>

namespace floatTetWild {

class InputFilesParser {
public:
    void get_meshes(const json &input_files, std::vector<std::string> &meshes)
    {
        meshes.clear();

        get_meshes_aux(input_files, meshes);
    }

    bool load_and_merge(const std::vector<std::string> &meshes, std::vector<Vector3> &V, std::vector<Vector3i> &F, GEO::Mesh &sf_mesh, std::vector<int> &tags);
    void merge(const std::vector<std::vector<Vector3>> &Vs, const std::vector<std::vector<Vector3i>> &Fs, std::vector<Vector3> &V, std::vector<Vector3i> &F, GEO::Mesh &sf_mesh, std::vector<int> &tags);

    std::vector<Vector3> bbox_mins;
    std::vector<Vector3> bbox_maxes;
    std::vector<Scalar>  bbox_diag_lengths;
    std::vector<Scalar>  target_edge_lengths;

private:
    void get_meshes_aux(const json &input_files_node, std::vector<std::string> &meshes);
};
}  // namespace floatTetWild
