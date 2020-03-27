// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Teseo Schneider <teseo.schneider@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include <floattetwild/InputFilesParser.hpp>

#include <floattetwild/MeshIO.hpp>
#include <floattetwild/Logger.hpp>
#include <floattetwild/AABBWrapper.h>

namespace floatTetWild {
    void InputFilesParser::get_meshes_aux(const json &input_files_node, std::vector<std::string> &meshes)
    {
        int index = 0;
        std::map<std::string, int> existings;

        auto files = input_files_node["files"];
        for(const auto &file: files) {
            const std::string name = file["name"];
            const auto iter = existings.find(name);
            if(iter == existings.end()) {
                meshes.push_back(name);

                Scalar target_edge_length = -1;
                if (file["target_edge_length"].is_number())
                    target_edge_length = file["target_edge_length"];
                
                target_edge_lengths.push_back(target_edge_length);

                existings[name] = index++;
            }
        }    
    }

    bool InputFilesParser::load_and_merge(const std::vector<std::string> &meshes, std::vector<Vector3> &V, std::vector<Vector3i> &F, GEO::Mesh &sf_mesh, std::vector<int> &tags)
    {
        bbox_mins.clear();
        bbox_maxes.clear();

        std::vector<std::vector<Vector3>> Vs;
        std::vector<std::vector<Vector3i>> Fs;

        Vs.resize(meshes.size());
        Fs.resize(meshes.size());

        GEO::Mesh tmp_mesh;
        std::vector<int> tmp_tags;

        for(int i = 0; i < meshes.size(); ++i) {
            const auto &m = meshes[i];
            if (!MeshIO::load_mesh(m, Vs[i], Fs[i], tmp_mesh, tmp_tags)) {
                logger().error("unable to open {} file", m);
                return false;
            }
            AABBWrapper tree(tmp_mesh);

            Vector3 bbox_min, bbox_max;
            tree.get_bbox(bbox_min, bbox_max);
            bbox_mins.push_back(bbox_min);
            bbox_maxes.push_back(bbox_max);
        }

        merge(Vs, Fs, V, F, sf_mesh, tags);
        return true;
    }

    void InputFilesParser::merge(const std::vector<std::vector<Vector3>> &Vs, const std::vector<std::vector<Vector3i>> &Fs, std::vector<Vector3> &V, std::vector<Vector3i> &F, GEO::Mesh &sf_mesh, std::vector<int> &tags)
    {
        V.clear();
        F.clear();

        for(const auto &vv : Vs)
            V.insert(V.end(), vv.begin(), vv.end());

        int offset = 0;
        int size = 0;
        for(int id = 0; id < Fs.size(); ++id) {
            const auto &ff = Fs[id];

            for(const auto fid : ff)
                F.push_back(Vector3i(fid(0)+offset, fid(1)+offset, fid(2)+offset));

            tags.insert(tags.begin()+size, Fs[id].size(), id);
            size += Fs[id].size();
            offset += Vs[id].size();
        }

        MeshIO::load_mesh(V, F, sf_mesh, tags);
    }
}
