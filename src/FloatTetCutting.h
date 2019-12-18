// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#ifndef FLOATTETWILD_FLOATTETCUTTING_H
#define FLOATTETWILD_FLOATTETCUTTING_H

#include <floattetwild/Mesh.hpp>
#include <floattetwild/AABBWrapper.h>

#ifdef FLOAT_TETWILD_USE_TBB
#include <tbb/concurrent_vector.h>
#endif

namespace floatTetWild {
    void cutting(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces, const std::vector<int> &input_tags,
                 Mesh &mesh, std::vector<bool> &is_face_inserted,
                 AABBWrapper& tree, bool is_again = false);

    bool insert_one(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
                    Mesh &sub_mesh, std::vector<std::array<std::vector<int>, 4>>& sub_cut_f_ids,
                    int f_id, int seed_v_id, bool is_parallel, bool is_again);

    void match_surface_fs(const Mesh &mesh, const std::vector<Vector3> &input_vertices,
                          const std::vector<Vector3i> &input_faces, std::vector<bool> &is_face_inserted,
                          std::vector<std::array<std::vector<int>, 4>> &cut_f_ids);

    void
    find_tets_for_cut(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces, int f_id,
                      int seed_v_id,
                      Mesh &mesh, std::vector<int> &intersection_results, std::vector<int> &oris);
    void
    find_tets_for_cut_again(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
            int f_id, Mesh &mesh, std::vector<int> &intersection_results, std::vector<int> &oris);


    int one_face_cut(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces, int f_id,
                     Mesh &mesh, std::vector<int> &intersection_results,
                     std::vector<std::array<std::vector<int>, 4>> &cut_f_ids, std::vector<int> &oris, bool is_parallel,
                     bool is_again);
    void snapping(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces, int f_id,
                  Mesh &mesh,
                  std::vector<int> &intersection_results, std::vector<int> &oris, std::vector<Scalar>& signed_dist);
    int tet_subdivision(int f_id,
                         Mesh& mesh, std::vector<int>& intersection_results, std::vector<int>& intersection_results_wn,
                         std::vector<std::array<std::vector<int>, 4>>& cut_f_ids, std::vector<Scalar>& signed_dist,std::vector<int>& oris,
                         std::map<std::array<int, 2>, int>& e_p_map, std::vector<MeshVertex>& tmp_tet_vertices,
                         bool is_again = false);
    Scalar get_min_volume(const Mesh &mesh, const std::vector<Vector4i>& new_tets, const std::vector<Vector3>& centroids,
            const std::vector<int>& global_v_ids, const std::vector<MeshVertex>& tmp_tet_vertices);
    void get_centroids(const Mesh &mesh, const std::vector<Vector4i>& new_tets, const std::vector<int>& global_v_ids,
            const std::vector<MeshVertex>& tmp_tet_vertices,  std::vector<Vector3>& cs);

    void preserve_cutting_open_boundary(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
        const std::vector<int>& current_inserted_cutting_fs,
        std::vector<std::array<std::vector<int>, 4>> &cut_f_ids, Mesh& mesh, AABBWrapper& tree,
        const std::vector<bool>& is_face_matched, std::vector<bool>& is_boundary_preserved, bool is_again);

    void track_surface(Mesh &mesh,
                       const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces, const std::vector<int> &input_tags,
                       std::vector<std::array<std::vector<int>, 4>> &cut_f_ids,
                       const AABBWrapper& tree, const std::vector<bool>& is_boundary_preserved);

    void get_current_open_boundary(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
            const std::vector<std::array<int, 2>>& edges, const std::vector<std::vector<int>>& conn_tris,
            std::vector<std::array<int, 2>>& f_id_le_id);
    void update_tmp_b_tree(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
                           const std::vector<bool> &is_face_inserted, AABBWrapper& tree);

    inline Vector3 get_normal(const Vector3& a, const Vector3& b, const Vector3& c) {
        return ((b - c).cross(a - c)).normalized();
    }
}
#endif //FLOATTETWILD_FLOATTETCUTTING_H
