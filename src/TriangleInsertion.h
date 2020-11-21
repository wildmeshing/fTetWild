// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

//
// Created by Yixin Hu on 2019-08-27.
//

#ifndef FLOATTETWILD_TRIANGLEINSERTION_H
#define FLOATTETWILD_TRIANGLEINSERTION_H

#include <floattetwild/Mesh.hpp>
#include <floattetwild/AABBWrapper.h>
#include <floattetwild/CutMesh.h>

#include <floattetwild/Rational.h>
namespace floatTetWild {
    void match_surface_fs(const Mesh &mesh,
                          const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
                          std::vector<bool> &is_face_inserted,
                          std::vector<std::array<std::vector<int>, 4>> &track_surface_fs);

    void insert_triangles(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
                          const std::vector<int> &input_tags, Mesh &mesh,
                          std::vector<bool> &is_face_inserted, AABBWrapper &tree, bool is_again);

    void insert_triangles_aux(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
                              const std::vector<int> &input_tags, Mesh &mesh,
                              std::vector<bool> &is_face_inserted, AABBWrapper &tree, bool is_again);

    void optimize_non_surface(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
                              const std::vector<int> &input_tags, std::vector<bool> &is_face_inserted,
                              const std::vector<std::array<std::vector<int>, 4 >> &track_surface_fs,
                              Mesh &mesh, AABBWrapper &tree, bool is_again);

    void sort_input_faces(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
                          const Mesh &mesh, std::vector<int> &sorted_f_ids);

    void push_new_tets(Mesh &mesh, std::vector<std::array<std::vector<int>, 4>> &track_surface_fs,
                       std::vector<Vector3> &points, std::vector<MeshTet> &new_tets,
                       std::vector<std::array<std::vector<int>, 4>> &new_track_surface_fs,
                       std::vector<int> &modified_t_ids, bool is_again);

    void simplify_subdivision_result(int insert_f_id, int input_v_size, Mesh &mesh, AABBWrapper &tree,
                                     std::vector<std::array<std::vector<int>, 4>> &track_surface_fs);

    ///face
    bool insert_multi_triangles(int f_id, const std::vector<Vector3> &input_vertices,
                               const std::vector<Vector3i> &input_faces, const std::vector<int> &input_tags,
                                const std::vector<std::vector<int>>& conn_fs,
                                const std::vector<Vector3> &ns, std::vector<bool> & is_visited, std::vector<int>& f_ids,
                               Mesh &mesh, std::vector<std::array<std::vector<int>, 4>> &track_surface_fs,
                               AABBWrapper &tree, bool is_again);

    bool insert_one_triangle(int f_id, const std::vector<Vector3> &input_vertices,
                             const std::vector<Vector3i> &input_faces, const std::vector<int> &input_tags,
                             Mesh &mesh, std::vector<std::array<std::vector<int>, 4>> &track_surface_fs,
                             AABBWrapper &tree, bool is_again);

    void
    find_cutting_tets(int f_id, const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
                      const std::array<Vector3, 3> &vs, Mesh &mesh, std::vector<int> &result, bool is_again);

    bool subdivide_tets(int insert_f_id, Mesh &mesh, CutMesh &cut_mesh, std::vector<Vector3> &points,
                        std::map<std::array<int, 2>, int> &map_edge_to_intersecting_point,
                        std::vector<std::array<std::vector<int>, 4>> &track_surface_fs,
                        std::vector<int> &subdivide_t_ids, std::vector<bool> &is_mark_surface,
                        std::vector<MeshTet> &new_tets,
                        std::vector<std::array<std::vector<int>, 4>> &new_track_surface_fs,
                        std::vector<int> &modified_t_ids);

    void pair_track_surface_fs(Mesh &mesh, std::vector<std::array<std::vector<int>, 4>> &track_surface_fs);

    ///edge
    void find_boundary_edges(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
                             const std::vector<bool> &is_face_inserted, const std::vector<bool> &old_is_face_inserted,
                             std::vector<std::pair<std::array<int, 2>, std::vector<int>>> &b_edge_infos,
                             std::vector<bool> &is_on_cut_edges,
                             std::vector<std::array<int, 2>> &b_edges);

    bool insert_boundary_edges(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
                               std::vector<std::pair<std::array<int, 2>, std::vector<int>>> &b_edge_infos,
                               std::vector<bool> &is_on_cut_edges,
                               std::vector<std::array<std::vector<int>, 4>> &track_surface_fs, Mesh &mesh,
                               AABBWrapper &tree,
                               std::vector<bool> &is_face_inserted, bool is_again,
                               std::vector<std::array<int, 3>> &known_surface_fs,
                               std::vector<std::array<int, 3>> &known_not_surface_fs);

    bool insert_boundary_edges_get_intersecting_edges_and_points(
            const std::vector<std::vector<std::pair<int, int>>> &covered_fs_infos,
            const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
            const std::array<int, 2> &e, const std::vector<int> &n_f_ids,
            std::vector<std::array<std::vector<int>, 4>> &track_surface_fs, Mesh &mesh,
            std::vector<Vector3> &points, std::map<std::array<int, 2>, int> &map_edge_to_intersecting_point,
            std::vector<int> &snapped_v_ids, std::vector<std::array<int, 3>> &cut_fs,
            bool is_again);

    ///other
    bool is_uninserted_face_covered(int uninserted_f_id, const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
                                    const std::vector<int>& cut_t_ids, Mesh &mesh);

    void mark_surface_fs(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
                         const std::vector<int> &input_tags,
                         std::vector<std::array<std::vector<int>, 4>> &track_surface_fs,
                         const std::vector<bool> &is_face_inserted,
                         const std::vector<std::array<int, 3>> &known_surface_fs,
                         const std::vector<std::array<int, 3>> &known_not_surface_fs,
                         std::vector<std::array<int, 2>> &b_edges,
                         Mesh &mesh, AABBWrapper &tree);

    int get_opp_t_id(int t_id, int j, const Mesh &mesh);

    void myassert(bool b, const std::string &s);

    void check_track_surface_fs(Mesh &mesh, std::vector<std::array<std::vector<int>, 4>> &track_surface_fs,
                                const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
                                const std::vector<int> &sorted_f_ids);
    Vector3 get_normal(const Vector3& a, const Vector3& b, const Vector3& c);

    //fortest
    typedef Eigen::Matrix<triwild::Rational, 3, 1> Vector3_r;

    int orient_rational(const Vector3_r &p1, const Vector3_r &p2, const Vector3_r &p3, const Vector3_r &p);
}


#endif //FLOATTETWILD_TRIANGLEINSERTION_H