// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#ifndef FLOATTETWILD_MESHIMPROVEMENT_H
#define FLOATTETWILD_MESHIMPROVEMENT_H

#include <floattetwild/Mesh.hpp>
#include <floattetwild/AABBWrapper.h>
#include <floattetwild/Types.hpp>

namespace floatTetWild {
    void init(Mesh &mesh, AABBWrapper& tree);
    void optimization(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces, const std::vector<int> &input_tags, std::vector<bool> &is_face_inserted,
            Mesh &mesh, AABBWrapper& tree, const std::array<int, 4> &ops = {{1, 1, 1, 1}});
    void cleanup_empty_slots(Mesh &mesh, double percentage = 0.7);
    void operation(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces, const std::vector<int> &input_tags, std::vector<bool> &is_face_inserted,
            Mesh &mesh, AABBWrapper& tree, const std::array<int, 5> &ops = {{1, 1, 1, 1, 1}});
    void operation(Mesh &mesh, AABBWrapper& tree, const std::array<int, 4> &ops = {{1, 1, 1, 1}});
    bool update_scaling_field(Mesh &mesh, Scalar max_energy);

    int get_max_p(const Mesh &mesh);

    void correct_tracked_surface_orientation(Mesh &mesh, AABBWrapper& tree);
    void get_tracked_surface(Mesh& mesh, Eigen::Matrix<Scalar, Eigen::Dynamic, 3> &V, Eigen::Matrix<int, Eigen::Dynamic, 3> &F, int c_id = 0);
    void boolean_operation(Mesh& mesh, int op);
    void boolean_operation(Mesh& mesh, const json& csg_tree_with_ids, const std::vector<std::string> &meshes);
    void boolean_operation(Mesh& mesh, const json &csg_tree_with_ids, const std::vector<Eigen::VectorXd> &w);
    void boolean_operation(Mesh& mesh, const json &csg_tree_with_ids);
    void filter_outside(Mesh& mesh, const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces);
    void filter_outside(Mesh& mesh, bool invert_faces = false);
    void filter_outside_floodfill(Mesh& mesh, bool invert_faces = false);
    void mark_outside(Mesh& mesh, bool invert_faces = false);
    void smooth_open_boundary(Mesh& mesh, const AABBWrapper& tree);
    void smooth_open_boundary_aux(Mesh& mesh, const AABBWrapper& tree);
    void get_surface(Mesh& mesh, Eigen::MatrixXd& V, Eigen::MatrixXi& F);
    void manifold_surface(Mesh& mesh, Eigen::MatrixXd& V, Eigen::MatrixXi& F);
    void manifold_edges(Mesh& mesh);
    void manifold_vertices(Mesh& mesh);

//    void apply_sizingfield(const Eigen::VectorXd& V_in, const Eigen::VectorXi& T_in, const Eigen::VectorXd& values, Mesh& mesh, AABBWrapper& tree);
    void apply_sizingfield(Mesh& mesh, AABBWrapper& tree);
    void apply_coarsening(Mesh& mesh, AABBWrapper& tree);

    void output_info(Mesh& mesh, const AABBWrapper& tree);
    void check_envelope(Mesh& mesh, const AABBWrapper& tree);
    void output_surface(Mesh& mesh, const std::string& filename);

    void untangle(Mesh &mesh);
}

#endif //FLOATTETWILD_MESHIMPROVEMENT_H
