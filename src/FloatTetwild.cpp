// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include <floattetwild/AABBWrapper.h>
#include <floattetwild/FloatTetDelaunay.h>
#include <floattetwild/MeshImprovement.h>
#include <floattetwild/Parameters.h>
#include <floattetwild/Simplification.h>
#include <floattetwild/TriangleInsertion.h>
#include <floattetwild/LocalOperations.h>
#include <floattetwild/Statistics.h>
#include <geogram/mesh/mesh_reorder.h>
#include <geogram/mesh/mesh_repair.h>
#include <igl/Timer.h>
#include <floattetwild/Logger.hpp>
#include <floattetwild/MeshIO.hpp>
#include <floattetwild/Types.hpp>
#include <fstream>

namespace floatTetWild {

int tetrahedralization(GEO::Mesh&       sf_mesh,
                       Parameters       params,
                       Eigen::MatrixXd& V,
                       Eigen::MatrixXi& T,
                       int              boolean_op,
                       bool             skip_simplify)
{
    if (!sf_mesh.facets.are_simplices()) {
        GEO::mesh_repair(
          sf_mesh, GEO::MeshRepairMode(GEO::MESH_REPAIR_TRIANGULATE | GEO::MESH_REPAIR_QUIET));
    }
    GEO::mesh_reorder(sf_mesh, GEO::MESH_ORDER_MORTON);

    std::vector<Vector3>  input_vertices(sf_mesh.vertices.nb());
    std::vector<Vector3i> input_faces(sf_mesh.facets.nb());
    for (size_t i = 0; i < input_vertices.size(); i++) {
        input_vertices[i] << (sf_mesh.vertices.point(i))[0], (sf_mesh.vertices.point(i))[1],
          (sf_mesh.vertices.point(i))[2];
    }
    for (size_t i = 0; i < input_faces.size(); i++) {
        input_faces[i] << sf_mesh.facets.vertex(i, 0), sf_mesh.facets.vertex(i, 1),
          sf_mesh.facets.vertex(i, 2);
    }

    if (input_vertices.empty() || input_faces.empty()) {
        return EXIT_FAILURE;
    }

    AABBWrapper      tree(sf_mesh);
#ifdef NEW_ENVELOPE
    tree.init_sf_tree(input_vertices, input_faces, params.eps);
#endif
    std::vector<int> input_tags(input_faces.size(), 0);

    if (!params.init(tree.get_sf_diag())) {
        return EXIT_FAILURE;
    }

    stats().record(StateInfo::init_id, 0, input_vertices.size(), input_faces.size(), -1, -1);

    /////////////////////////////////////////////////
    // STEP 1: Preprocessing (mesh simplification) //
    /////////////////////////////////////////////////
    Mesh mesh;
    mesh.params = params;

    igl::Timer timer;

    timer.start();
    simplify(input_vertices, input_faces, input_tags, tree, mesh.params, skip_simplify);
    tree.init_b_mesh_and_tree(input_vertices, input_faces, mesh);
    logger().info("preprocessing {}s", timer.getElapsedTimeInSec());
    logger().info("");
    stats().record(StateInfo::preprocessing_id,
                   timer.getElapsedTimeInSec(),
                   input_vertices.size(),
                   input_faces.size(),
                   -1,
                   -1);
    if (mesh.params.log_level <= 1) {
        output_component(input_vertices, input_faces, input_tags);
    }

    ///////////////////////////////////////
    // STEP 2: Volume tetrahedralization //
    ///////////////////////////////////////

    timer.start();
    std::vector<bool> is_face_inserted(input_faces.size(), false);

    FloatTetDelaunay::tetrahedralize(input_vertices, input_faces, tree, mesh, is_face_inserted);

    logger().info("#v = {}", mesh.get_v_num());
    logger().info("#t = {}", mesh.get_t_num());
    logger().info("tetrahedralizing {}s", timer.getElapsedTimeInSec());
    logger().info("");
    stats().record(StateInfo::tetrahedralization_id,
                   timer.getElapsedTimeInSec(),
                   mesh.get_v_num(),
                   mesh.get_t_num(),
                   -1,
                   -1);

    /////////////////////
    // STEP 3: Cutting //
    /////////////////////

    timer.start();
    insert_triangles(input_vertices, input_faces, input_tags, mesh, is_face_inserted, tree, false);
    logger().info("cutting {}s", timer.getElapsedTimeInSec());
    logger().info("");
    stats().record(StateInfo::cutting_id,
                   timer.getElapsedTimeInSec(),
                   mesh.get_v_num(),
                   mesh.get_t_num(),
                   mesh.get_max_energy(),
                   mesh.get_avg_energy(),
                   std::count(is_face_inserted.begin(), is_face_inserted.end(), false));

    //////////////////////////////////////
    // STEP 4: Volume mesh optimization //
    //////////////////////////////////////

    timer.start();
    optimization(input_vertices, input_faces, input_tags, is_face_inserted, mesh, tree, {{1, 1, 1, 1}});
    logger().info("mesh optimization {}s", timer.getElapsedTimeInSec());
    logger().info("");
    stats().record(StateInfo::optimization_id,
                   timer.getElapsedTimeInSec(),
                   mesh.get_v_num(),
                   mesh.get_t_num(),
                   mesh.get_max_energy(),
                   mesh.get_avg_energy());

    /////////////////////////////////
    // STEP 5: Interior extraction //
    /////////////////////////////////

    timer.start();
    if (boolean_op < 0) {
//        filter_outside(mesh);
        if (params.smooth_open_boundary) {
            smooth_open_boundary(mesh, tree);
            for (auto &t: mesh.tets) {
                if (t.is_outside)
                    t.is_removed = true;
            }
        } else {
            if(!params.disable_filtering) {
                if(params.use_floodfill) {
                    filter_outside_floodfill(mesh);
                } else if(params.use_input_for_wn){
                    filter_outside(mesh, input_vertices, input_faces);
                } else
                    filter_outside(mesh);
            }
        }
    } else {
        boolean_operation(mesh, boolean_op);
    }
    stats().record(StateInfo::wn_id,
                   timer.getElapsedTimeInSec(),
                   mesh.get_v_num(),
                   mesh.get_t_num(),
                   mesh.get_max_energy(),
                   mesh.get_avg_energy());
    logger().info("after winding number");
    logger().info("#v = {}", mesh.get_v_num());
    logger().info("#t = {}", mesh.get_t_num());
    logger().info("winding number {}s", timer.getElapsedTimeInSec());
    logger().info("");

    MeshIO::extract_volume_mesh(mesh, V, T, false);

    if (!params.log_path.empty()) {
        std::ofstream fout(params.log_path + "_" + params.postfix + ".csv");
        if (fout) {
            fout << stats();
        }
    }

    if (!params.envelope_log.empty()) {
        std::ofstream fout(params.envelope_log);
        fout << envelope_log_csv;
    }

    return 0;
}

}  // namespace floatTetWild
