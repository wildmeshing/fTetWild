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

#include <floattetwild/TriangleInsertion.h>

#include <floattetwild/LocalOperations.h>
#include <floattetwild/Predicates.hpp>
#include <floattetwild/MeshIO.hpp>
#include <floattetwild/auto_table.hpp>
#include <floattetwild/Logger.hpp>
#include <floattetwild/intersections.h>

#include <floattetwild/MeshImprovement.h>//fortest

#include <igl/writeSTL.h>
#include <igl/writeOFF.h>
#include <igl/Timer.h>

#ifdef FLOAT_TETWILD_USE_TBB
#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include <tbb/atomic.h>
//#include <floattetwild/FloatTetCuttingParallel.h>
#include <tbb/concurrent_queue.h>
#include <tbb/concurrent_vector.h>
#endif

#include <bitset>
#include <numeric>
#include <unordered_map>


#define III -1

//fortest
#include <floattetwild/Rational.h>
double time_find_cutting_tets = 0;
double time_find_cutting_tets1 = 0;
double time_find_cutting_tets2 = 0;
double time_find_cutting_tets3 = 0;
double time_find_cutting_tets4 = 0;
double time_cut_mesh = 0;
double time_cut_mesh1 = 0;
double time_cut_mesh2 = 0;
double time_get_intersecting_edges_and_points = 0;
double time_subdivide_tets = 0;
double time_push_new_tets = 0;
double time_push_new_tets1 = 0;
double time_push_new_tets2 = 0;
double time_push_new_tets3 = 0;
double time_simplify_subdivision_result = 0;
int cnt_snapped = 0;

double old_time_find_cutting_tets = 0;
double old_time_cut_mesh = 0;
double old_time_get_intersecting_edges_and_points = 0;
double old_time_subdivide_tets = 0;
double old_time_push_new_tets = 0;
double old_time_simplify_subdivision_result = 0;

std::vector<std::array<int, 3>> covered_tet_fs;//fortest
//fortest

///two places to update quality: snapping tet vertices to plane, push new tets

floatTetWild::Vector3 floatTetWild::get_normal(const Vector3& a, const Vector3& b, const Vector3& c) {
    return ((b - c).cross(a - c)).normalized();
}

void floatTetWild::match_surface_fs(const Mesh &mesh,
                                    const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
                                    std::vector<bool> &is_face_inserted,
                                    std::vector<std::array<std::vector<int>, 4>> &track_surface_fs){
    auto comp = [](const std::array<int, 4> &a, const std::array<int, 4> &b) {
        return std::tuple<int, int, int>(a[0], a[1], a[2]) < std::tuple<int, int, int>(b[0], b[1], b[2]);
    };

    std::vector<std::array<int, 4>> input_fs(input_faces.size());
    for (int i = 0; i < input_faces.size(); i++) {
        input_fs[i] = {{input_faces[i][0], input_faces[i][1], input_faces[i][2], i}};
        std::sort(input_fs[i].begin(), input_fs[i].begin() + 3);
    }
    std::sort(input_fs.begin(), input_fs.end(), comp);

    for (int i = 0; i < mesh.tets.size(); i++) {
        auto &t = mesh.tets[i];
        for (int j = 0; j < 4; j++) {
            std::array<int, 3> f = {{t[mod4(j + 1)], t[mod4(j + 2)], t[mod4(j + 3)]}};
            std::sort(f.begin(), f.end());
            auto bounds = std::equal_range(input_fs.begin(), input_fs.end(),
                                           std::array<int, 4>({{f[0], f[1], f[2], -1}}),
                                           comp);
            for (auto it = bounds.first; it != bounds.second; ++it) {
                int f_id = (*it)[3];
                is_face_inserted[f_id] = true;
                track_surface_fs[i][j].push_back(f_id);
            }
        }
    }
}

void floatTetWild::sort_input_faces(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
                                    const Mesh &mesh, std::vector<int> &sorted_f_ids) {///use 38416, 232368 as example //todo: think why
    std::vector<Scalar> weights(input_faces.size());
    sorted_f_ids.resize(input_faces.size());
    for (int i = 0; i < input_faces.size(); i++) {
        sorted_f_ids[i] = i;

//        //fortest
//        Vector3 u = input_vertices[input_faces[i][1]] - input_vertices[input_faces[i][0]];
//        Vector3 v = input_vertices[input_faces[i][2]] - input_vertices[input_faces[i][0]];
//        if(u.cross(v).norm()/2 < SCALAR_ZERO_2) {
//            cout << "degenerate input triangle!!" << endl;
//            pausee();
//        }
//        //fortest

//        for (int j = 0; j < 3; j++) {
//            Scalar dis =
//                    (input_vertices[input_faces[i][j]] - input_vertices[input_faces[i][(j + 1) % 3]]).squaredNorm();
//            if (j == 0)
//                weights[i] = dis;
//            else if (dis > weights[i])
//                weights[i] = dis;
//        }
        Vector3 u = input_vertices[input_faces[i][1]] - input_vertices[input_faces[i][0]];
        Vector3 v = input_vertices[input_faces[i][2]] - input_vertices[input_faces[i][0]];
        weights[i] = u.cross(v).squaredNorm();
    }

    if (mesh.params.not_sort_input)
        return;

    std::random_shuffle(sorted_f_ids.begin(), sorted_f_ids.end());
//    std::sort(sorted_f_ids.begin(), sorted_f_ids.end(), [&weights](int a, int b) {
//        return weights[a] < weights[b];
//    });
}

void floatTetWild::insert_triangles(const std::vector<Vector3> &input_vertices,
                                    const std::vector<Vector3i> &input_faces, const std::vector<int> &input_tags,
                                    Mesh &mesh, std::vector<bool> &is_face_inserted, AABBWrapper &tree, bool is_again) {
    insert_triangles_aux(input_vertices, input_faces, input_tags, mesh, is_face_inserted, tree, is_again);
    return;

    if (mesh.is_input_all_inserted)
        return;

    for (auto &v: mesh.tet_vertices) {
        if (v.is_removed)
            continue;
        v.is_on_surface = false;
        v.is_on_bbox = false;
    }
    //
    for (auto &t: mesh.tets) {
        if (t.is_removed)
            continue;
        for (int j = 0; j < 4; j++) {
            if (t.is_surface_fs[j] <= 0) {
                for (int k = 0; k < 3; k++) {
                    mesh.tet_vertices[t[(j + 1 + k) % 4]].is_on_surface = true;
                    mesh.tet_vertices[t[(j + 1 + k) % 4]].is_freezed = true;
                }
            }
            if (t.is_bbox_fs[j] != NOT_BBOX) {
                mesh.tet_vertices[t[mod4(j + 1)]].is_on_bbox = true;
                mesh.tet_vertices[t[mod4(j + 2)]].is_on_bbox = true;
                mesh.tet_vertices[t[mod4(j + 3)]].is_on_bbox = true;
            }
        }
    }
    //
    operation(input_vertices, input_faces, input_tags, is_face_inserted, mesh, tree,
              std::array<int, 5>({{0, 1, 0, 1, 0}}));
    //
    for (auto &v: mesh.tet_vertices) {
        v.is_freezed = false;
    }

    insert_triangles_aux(input_vertices, input_faces, input_tags, mesh, is_face_inserted, tree, is_again);
}

void floatTetWild::optimize_non_surface(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
                                        const std::vector<int> &input_tags, std::vector<bool> &is_face_inserted,
                                        const std::vector<std::array<std::vector<int>, 4 >>& track_surface_fs,
                                        Mesh &mesh, AABBWrapper &tree, bool is_again) {
    if (!is_again) {
        for (int i = 0; i < mesh.tet_vertices.size(); i++)
            if (i < input_vertices.size())
                mesh.tet_vertices[i].is_freezed = true;
    }
    //
    for (auto &v: mesh.tet_vertices) {
        if (v.is_removed)
            continue;
        v.is_on_surface = false;
        v.is_on_bbox = false;
    }
    //
    for (int i = 0; i < mesh.tets.size(); i++) {
        auto &t = mesh.tets[i];
        if (t.is_removed)
            continue;
        for (int j = 0; j < 4; j++) {
            if (t.is_surface_fs[j] <= 0) {
                for (int k = 0; k < 3; k++) {
                    mesh.tet_vertices[t[(j + 1 + k) % 4]].is_on_surface = true;
                    mesh.tet_vertices[t[(j + 1 + k) % 4]].is_freezed = true;
                }
            }
            if (t.is_bbox_fs[j] != NOT_BBOX) {
                for (int k = 0; k < 3; k++)
                    mesh.tet_vertices[t[mod4(j + 1 + k)]].is_on_bbox = true;
            }
            if (!track_surface_fs[i][j].empty()) {
                for (int k = 0; k < 3; k++)
                    mesh.tet_vertices[t[(j + 1 + k) % 4]].is_freezed = true;
            }
        }
//        if (t.quality == 0)
//            t.quality = get_quality(mesh, t);
    }
    //
    operation(input_vertices, input_faces, input_tags, is_face_inserted, mesh, tree,
              std::array<int, 5>({{0, 1, 1, 1, 0}}));
    //
    for (auto &v: mesh.tet_vertices)
        v.is_freezed = false;
}

void floatTetWild::insert_triangles_aux(const std::vector<Vector3> &input_vertices,
        const std::vector<Vector3i> &input_faces, const std::vector<int> &input_tags,
        Mesh &mesh, std::vector<bool> &is_face_inserted,
        AABBWrapper &tree, bool is_again) {
    std::vector<bool> old_is_face_inserted = is_face_inserted;///is_face_inserted has been intialized in main

    logger().info("triangle insertion start, #f = {}, #v = {}, #t = {}",
                  input_faces.size(), mesh.tet_vertices.size(), mesh.tets.size());
    /////
    std::vector<std::array<std::vector<int>, 4 >> track_surface_fs(mesh.tets.size());
    if (!is_again) {
        match_surface_fs(mesh, input_vertices, input_faces, is_face_inserted, track_surface_fs);
    }
    int cnt_matched = std::count(is_face_inserted.begin(), is_face_inserted.end(), true);
    logger().info("matched #f = {}, uninserted #f = {}", cnt_matched, is_face_inserted.size() - cnt_matched);

    std::vector<int> sorted_f_ids;
    sort_input_faces(input_vertices, input_faces, mesh, sorted_f_ids);

    /////
    std::vector<Vector3> new_vertices;
    std::vector<std::array<int, 4>> new_tets;
    int cnt_fail = 0;
    int cnt_total = 0;

    /////
//    if(!is_again) {
//        std::vector<std::vector<int>> conn_fs(input_vertices.size());
//        for (int i = 0; i < input_faces.size(); i++) {
//            for (int j = 0; j < 3; j++)
//                conn_fs[input_faces[i][j]].push_back(i);
//        }
//        std::vector<bool> is_visited(input_faces.size(), false);
//        std::vector<Vector3> ns(input_faces.size());
//        for (int i = 0; i < input_faces.size(); i++) {
//            ns[i] = (input_vertices[input_faces[i][1]] - input_vertices[input_faces[i][0]]).cross(
//                    input_vertices[input_faces[i][2]] - input_vertices[input_faces[i][0]]).normalized();
//        }
//        for (int i = 0; i < sorted_f_ids.size(); i++) {
//            int f_id = sorted_f_ids[i];
//            if (is_face_inserted[f_id])
//                continue;
//
//            std::vector<int> f_ids;
//            if (insert_multi_triangles(f_id, input_vertices, input_faces, input_tags,
//                                       conn_fs, ns, is_visited, f_ids,
//                                       mesh, track_surface_fs, tree, is_again)) {
//                for (int inserted_f_id: f_ids)
//                    is_face_inserted[inserted_f_id] = true;
//            } else
//                cnt_total += f_ids.size();
//            cnt_total += f_ids.size();
//        }
//        logger().info("insert_multi_triangles * n done, #v = {}, #t = {}", mesh.tet_vertices.size(), mesh.tets.size());
//        logger().info("uninserted #f = {}/{}", std::count(is_face_inserted.begin(), is_face_inserted.end(), false),
//                      is_face_inserted.size() - cnt_matched);
//    }

    //////
    for (int i = 0; i < sorted_f_ids.size(); i++) {
        //fortest
        if (!is_again && i > 0 && i % 1000 == 0) {
            logger().debug("inserting f{}... {} failed", i, cnt_fail);
            logger().debug("snapped {}/{}", cnt_snapped, cnt_total);
            logger().debug("\t- time_find_cutting_tets = {}s (total {}s)",
                          time_find_cutting_tets - old_time_find_cutting_tets, time_find_cutting_tets);
//            logger().info("\t\t- time_find_cutting_tets1 = {}s", time_find_cutting_tets1);
//            logger().info("\t\t- time_find_cutting_tets2 = {}s", time_find_cutting_tets2);
//            logger().info("\t\t- time_find_cutting_tets3 = {}s", time_find_cutting_tets3);
//            logger().info("\t\t- time_find_cutting_tets4 = {}s", time_find_cutting_tets4);
            logger().debug("\t- time_cut_mesh = {}s (total {}s)",
                          time_cut_mesh - old_time_cut_mesh, time_cut_mesh);
//            logger().info("\t\t- time_cut_mesh1 = {}s", time_cut_mesh1);
//            logger().info("\t\t- time_cut_mesh2 = {}s", time_cut_mesh2);
//            print_times1();
            logger().debug("\t- time_get_intersecting_edges_and_points = {}s (total {}s)",
                          time_get_intersecting_edges_and_points - old_time_get_intersecting_edges_and_points,
                          time_get_intersecting_edges_and_points);
            print_times1();
            logger().debug("\t- time_subdivide_tets = {}s (total {}s)",
                          time_subdivide_tets - old_time_subdivide_tets, time_subdivide_tets);
            logger().debug("\t- time_push_new_tets = {}s (total {}s)",
                          time_push_new_tets - old_time_push_new_tets, time_push_new_tets);
//            logger().info("\t\t- time_push_new_tets1 = {}s", time_push_new_tets1);
//            logger().info("\t\t- time_push_new_tets2 = {}s", time_push_new_tets2);
//            logger().info("\t\t- time_push_new_tets3 = {}s", time_push_new_tets3);
            logger().debug("\t- time_simplify_subdivision_result = {}s (total {}s)",
                          time_simplify_subdivision_result - old_time_simplify_subdivision_result,
                          time_simplify_subdivision_result);

            old_time_find_cutting_tets = time_find_cutting_tets;
            old_time_cut_mesh = time_cut_mesh;
            old_time_get_intersecting_edges_and_points = time_get_intersecting_edges_and_points;
            old_time_subdivide_tets = time_subdivide_tets;
            old_time_push_new_tets = time_push_new_tets;
            old_time_simplify_subdivision_result = time_simplify_subdivision_result;
            logger().debug("#v = {}/{}", mesh.get_v_num(), mesh.tet_vertices.size());
            logger().debug("#t = {}/{}", mesh.get_t_num(), mesh.tets.size());
        }
        //fortest

//        //fortest
//        if(i>0 && i%10000 == 0) {
//            logger().info("before opt");
//            logger().info("#v = {}/{}", mesh.get_v_num(), mesh.tet_vertices.size());
//            logger().info("#t = {}/{}", mesh.get_t_num(), mesh.tets.size());
//            optimize_non_surface(input_vertices, input_faces, input_tags, is_face_inserted, track_surface_fs,
//                                 mesh, tree, is_again);
//            logger().info("after opt");
//            logger().info("#v = {}/{}", mesh.get_v_num(), mesh.tet_vertices.size());
//            logger().info("#t = {}/{}", mesh.get_t_num(), mesh.tets.size());
//        }
//        //fortest

        int f_id = sorted_f_ids[i];
        if (is_face_inserted[f_id])
            continue;

        cnt_total++;
        if (insert_one_triangle(f_id, input_vertices, input_faces, input_tags, mesh, track_surface_fs,
                                tree, is_again))
            is_face_inserted[f_id] = true;
        else
            cnt_fail++;

//        pausee();//fortest
        if (f_id == III)
            break;//fortest
    }
    logger().info("insert_one_triangle * n done, #v = {}, #t = {}", mesh.tet_vertices.size(), mesh.tets.size());
    logger().info("uninserted #f = {}/{}", std::count(is_face_inserted.begin(), is_face_inserted.end(), false),
                  is_face_inserted.size() - cnt_matched);
    logger().info("total timing: {}s", time_find_cutting_tets + time_cut_mesh + time_get_intersecting_edges_and_points +
                                       time_subdivide_tets + time_push_new_tets);

    pair_track_surface_fs(mesh, track_surface_fs);
    logger().info("pair_track_surface_fs done");


    /////
    std::vector<std::array<int, 2>> b_edges1;
    std::vector<std::pair<std::array<int, 2>, std::vector<int>>> b_edge_infos;
    std::vector<bool> is_on_cut_edges;
    find_boundary_edges(input_vertices, input_faces, is_face_inserted, old_is_face_inserted,
            b_edge_infos, is_on_cut_edges, b_edges1);
    logger().info("find_boundary_edges done");
    std::vector<std::array<int, 3>> known_surface_fs;
    std::vector<std::array<int, 3>> known_not_surface_fs;
    insert_boundary_edges(input_vertices, input_faces, b_edge_infos, is_on_cut_edges, track_surface_fs, mesh, tree,
                          is_face_inserted, is_again, known_surface_fs, known_not_surface_fs);
    logger().info("uninserted #f = {}/{}", std::count(is_face_inserted.begin(), is_face_inserted.end(), false),
                  is_face_inserted.size() - cnt_matched);

    //fortest
    check_track_surface_fs(mesh, track_surface_fs, input_vertices, input_faces, sorted_f_ids);
    //fortest

    /////
    std::vector<std::array<int, 2>> b_edges2;
    mark_surface_fs(input_vertices, input_faces, input_tags, track_surface_fs, is_face_inserted,
                    known_surface_fs, known_not_surface_fs, b_edges2, mesh, tree);
    //fortest: output and check
//    output_surface(mesh, mesh.params.output_path+"_"+mesh.params.postfix+"_surface.stl");
    logger().info("mark_surface_fs done");

    /////
    //build b_tree using b_edges
    tree.init_tmp_b_mesh_and_tree(input_vertices, input_faces, b_edges1, mesh, b_edges2);
//    for (int v_id = 0; v_id < mesh.tet_vertices.size(); v_id++) {
//        if (mesh.tet_vertices[v_id].is_removed)
//            continue;
//        if (!mesh.tet_vertices[v_id].is_on_boundary)
//            continue;
//
//        GEO::index_t prev_facet;
//        if (tree.is_out_tmp_b_envelope(mesh.tet_vertices[v_id].pos, mesh.params.eps_2, prev_facet))
//            mesh.tet_vertices[v_id].is_on_boundary = false;
//    }

#ifdef FLOAT_TETWILD_USE_TBB
    tbb::parallel_for(size_t(0), mesh.tets.size(), [&](size_t i){
        auto &t = mesh.tets[i];
#else
    for (auto &t:mesh.tets) {
#endif
        if (!t.is_removed) {
            t.quality = get_quality(mesh, t);
        }
#ifdef FLOAT_TETWILD_USE_TBB
    });
#else
    }
#endif

    if (std::count(is_face_inserted.begin(), is_face_inserted.end(), false) == 0)
        mesh.is_input_all_inserted = true;
    logger().info("#b_edge1 = {}, #b_edges2 = {}", b_edges1.size(), b_edges2.size());

//    ///fortest
//    Eigen::MatrixXd V(input_vertices.size(), 3);
//    Eigen::MatrixXi F(std::count(is_face_inserted.begin(), is_face_inserted.end(), false), 3);
//    for (int i = 0; i < input_vertices.size(); i++)
//        V.row(i) = input_vertices[i];
//    int cnt = 0;
//    for (int i = 0; i < input_faces.size(); i++) {
//        if (is_face_inserted[i])
//            continue;
//        F.row(cnt) << input_faces[i][0], input_faces[i][1], input_faces[i][2];
//        cnt++;
//    }
//    igl::writeSTL(mesh.params.output_path+"_"+mesh.params.postfix+"_uninserted.stl", V, F);
//    //
//    std::ofstream fout(mesh.params.output_path+"_"+mesh.params.postfix+"_b_es.obj");
//    for(int i=0;i<tree.tmp_b_mesh.vertices.nb();i++){
//        fout<<"v "<<tree.tmp_b_mesh.vertices.point(i)[0]<<" "
//                <<tree.tmp_b_mesh.vertices.point(i)[1]<<" "
//                <<tree.tmp_b_mesh.vertices.point(i)[2]<<endl;
//    }
//    for(int i=0;i<tree.tmp_b_mesh.facets.nb();i++) {
//        fout << "l " << tree.tmp_b_mesh.facets.vertex(i, 1) + 1 << " "
//             << tree.tmp_b_mesh.facets.vertex(i, 2) + 1 << endl;
//    }
//    fout.close();
//    //fortest

    pausee();
}

bool floatTetWild::insert_multi_triangles(int insert_f_id, const std::vector<Vector3> &input_vertices,
                            const std::vector<Vector3i> &input_faces, const std::vector<int> &input_tags,
                            const std::vector<std::vector<int>>& conn_fs,
                            const std::vector<Vector3> &ns, std::vector<bool> & is_visited, std::vector<int>& f_ids,
                            Mesh &mesh, std::vector<std::array<std::vector<int>, 4>> &track_surface_fs,
                            AABBWrapper &tree, bool is_again){
    std::array<Vector3, 3> vs = {{input_vertices[input_faces[insert_f_id][0]],
                                         input_vertices[input_faces[insert_f_id][1]],
                                         input_vertices[input_faces[insert_f_id][2]]}};
    const auto& n = ns[insert_f_id];
    int t = get_t(vs[0], vs[1], vs[2]);

    f_ids.push_back(insert_f_id);
    std::queue<int> f_queue;
    f_queue.push(insert_f_id);
    while(!f_queue.empty()){
        int f_id = f_queue.front();
        f_queue.pop();
        for(int j=0;j<3;j++){
            for(int n_f_id: conn_fs[input_faces[f_id][j]]) {
                if (is_visited[n_f_id] || abs(ns[n_f_id].dot(n) - 1) > 1e-6)
                    continue;
                is_visited[n_f_id] = true;
                f_queue.push(n_f_id);
                f_ids.push_back(n_f_id);
            }
        }
    }
    vector_unique(f_ids);
    if(f_ids.size()<2) {
        return false;
    }
//    cout<<"f_ids.size() = "<<f_ids.size()<<endl;

    std::vector<int> cut_t_ids;
    for(int f_id: f_ids) {
        std::array<Vector3, 3> tmp_vs = {{input_vertices[input_faces[f_id][0]],
                                             input_vertices[input_faces[f_id][1]],
                                             input_vertices[input_faces[f_id][2]]}};
        std::vector<int> tmp_cut_t_ids;
        find_cutting_tets(f_id, input_vertices, input_faces, tmp_vs, mesh, tmp_cut_t_ids, is_again);

        CutMesh cut_mesh(mesh, n, tmp_vs);
        cut_mesh.construct(tmp_cut_t_ids);
        if (cut_mesh.snap_to_plane()) {
            cnt_snapped++;
            cut_mesh.project_to_plane(input_vertices.size());
            cut_mesh.expand_new(tmp_cut_t_ids);
            cut_mesh.project_to_plane(input_vertices.size());
        }
        cut_t_ids.insert(cut_t_ids.end(), tmp_cut_t_ids.begin(), tmp_cut_t_ids.end());
    }
    vector_unique(cut_t_ids);
    if (cut_t_ids.empty()) {
//        cout<<"fail 0"<<endl;
        return false;
    }

    CutMesh cut_mesh(mesh, n, vs);
    cut_mesh.construct(cut_t_ids);
    cut_mesh.snap_to_plane();
//    if (cut_mesh.snap_to_plane()) {
//        cnt_snapped++;
//        cut_mesh.project_to_plane(input_vertices.size());
////        cut_mesh.expand_new(cut_t_ids);
////        cut_mesh.project_to_plane(input_vertices.size());
//    }

    std::vector<Vector3> points;
    std::map<std::array<int, 2>, int> map_edge_to_intersecting_point;
    std::vector<int> subdivide_t_ids;
    if (!cut_mesh.get_intersecting_edges_and_points(points, map_edge_to_intersecting_point, subdivide_t_ids)) {
//        cout<<"fail 1"<<endl;
        return false;
    }

    vector_unique(cut_t_ids);
    std::vector<int> tmp;
    std::set_difference(subdivide_t_ids.begin(), subdivide_t_ids.end(), cut_t_ids.begin(), cut_t_ids.end(),
                        std::back_inserter(tmp));
    std::vector<bool> is_mark_surface(cut_t_ids.size(), true);
    cut_t_ids.insert(cut_t_ids.end(), tmp.begin(), tmp.end());
    is_mark_surface.resize(is_mark_surface.size() + tmp.size(), false);

    std::vector<MeshTet> new_tets;
    std::vector<std::array<std::vector<int>, 4>> new_track_surface_fs;
    std::vector<int> modified_t_ids;
    if (!subdivide_tets(insert_f_id, mesh, cut_mesh, points, map_edge_to_intersecting_point, track_surface_fs,
                        cut_t_ids, is_mark_surface,
                        new_tets, new_track_surface_fs, modified_t_ids)) {
//        cout<<"fail 2"<<endl;
        return false;
    }
    push_new_tets(mesh, track_surface_fs, points, new_tets, new_track_surface_fs, modified_t_ids, is_again);
    simplify_subdivision_result(insert_f_id, input_vertices.size(), mesh, tree, track_surface_fs);

    return true;
}

bool floatTetWild::insert_one_triangle(int insert_f_id, const std::vector<Vector3> &input_vertices,
        const std::vector<Vector3i> &input_faces, const std::vector<int> &input_tags,
        Mesh &mesh, std::vector<std::array<std::vector<int>, 4>>& track_surface_fs,
        AABBWrapper &tree, bool is_again) {
//    igl::Timer timer;
    std::array<Vector3, 3> vs = {{input_vertices[input_faces[insert_f_id][0]],
                                         input_vertices[input_faces[insert_f_id][1]],
                                         input_vertices[input_faces[insert_f_id][2]]}};
    Vector3 n = (vs[1] - vs[0]).cross(vs[2] - vs[0]);
    n.normalize();
    int t = get_t(vs[0], vs[1], vs[2]);

    /////
//    timer.start();
    std::vector<int> cut_t_ids;
    find_cutting_tets(insert_f_id, input_vertices, input_faces, vs, mesh, cut_t_ids, is_again);
//    time_find_cutting_tets += timer.getElapsedTime();

    //fortest
    myassert(!cut_t_ids.empty(), "cut_t_ids.empty()!!!");
    if (cut_t_ids.empty()) {
        cout << get_area(vs[0], vs[1], vs[2]) << endl;
        cout << "f" << insert_f_id << ": " << input_faces[insert_f_id][0] << " " << input_faces[insert_f_id][1]
             << " " << input_faces[insert_f_id][2] << endl;
        pausee();
        return false;
    }
    //fortest

    /////
//    timer.start();
//    igl::Timer timer1;
//    timer1.start();
    CutMesh cut_mesh(mesh, n, vs);
    cut_mesh.construct(cut_t_ids);
//    time_cut_mesh1 += timer1.getElapsedTime();
//    timer1.start();
//    bool is_expanded = false;//fortest

    if (cut_mesh.snap_to_plane()) {
        cnt_snapped++;
        cut_mesh.project_to_plane(input_vertices.size());
        cut_mesh.expand_new(cut_t_ids);
        cut_mesh.project_to_plane(input_vertices.size());
        //fortest
//        int cnt_proj = cut_mesh.project_to_plane(input_vertices.size());
//        int cnt_all = std::count(cut_mesh.is_snapped.begin(), cut_mesh.is_snapped.end(), true);
//        if (cnt_proj != cnt_all)
//            cout << cnt_proj << "/" << cnt_all << endl;
        //fortest
    }
//    time_cut_mesh2 += timer1.getElapsedTime();
//    time_cut_mesh += timer.getElapsedTime();

    /////
//    timer.start();
    std::vector<Vector3> points;
    std::map<std::array<int, 2>, int> map_edge_to_intersecting_point;
    std::vector<int> subdivide_t_ids;
    if (!cut_mesh.get_intersecting_edges_and_points(points, map_edge_to_intersecting_point, subdivide_t_ids)) {
//        time_get_intersecting_edges_and_points += timer.getElapsedTime();
        if(is_again){
            if(is_uninserted_face_covered(insert_f_id, input_vertices, input_faces, cut_t_ids, mesh))
                return true;
        }
        cout<<"FAIL get_intersecting_edges_and_points"<<endl;
        return false;
    }
//    time_get_intersecting_edges_and_points += timer.getElapsedTime();
    //have to add all cut_t_ids
    vector_unique(cut_t_ids);
    std::vector<int> tmp;
    std::set_difference(subdivide_t_ids.begin(), subdivide_t_ids.end(), cut_t_ids.begin(), cut_t_ids.end(),
                        std::back_inserter(tmp));
    std::vector<bool> is_mark_surface(cut_t_ids.size(), true);
    cut_t_ids.insert(cut_t_ids.end(), tmp.begin(), tmp.end());
    is_mark_surface.resize(is_mark_surface.size() + tmp.size(), false);
//    cout << "cut_mesh.get_intersecting_edges_and_points OK" << endl;
//    time_get_intersecting_edges_and_points += timer.getElapsedTime();

    /////
//    timer.start();
    std::vector<MeshTet> new_tets;
    std::vector<std::array<std::vector<int>, 4>> new_track_surface_fs;
    std::vector<int> modified_t_ids;
    if (!subdivide_tets(insert_f_id, mesh, cut_mesh, points, map_edge_to_intersecting_point, track_surface_fs,
                        cut_t_ids, is_mark_surface,
                        new_tets, new_track_surface_fs, modified_t_ids)) {
//        time_subdivide_tets += timer.getElapsedTime();
        if(is_again){
            if(is_uninserted_face_covered(insert_f_id, input_vertices, input_faces, cut_t_ids, mesh))
                return true;
        }
        cout<<"FAIL subdivide_tets"<<endl;
        return false;
    }
//    time_subdivide_tets += timer.getElapsedTime();

//    timer.start();
    push_new_tets(mesh, track_surface_fs, points, new_tets, new_track_surface_fs, modified_t_ids, is_again);
//    time_push_new_tets += timer.getElapsedTime();

//    timer.start();
    simplify_subdivision_result(insert_f_id, input_vertices.size(), mesh, tree, track_surface_fs);
//    time_simplify_subdivision_result += timer.getElapsedTime();

    return true;
}

void floatTetWild::push_new_tets(Mesh &mesh, std::vector<std::array<std::vector<int>, 4>> &track_surface_fs,
                                 std::vector<Vector3> &points, std::vector<MeshTet> &new_tets,
                                 std::vector<std::array<std::vector<int>, 4>> &new_track_surface_fs,
                                 std::vector<int> &modified_t_ids, bool is_again) {
//    igl::Timer timer;
//    timer.start();
    ///vs
    const int old_v_size = mesh.tet_vertices.size();
    mesh.tet_vertices.resize(mesh.tet_vertices.size() + points.size());
    for (int i = 0; i < points.size(); i++) {
        mesh.tet_vertices[old_v_size + i].pos = points[i];
        //todo: tags???
    }
//    time_push_new_tets1 += timer.getElapsedTime();

    ///tets
//    timer.start();
//    mesh.tets.reserve(mesh.tets.size() + new_tets.size() - modified_t_ids.size());
//    time_push_new_tets2 += timer.getElapsedTime();

//    timer.start();
    for (int i = 0; i < new_tets.size(); i++) {
//        if (is_again)
//            new_tets[i].quality = get_quality(mesh, new_tets[i]);

        if (i < modified_t_ids.size()) {
            for (int j = 0; j < 4; j++) {
                vector_erase(mesh.tet_vertices[mesh.tets[modified_t_ids[i]][j]].conn_tets, modified_t_ids[i]);
            }
            mesh.tets[modified_t_ids[i]] = new_tets[i];
            track_surface_fs[modified_t_ids[i]] = new_track_surface_fs[i];
            for (int j = 0; j < 4; j++) {
                mesh.tet_vertices[mesh.tets[modified_t_ids[i]][j]].conn_tets.push_back(modified_t_ids[i]);
            }
        } else {
//            mesh.tets.push_back(new_tets[i]);
//            track_surface_fs.push_back(new_track_surface_fs[i]);
//            for (int j = 0; j < 4; j++) {
//                mesh.tet_vertices[mesh.tets.back()[j]].conn_tets.push_back(mesh.tets.size() - 1);
//            }
            for (int j = 0; j < 4; j++) {
                mesh.tet_vertices[new_tets[i][j]].conn_tets.push_back(mesh.tets.size() + i - modified_t_ids.size());
            }
        }
        //todo: tags???
    }
//    time_push_new_tets2 += timer.getElapsedTime();

//    timer.start();
    mesh.tets.insert(mesh.tets.end(), new_tets.begin() + modified_t_ids.size(), new_tets.end());
    track_surface_fs.insert(track_surface_fs.end(), new_track_surface_fs.begin() + modified_t_ids.size(),
                            new_track_surface_fs.end());
    modified_t_ids.clear();
//    time_push_new_tets3 += timer.getElapsedTime();
}

#include <floattetwild/EdgeCollapsing.h>
void floatTetWild::simplify_subdivision_result(int insert_f_id, int input_v_size, Mesh &mesh, AABBWrapper &tree,
        std::vector<std::array<std::vector<int>, 4>> &track_surface_fs) {
    if(covered_tet_fs.empty())
        return;

    for (int i = 0; i < covered_tet_fs.size(); i++)
        std::sort(covered_tet_fs[i].begin(), covered_tet_fs[i].end());
    vector_unique(covered_tet_fs);

    std::vector<std::array<int, 3>> edges;
    for (int i = 0; i < covered_tet_fs.size(); i++) {
        const auto f = covered_tet_fs[i];
        for (int j = 0; j < 3; j++) {
            if (f[j] < f[(j + 1) % 3])
                edges.push_back({{f[j], f[(j + 1) % 3], i}});
            else
                edges.push_back({{f[(j + 1) % 3], f[j], i}});
        }
    }
    std::sort(edges.begin(), edges.end(), [](const std::array<int, 3> &a, const std::array<int, 3> &b) {
        return std::make_tuple(a[0], a[1]) < std::make_tuple(b[0], b[1]);
    });
    //
    std::unordered_set<int> freezed_v_ids;
    bool is_duplicated = false;
    for (int i = 0; i < edges.size() - 1; i++) {
        if(edges[i][0]<input_v_size)
            freezed_v_ids.insert(edges[i][0]);
        if(edges[i][1]<input_v_size)
            freezed_v_ids.insert(edges[i][1]);

        if (edges[i][0] == edges[i + 1][0] && edges[i][1] == edges[i + 1][1]) {
            is_duplicated = true;
            edges.erase(edges.begin() + i);
            i--;
        } else {
            if (!is_duplicated) {///boundary edges
                freezed_v_ids.insert(edges[i][0]);
                freezed_v_ids.insert(edges[i][1]);
            }
            is_duplicated = false;
        }
    }
    if (!is_duplicated) {
        freezed_v_ids.insert(edges.back()[0]);
        freezed_v_ids.insert(edges.back()[1]);
    }

    std::priority_queue<ElementInQueue, std::vector<ElementInQueue>, cmp_s> ec_queue;
    for (const auto &e:edges) {
        Scalar l_2 = get_edge_length_2(mesh, e[0], e[1]);
        if (freezed_v_ids.find(e[0]) == freezed_v_ids.end())
            ec_queue.push(ElementInQueue({{e[0], e[1]}}, l_2));
        if (freezed_v_ids.find(e[1]) == freezed_v_ids.end())
            ec_queue.push(ElementInQueue({{e[1], e[0]}}, l_2));
    }
    //
    if(ec_queue.empty())
        return;
    //
    std::unordered_set<int> all_v_ids;
    for (const auto &e:edges) {
        all_v_ids.insert(e[0]);
        all_v_ids.insert(e[1]);
    }
    //
    int _ts = 0;
    std::vector<int> _tet_tss;
    bool is_update_tss = false;
    int cnt_suc = 0;
    while (!ec_queue.empty()) {
        std::array<int, 2> v_ids = ec_queue.top().v_ids;
        Scalar old_weight = ec_queue.top().weight;
        ec_queue.pop();

        while (!ec_queue.empty()) {
            if (ec_queue.top().v_ids == v_ids)
                ec_queue.pop();
            else
                break;
        }

        if (!is_valid_edge(mesh, v_ids[0], v_ids[1]))
            continue;
        if (freezed_v_ids.find(v_ids[0]) != freezed_v_ids.end())
            continue;
        Scalar weight = get_edge_length_2(mesh, v_ids[0], v_ids[1]);
        if (weight != old_weight)
            continue;

        //check track_surface_fs
        int v1_id = v_ids[0];
        int v2_id = v_ids[1];
        bool is_valid = true;
        for (int t_id: mesh.tet_vertices[v1_id].conn_tets) {
            for (int j = 0; j < 4; j++) {
                if ((!track_surface_fs[t_id][j].empty() && !vector_contains(track_surface_fs[t_id][j], insert_f_id))
                    || (mesh.tets[t_id][j] != v1_id && (mesh.tets[t_id].is_surface_fs[j] != NOT_SURFACE || mesh.tets[t_id].is_bbox_fs[j] != NOT_BBOX))) {
                    is_valid = false;
                    break;
                }
            }
            if (!is_valid)
                break;
        }
        if (!is_valid)
            continue;

        std::vector<std::array<int, 2>> new_edges;
        static const bool is_check_quality = true;
        auto v1_conn_tets = mesh.tet_vertices[v1_id].conn_tets;
        for (int t_id: mesh.tet_vertices[v1_id].conn_tets) {
//            if(mesh.tets[t_id].quality == 0)
                mesh.tets[t_id].quality = get_quality(mesh, t_id);
        }
        int result = collapse_an_edge(mesh, v_ids[0], v_ids[1], tree, new_edges, _ts, _tet_tss,
                                      is_check_quality, is_update_tss);
        if (result > 0) {
            for(const auto& e: new_edges){
                if(all_v_ids.find(e[0]) == all_v_ids.end() || all_v_ids.find(e[1]) == all_v_ids.end())
                    continue;
                Scalar l_2 = get_edge_length_2(mesh, e[0], e[1]);
                if (freezed_v_ids.find(e[0]) == freezed_v_ids.end())
                    ec_queue.push(ElementInQueue({{e[0], e[1]}}, l_2));
                if (freezed_v_ids.find(e[1]) == freezed_v_ids.end())
                    ec_queue.push(ElementInQueue({{e[1], e[0]}}, l_2));
            }
            cnt_suc++;
        }

//        ////
//        std::vector<int> n12_t_ids;
//        set_intersection(mesh.tet_vertices[v1_id].conn_tets, mesh.tet_vertices[v2_id].conn_tets, n12_t_ids);
//        std::vector<int> n1_t_ids;//v1.conn_tets - n12_t_ids
//        std::sort(mesh.tet_vertices[v1_id].conn_tets.begin(), mesh.tet_vertices[v1_id].conn_tets.end());
//        std::sort(n12_t_ids.begin(), n12_t_ids.end());
//        std::set_difference(mesh.tet_vertices[v1_id].conn_tets.begin(), mesh.tet_vertices[v1_id].conn_tets.end(),
//                            n12_t_ids.begin(), n12_t_ids.end(), std::back_inserter(n1_t_ids));
//
//        //inversion
//        std::vector<int> js_n1_t_ids;
//        for (int t_id:n1_t_ids) {
//            int j = mesh.tets[t_id].find(v1_id);
//            js_n1_t_ids.push_back(j);
//            if (is_inverted(mesh, t_id, j, mesh.tet_vertices[v2_id].pos)) {
//                is_valid = false;
//                break;
//            }
//        }
//        if (!is_valid)
//            continue;
//
//        //check quality
//        //todo
//
//        //real update
//        mesh.tet_vertices[v2_id].pos = mesh.tet_vertices[v1_id].pos;
//        int ii = 0;
//        for (int t_id:n1_t_ids) {
//            int j = js_n1_t_ids[ii++];
//            mesh.tets[t_id][j] = v2_id;
//            mesh.tet_vertices[v2_id].conn_tets.push_back(t_id);
//        }
//        for (int t_id: n12_t_ids) {
//            mesh.tets[t_id].is_removed = true;
//            for (int j = 0; j < 4; j++) {
//                if (mesh.tets[t_id][j] != v1_id)
//                    vector_erase(mesh.tet_vertices[mesh.tets[t_id][j]].conn_tets, t_id);
//            }
//        }
//        mesh.tet_vertices[v1_id].is_removed = true;
//        mesh.tet_vertices[v1_id].conn_tets.clear();
//
//        cnt_suc++;
    }
//    if (cnt_suc > 0)
//        cout << "cnt_suc = " << cnt_suc << endl;
}

void floatTetWild::find_cutting_tets(int f_id, const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
                                     const std::array<Vector3, 3>& vs, Mesh &mesh, std::vector<int> &cut_t_ids, bool is_again) {
    std::vector<bool> is_visited(mesh.tets.size(), false);
    std::queue<int> queue_t_ids;

    if (!is_again) {
        std::vector<int> n_t_ids;
        for (int j = 0; j < 3; j++) {
            n_t_ids.insert(n_t_ids.end(), mesh.tet_vertices[input_faces[f_id][j]].conn_tets.begin(),
                           mesh.tet_vertices[input_faces[f_id][j]].conn_tets.end());
        }
        vector_unique(n_t_ids);

        for (int t_id: n_t_ids) {
            is_visited[t_id] = true;
            queue_t_ids.push(t_id);
        }
    } else {
        Vector3 min_f, max_f;
        get_bbox_face(input_vertices[input_faces[f_id][0]], input_vertices[input_faces[f_id][1]],
                      input_vertices[input_faces[f_id][2]], min_f, max_f);
#ifdef FLOAT_TETWILD_USE_TBB
        tbb::concurrent_vector<int> tbb_t_ids;
        tbb::parallel_for(size_t(0), mesh.tets.size(), [&](size_t t_id){
            if (mesh.tets[t_id].is_removed)
                return;

            Vector3 min_t, max_t;
            get_bbox_tet(mesh.tet_vertices[mesh.tets[t_id][0]].pos, mesh.tet_vertices[mesh.tets[t_id][1]].pos,
                         mesh.tet_vertices[mesh.tets[t_id][2]].pos, mesh.tet_vertices[mesh.tets[t_id][3]].pos,
                         min_t, max_t);
            if (!is_bbox_intersected(min_f, max_f, min_t, max_t))
                return;

            tbb_t_ids.push_back(t_id);

        });
        for(int t_id: tbb_t_ids) {
            queue_t_ids.push(t_id);
            is_visited[t_id] = true;
        }
#else
        for (int t_id = 0; t_id < mesh.tets.size(); t_id++) {
            if (mesh.tets[t_id].is_removed)
                continue;

            Vector3 min_t, max_t;
            get_bbox_tet(mesh.tet_vertices[mesh.tets[t_id][0]].pos, mesh.tet_vertices[mesh.tets[t_id][1]].pos,
                         mesh.tet_vertices[mesh.tets[t_id][2]].pos, mesh.tet_vertices[mesh.tets[t_id][3]].pos,
                         min_t, max_t);
            if (!is_bbox_intersected(min_f, max_f, min_t, max_t))
                continue;

            queue_t_ids.push(t_id);
            is_visited[t_id] = true;
        }
#endif
    }

//    const int CUT_UNKNOWN = INT_MIN;
//    std::vector<std::array<int, 4>> visited_results(mesh.tets.size(), {{CUT_UNKNOWN, CUT_UNKNOWN, CUT_UNKNOWN, CUT_UNKNOWN}});
//    const int test_f_id = 771;
//    const int test_f_id = -1;
//    const int test_j = 1;
//    const int test_v_id = input_faces[test_f_id][test_j];
    while (!queue_t_ids.empty()) {
        int t_id = queue_t_ids.front();
        queue_t_ids.pop();

        if (is_again) {
            bool is_cut = false;
            for (int j = 0; j < 3; j++) {
                if (is_point_inside_tet(input_vertices[input_faces[f_id][j]],
                                        mesh.tet_vertices[mesh.tets[t_id][0]].pos,
                                        mesh.tet_vertices[mesh.tets[t_id][1]].pos,
                                        mesh.tet_vertices[mesh.tets[t_id][2]].pos,
                                        mesh.tet_vertices[mesh.tets[t_id][3]].pos)) {
                    is_cut = true;
                    break;
                }
            }
            if(is_cut) {
                cut_t_ids.push_back(t_id);
                for (int j = 0; j < 4; j++) {
                    for (int n_t_id: mesh.tet_vertices[mesh.tets[t_id][j]].conn_tets) {
                        if (!is_visited[n_t_id]) {
                            is_visited[n_t_id] = true;
                            queue_t_ids.push(n_t_id);
                        }
                    }
                }
                continue;
            }
        }


        std::array<int, 4> oris;
//        int cnt_pos = 0;
//        int cnt_neg = 0;
//        int cnt_on = 0;
        for (int j = 0; j < 4; j++) {
            oris[j] = Predicates::orient_3d(vs[0], vs[1], vs[2], mesh.tet_vertices[mesh.tets[t_id][j]].pos);
//            if (oris[j] == Predicates::ORI_ZERO)
//                cnt_on++;
//            else if (oris[j] == Predicates::ORI_POSITIVE)
//                cnt_pos++;
//            else
//                cnt_neg++;
        }
//        if (cnt_pos == 0 && cnt_neg == 0 && cnt_on < 3)
//            continue;

        bool is_cutted = false;
        std::vector<bool> is_cut_vs = {{false, false, false, false}}; /// is v on cut face
        for (int j = 0; j < 4; j++) {
            int cnt_pos = 0;
            int cnt_neg = 0;
            int cnt_on = 0;
            for (int k = 0; k < 3; k++) {
                if (oris[(j + k + 1) % 4] == Predicates::ORI_ZERO)
                    cnt_on++;
                else if (oris[(j + k + 1) % 4] == Predicates::ORI_POSITIVE)
                    cnt_pos++;
                else
                    cnt_neg++;
            }

            int result = CUT_EMPTY;
            auto &tp1 = mesh.tet_vertices[mesh.tets[t_id][(j + 1) % 4]].pos;
            auto &tp2 = mesh.tet_vertices[mesh.tets[t_id][(j + 2) % 4]].pos;
            auto &tp3 = mesh.tet_vertices[mesh.tets[t_id][(j + 3) % 4]].pos;

//            if (cnt_on == 3) {
//                result = is_tri_tri_cutted_hint(vs[0], vs[1], vs[2], tp1, tp2, tp3, CUT_COPLANAR);
//            } else if (cnt_pos > 0 && cnt_neg > 0) {
//                result = is_tri_tri_cutted_hint(vs[0], vs[1], vs[2], tp1, tp2, tp3, CUT_FACE);
//            }
//            if (result == CUT_EMPTY)
//                continue;
//
//            is_cutted = true;
//            is_cut_vs[(j + 1) % 4] = true;
//            is_cut_vs[(j + 2) % 4] = true;
//            is_cut_vs[(j + 3) % 4] = true;

//            //fortest
//            if(t_id == 1016 && f_id == test_f_id) {
//                cout<<"1016"<<endl;
//                mesh.tets[t_id].print();
//                cout<<oris[0]<<" "<<oris[1]<<" "<<oris[2]<<" "<<oris[3]<<endl;
//                is_tri_tri_cutted_hint(vs[0], vs[1], vs[2], tp1, tp2, tp3, CUT_FACE, true);
//                pausee();
//            }
//            //fortest

            if (cnt_on == 3) {
                if (is_tri_tri_cutted_hint(vs[0], vs[1], vs[2], tp1, tp2, tp3, CUT_COPLANAR) == CUT_COPLANAR) {
                    result = CUT_COPLANAR;
                    is_cutted = true;
                    is_cut_vs[(j + 1) % 4] = true;
                    is_cut_vs[(j + 2) % 4] = true;
                    is_cut_vs[(j + 3) % 4] = true;
                }
            } else if (cnt_pos > 0 && cnt_neg > 0) {
//                //fortest
//                if(t_id == 3976 && f_id == test_f_id) {
//                    mesh.tets[t_id].print();
//                    cout << oris[0] << " " << oris[1] << " " << oris[2] << " " << oris[3] << endl;
//                    std::vector<int> tmp;
//                    set_intersection(mesh.tet_vertices[mesh.tets[t_id][1]].conn_tets,
//                                     mesh.tet_vertices[mesh.tets[t_id][2]].conn_tets,
//                                     mesh.tet_vertices[mesh.tets[t_id][3]].conn_tets, tmp);
//                    vector_print(tmp);
//                    is_tri_tri_cutted_hint(vs[0], vs[1], vs[2], tp1, tp2, tp3, CUT_FACE, true);
//                    is_tri_tri_cutted_hint(vs[0], vs[1], vs[2], tp1, tp3, tp2, CUT_FACE, true);
//                    is_tri_tri_cutted_hint(vs[0], vs[2], vs[1], tp1, tp2, tp3, CUT_FACE, true);
//                    is_tri_tri_cutted_hint(vs[0], vs[2], vs[1], tp1, tp3, tp2, CUT_FACE, true);
//                    for (int k = 0; k < 3; k++) {
//                        cout << Predicates::orient_3d(tp1, tp2, tp3, vs[k]) << " "
//                             << Predicates::orient_3d(tp1, tp3, tp2, vs[k]) << " "
//                             << Predicates::orient_3d(tp3, tp2, tp1, vs[k]) << endl;
//                    }
//                    pausee();
//                }
//                //fortest

                if (is_tri_tri_cutted_hint(vs[0], vs[1], vs[2], tp1, tp2, tp3, CUT_FACE) == CUT_FACE) {
                    result = CUT_FACE;
                    is_cutted = true;
                    is_cut_vs[(j + 1) % 4] = true;
                    is_cut_vs[(j + 2) % 4] = true;
                    is_cut_vs[(j + 3) % 4] = true;
                }
            } else if (cnt_on == 2 && oris[(j + 1) % 4] == Predicates::ORI_ZERO
                       && oris[(j + 2) % 4] == Predicates::ORI_ZERO) {
                if (is_tri_tri_cutted_hint(vs[0], vs[1], vs[2], tp1, tp2, tp3, CUT_EDGE_0) == CUT_EDGE_0) {
                    result = CUT_EDGE_0;
                    is_cut_vs[(j + 1) % 4] = true;
                    is_cut_vs[(j + 2) % 4] = true;
                }
            } else if (cnt_on == 2 && oris[(j + 2) % 4] == Predicates::ORI_ZERO
                       && oris[(j + 3) % 4] == Predicates::ORI_ZERO) {
                if (is_tri_tri_cutted_hint(vs[0], vs[1], vs[2], tp1, tp2, tp3, CUT_EDGE_1) == CUT_EDGE_1) {
                    result = CUT_EDGE_1;
                    is_cut_vs[(j + 2) % 4] = true;
                    is_cut_vs[(j + 3) % 4] = true;
                }
            } else if (cnt_on == 2 && oris[(j + 3) % 4] == Predicates::ORI_ZERO
                       && oris[(j + 1) % 4] == Predicates::ORI_ZERO) {
                if (is_tri_tri_cutted_hint(vs[0], vs[1], vs[2], tp1, tp2, tp3, CUT_EDGE_2) == CUT_EDGE_2) {
                    result = CUT_EDGE_2;
                    is_cut_vs[(j + 3) % 4] = true;
                    is_cut_vs[(j + 1) % 4] = true;
                }
            }
//            if (is_cut_vs[0] && is_cut_vs[1] && is_cut_vs[2] && is_cut_vs[3])
//                break;

//            //fortest
//            if (f_id == test_f_id && mesh.tets[t_id].find(test_v_id) >= 0) {
//                cout << "input_f " << input_faces[f_id].transpose() << endl;
////                cout << input_vertices[input_faces[f_id][0]].transpose() << endl;
////                cout << input_vertices[input_faces[f_id][1]].transpose() << endl;
////                cout << input_vertices[input_faces[f_id][2]].transpose() << endl;
//                cout << "t " << t_id << ": ";
//                mesh.tets[t_id].print();
//                cout << "j " << j << endl;
//                cout << "cnt_on = " << cnt_on << endl;
//                cout << "cnt_pos = " << cnt_pos << endl;
//                cout << "cnt_neg = " << cnt_neg << endl;
//                cout << "result = " << result << endl;
////                if (cnt_pos > 0 && cnt_neg > 0) {
//                if(t_id == 3976 || t_id == 1016){
//                    cout<<"//////"<<endl;
//                    std::array<Vector3_r, 4> tet_vr;
//                    std::array<Vector3_r, 3> tri_vr;
//                    std::array<Vector3, 4> tet_vf;
//                    std::array<Vector3, 3> tri_vf;
//                    std::array<int, 4> tet_vids;
//                    std::array<int, 3> tri_vids;
//                    for(int k=0;k<4;k++){
//                        for(int r=0;r<3;r++) {
//                            tet_vr[k][r] = mesh.tet_vertices[mesh.tets[t_id][k]].pos[r];
//                            tet_vf[k][r] = mesh.tet_vertices[mesh.tets[t_id][k]].pos[r];
//                        }
//                        tet_vids[k] = mesh.tets[t_id][k];
//                    }
//                    for(int k=0;k<3;k++){
//                        for(int r=0;r<3;r++) {
//                            tri_vr[k][r] = input_vertices[input_faces[f_id][k]][r];
//                            tri_vf[k][r] = input_vertices[input_faces[f_id][k]][r];
//                        }
//                        tri_vids[k] = input_faces[f_id][k];
//                    }
//                    for(int k=0;k<4;k++) {
//                        cout << "tet " << t_id << " face" << k << endl;
//                        cout<<"plane of tet face:"<<endl;
//                        for (int r = 0; r < 3; r++) {
////                            cout<<mesh.tets[t_id][k]<<" "<<mesh.tets[t_id][(k + 1) % 4]<<" "<<mesh.tets[t_id][(k + 2) % 4]<<endl;
////                            cout<<input_faces[f_id][r]<<endl;
//                            cout<<tet_vids[k]<<" "<<tet_vids[(k + 1) % 4]<<" "<<tet_vids[(k + 2) % 4]<<endl;
//                            cout<<tri_vids[r]<<endl;
//                            if(tri_vids[r] == 406 && tet_vids[(k + 1) % 4] == 406){
//                                auto v = tet_vr[(k + 1) % 4] - tri_vr[r];
//                                cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<endl;
//                                if(v[0] != 0)
//                                    cout<<"v0!=0"<<endl;
//                                if(v[1] != 0)
//                                    cout<<"v1!=0"<<endl;
//                                if(v[2] != 0)
//                                    cout<<"v2!=0"<<endl;
//                                pausee();
//                            }
//                            cout << orient_rational(tet_vr[k], tet_vr[(k + 1) % 4], tet_vr[(k + 2) % 4], tri_vr[r])
//                                 << "/ f:";
//                            cout << Predicates::orient_3d(tet_vf[k], tet_vf[(k + 1) % 4], tet_vf[(k + 2) % 4], tri_vf[r])
//                                 << endl;
//                        }
//                        cout<<"plane of tri:"<<endl;
//                        for (int r = 0; r < 3; r++) {
////                            cout<<input_faces[f_id][0]<<" "<<input_faces[f_id][1]<<" "<<input_faces[f_id][2]<<endl;
////                            cout<<mesh.tets[t_id][(k + r) % 4]<<endl;
//                            cout<<tri_vids[0]<<" "<<tri_vids[1]<<" "<<tri_vids[2]<<endl;
//                            cout<<tet_vids[(k + r) % 4]<<endl;
//                            cout << orient_rational(tri_vr[0], tri_vr[1], tri_vr[2], tet_vr[(k + r) % 4]) << "/ f:";
//                            cout << Predicates::orient_3d(tri_vf[0], tri_vf[1], tri_vf[2], tet_vf[(k + r) % 4]) << endl;
//                        }
//                    }
//                    cout<<"//////"<<endl;
//                }
//
//                if (t_id == 3976){
//                    {
//                        Eigen::MatrixXd V(4, 3);
//                        Eigen::MatrixXi F(4, 3);
//                        for (int k = 0; k < 4; k++) {
//                            V.row(k) = mesh.tet_vertices[mesh.tets[t_id][k]].pos;
//                            F.row(k) << (k + 1) % 4, (k + 2) % 4, (k + 3) % 4;
//                        }
//                        igl::writeOFF("test_cut_t_ids1_"+std::to_string(t_id)+".off", V, F);
//                    }
//                }
//                if (t_id == 1016) {
//                    {
//                        Eigen::MatrixXd V(4, 3);
//                        Eigen::MatrixXi F(4, 3);
//                        for (int k = 0; k < 4; k++) {
//                            V.row(k) = mesh.tet_vertices[mesh.tets[t_id][k]].pos;
//                            F.row(k) << (k + 1) % 4, (k + 2) % 4, (k + 3) % 4;
//                        }
//                        igl::writeOFF("test_cut_t_ids1_"+std::to_string(t_id)+".off", V, F);
//                    }
//                    {
//                        Eigen::MatrixXd V(3, 3);
//                        Eigen::MatrixXi F(1, 3);
//                        for (int k = 0; k < 3; k++)
//                            V.row(k) = input_vertices[input_faces[f_id][k]];
//                        F.row(0) << 0, 1, 2;
//                        igl::writeOFF("test_cut_t_ids2_"+std::to_string(t_id)+".off", V, F);
//                    }
////                    pausee();
//                }
//            }
//            //fortest
        }
        if (is_cutted)
            cut_t_ids.push_back(t_id);

        for (int j = 0; j < 4; j++) {
            if (!is_cut_vs[j])
                continue;
            for (int n_t_id: mesh.tet_vertices[mesh.tets[t_id][j]].conn_tets) {
                if (!is_visited[n_t_id]) {
                    is_visited[n_t_id] = true;
                    queue_t_ids.push(n_t_id);
                }
            }
        }
    }

//    //fortest
//    if (cut_t_ids.empty()) {
//        cout << "cut_t_ids.empty()" << endl;
//        cout << "f" << f_id << ": " << input_faces[f_id][0] << " " << input_faces[f_id][1] << " "
//             << input_faces[f_id][2] << endl;
//        std::vector<int> all_cut_t_ids;
//        for (int t_id = 0; t_id < mesh.tets.size(); t_id++) {
//            std::array<int, 4> oris;
//            for (int j = 0; j < 4; j++) {
//                oris[j] = Predicates::orient_3d(vs[0], vs[1], vs[2], mesh.tet_vertices[mesh.tets[t_id][j]].pos);
//            }
//
//            for (int j = 0; j < 4; j++) {
//                int cnt_pos = 0;
//                int cnt_neg = 0;
//                int cnt_on = 0;
//                for (int k = 0; k < 3; k++) {
//                    if (oris[(j + k + 1) % 4] == Predicates::ORI_ZERO)
//                        cnt_on++;
//                    else if (oris[(j + k + 1) % 4] == Predicates::ORI_POSITIVE)
//                        cnt_pos++;
//                    else
//                        cnt_neg++;
//                }
//
//                int result = CUT_EMPTY;
//                auto &tp1 = mesh.tet_vertices[mesh.tets[t_id][(j + 1) % 4]].pos;
//                auto &tp2 = mesh.tet_vertices[mesh.tets[t_id][(j + 2) % 4]].pos;
//                auto &tp3 = mesh.tet_vertices[mesh.tets[t_id][(j + 3) % 4]].pos;
//                if (cnt_on == 3) {
//                    result = is_tri_tri_cutted_hint(vs[0], vs[1], vs[2], tp1, tp2, tp3, CUT_COPLANAR);
//                } else if (cnt_pos > 0 && cnt_neg > 0) {
//                    result = is_tri_tri_cutted_hint(vs[0], vs[1], vs[2], tp1, tp2, tp3, CUT_FACE);
//                }
//                if (result == CUT_EMPTY)
//                    continue;
//
//                all_cut_t_ids.push_back(t_id);
//                break;
//            }
//        }
//        cout << all_cut_t_ids.size() << endl;
//        for (int t_id:all_cut_t_ids) {
//            cout << "t" << t_id << ": ";
//            mesh.tets[t_id].print();
//        }
//        pausee();
//    }
//    //fortest
}

bool floatTetWild::subdivide_tets(int insert_f_id, Mesh& mesh, CutMesh& cut_mesh,
        std::vector<Vector3>& points, std::map<std::array<int, 2>, int>& map_edge_to_intersecting_point,
        std::vector<std::array<std::vector<int>, 4>>& track_surface_fs,
        std::vector<int>& subdivide_t_ids, std::vector<bool>& is_mark_surface,
        std::vector<MeshTet>& new_tets, std::vector<std::array<std::vector<int>, 4>>& new_track_surface_fs,
        std::vector<int>& modified_t_ids) {

    static const std::array<std::array<int, 2>, 6> t_es = {{{{0, 1}}, {{1, 2}}, {{2, 0}}, {{0, 3}}, {{1, 3}}, {{2, 3}}}};
    static const std::array<std::array<int, 3>, 4> t_f_es = {{{{1, 5, 4}}, {{5, 3, 2}}, {{3, 0, 4}}, {{0, 1, 2}}}};
    static const std::array<std::array<int, 3>, 4> t_f_vs = {{{{3, 1, 2}}, {{0, 2, 3}}, {{1, 3, 0}}, {{2, 0, 1}}}};

    //fortest
//    for (auto m: map_edge_to_intersecting_point)
//        cout << (m.first[0]) << " " << (m.first[1]) << ": " << mesh.tet_vertices.size() + m.second << endl;
    //fortest

    covered_tet_fs.clear();
    for (int I = 0; I < subdivide_t_ids.size(); I++) {
        int t_id = subdivide_t_ids[I];
        bool is_mark_sf = is_mark_surface[I];

        //fortest
//        cout << endl << "t_id = " << t_id << endl;
//        cout << "is_mark_sf = " << is_mark_sf << endl;
//        cout << mesh.tets[t_id][0] << " " << mesh.tets[t_id][1] << " " << mesh.tets[t_id][2] << " "
//             << mesh.tets[t_id][3] << endl;
        //fortest

        /////
        std::bitset<6> config_bit;
        std::array<std::pair<int, int>, 6> on_edge_p_ids;
        int cnt = 4;
        for (int i = 0; i < t_es.size(); i++) {
            const auto &le = t_es[i];
            std::array<int, 2> e = {{mesh.tets[t_id][le[0]], mesh.tets[t_id][le[1]]}};
            if (e[0] > e[1])
                std::swap(e[0], e[1]);
            if (map_edge_to_intersecting_point.find(e) == map_edge_to_intersecting_point.end()) {
                on_edge_p_ids[i].first = -1;
                on_edge_p_ids[i].second = -1;
            } else {
                on_edge_p_ids[i].first = cnt++;
                on_edge_p_ids[i].second = map_edge_to_intersecting_point[e];
                config_bit.set(i);
            }
        }
        int config_id = config_bit.to_ulong();
//        cout << "config_id = " << config_id << endl;//fortest
        if (config_id == 0) { //no intersection
            if (is_mark_sf) {
//                cout << "is_mark_surface" << endl;//fortest
                for (int j = 0; j < 4; j++) {
                    int cnt_on = 0;
                    for (int k = 0; k < 3; k++) {
                        myassert(cut_mesh.map_v_ids.find(mesh.tets[t_id][(j + k + 1) % 4]) != cut_mesh.map_v_ids.end(),
                                "cut_mesh.map_v_ids.find(mesh.tets[t_id][(j + k + 1) % 4]) != cut_mesh.map_v_ids.end()!!");//fortest
                        if (cut_mesh.is_v_on_plane(cut_mesh.map_v_ids[mesh.tets[t_id][(j + k + 1) % 4]])) {
                            cnt_on++;
//                            cout << (j + k + 1) % 4 << endl;//fortest
                        }
                    }
//                    cout << "cnt_on = " << cnt_on << endl;//fortest
                    //fortest
                    if(cnt_on == 4){
                        cout<<"cnt_on==4!!"<<endl;
                    }
                    //fortest

                    if (cnt_on == 3) {
                        new_tets.push_back(mesh.tets[t_id]);
                        new_track_surface_fs.push_back(track_surface_fs[t_id]);
                        (new_track_surface_fs.back())[j].push_back(insert_f_id);
                        modified_t_ids.push_back(t_id);

                        covered_tet_fs.push_back({{mesh.tets[t_id][(j+1)%4],mesh.tets[t_id][(j+2)%4],
                                                          mesh.tets[t_id][(j+3)%4]}});

//                        //fortest
//                        int opp_t_id = get_opp_t_id(t_id, j, mesh);
//                        if(std::find(subdivide_t_ids.begin(), subdivide_t_ids.end(), opp_t_id) == subdivide_t_ids.end()) {
//                            cout << "cannot find opp_t_id in subdivide_t_ids!" << endl;
//                            cout << "t_id j: " << t_id << " " << j << endl;
//                            mesh.tets[t_id].print();
//                            cout << "opp_t_id: " << opp_t_id << endl;
//                            mesh.tets[opp_t_id].print();
//                            for (int k = 0; k < 3; k++) {
//                                int lv_id = cut_mesh.map_v_ids[mesh.tets[t_id][(j + k + 1) % 4]];
//                                cout << lv_id << ": " << cut_mesh.to_plane_dists[lv_id] << " "
//                                     << cut_mesh.is_snapped[lv_id] << endl;
//                            }
//                            pausee();
//                        }
//                        //fortest

                        break;
                    }
                }
            }
            continue;
        }

	const auto &configs = CutTable::get_tet_confs(config_id);
	if(configs.empty())
		continue;

        /////
        std::vector<Vector2i> my_diags;
        for (int j = 0; j < 4; j++) {
            std::vector<int> le_ids;
            for (int k = 0; k < 3; k++) {
                if (on_edge_p_ids[t_f_es[j][k]].first < 0)
                    continue;
                le_ids.push_back(k);
            }
            if (le_ids.size() != 2)//no ambiguity
                continue;

            my_diags.emplace_back();
            auto &diag = my_diags.back();
            if (on_edge_p_ids[t_f_es[j][le_ids[0]]].second > on_edge_p_ids[t_f_es[j][le_ids[1]]].second) {
                diag << on_edge_p_ids[t_f_es[j][le_ids[0]]].first, t_f_vs[j][le_ids[0]];
            } else {
                diag << on_edge_p_ids[t_f_es[j][le_ids[1]]].first, t_f_vs[j][le_ids[1]];
            }
            if (diag[0] > diag[1])
                std::swap(diag[0], diag[1]);
        }
        std::sort(my_diags.begin(), my_diags.end(), [](const Vector2i &a, const Vector2i &b) {
            return std::make_tuple(a[0], a[1]) < std::make_tuple(b[0], b[1]);
        });

        /////
        std::map<int, int> map_lv_to_v_id;
        const int v_size = mesh.tet_vertices.size();
        const int vp_size = mesh.tet_vertices.size() + points.size();
        for (int i = 0; i < 4; i++)
            map_lv_to_v_id[i] = mesh.tets[t_id][i];
        cnt = 0;
        for (int i = 0; i < t_es.size(); i++) {
            if (config_bit[i] == 0)
                continue;
            map_lv_to_v_id[4 + cnt] = v_size + on_edge_p_ids[i].second;
            cnt++;
        }

        /////
        auto get_centroid = [&](const std::vector<Vector4i> &config, int lv_id, Vector3 &c) {
            std::vector<int> n_ids;
            for (const auto &tet: config) {
                std::vector<int> tmp;
                for (int j = 0; j < 4; j++) {
                    if (tet[j] != lv_id)
                        tmp.push_back(tet[j]);
                }
                if (tmp.size() == 4)
                    continue;
                n_ids.insert(n_ids.end(), tmp.begin(), tmp.end());
            }
            vector_unique(n_ids);
            c << 0, 0, 0;
            for (int n_id:n_ids) {
                int v_id = map_lv_to_v_id[n_id];
                if (v_id < v_size)
                    c += mesh.tet_vertices[v_id].pos;
                else
                    c += points[v_id - v_size];
            }
            c /= n_ids.size();
        };

        auto check_config = [&](int diag_config_id, std::vector<std::pair<int, Vector3>> &centroids) {
            const std::vector<Vector4i> &config = CutTable::get_tet_conf(config_id, diag_config_id);
            Scalar min_q = -666;
            int cnt = 0;
            std::map<int, int> map_lv_to_c;
            for (const auto &tet: config) {
                std::array<Vector3, 4> vs;
                for (int j = 0; j < 4; j++) {
                    if (map_lv_to_v_id.find(tet[j]) == map_lv_to_v_id.end()) {
                        if (map_lv_to_c.find(tet[j]) == map_lv_to_c.end()) {
                            get_centroid(config, tet[j], vs[j]);
                            map_lv_to_c[tet[j]] = centroids.size();
                            centroids.push_back(std::make_pair(tet[j], vs[j]));
                        }
                        vs[j] = centroids[map_lv_to_c[tet[j]]].second;
                    } else {
                        int v_id = map_lv_to_v_id[tet[j]];
                        if (v_id < v_size)
                            vs[j] = mesh.tet_vertices[v_id].pos;
                        else
                            vs[j] = points[v_id - v_size];
                    }
                }

                Scalar volume = Predicates::orient_3d_volume(vs[0], vs[1], vs[2], vs[3]);

//                //fortest
//                if(volume==0) {
//                    cout<<std::setprecision(16)<<endl;
//                    cout << "volume = " << volume << ",ori = " << Predicates::orient_3d(vs[0], vs[1], vs[2], vs[3]) << endl;
//                    cout << "centroids.size = " << centroids.size() << endl;
//                    cout << "config_id = " << config_id << endl;
//
//                    cout<<"vertices"<<endl;
//                    for(int k=0;k<4;k++){
//                        cout<<"v"<<mesh.tets[t_id][k]<<": "<<mesh.tet_vertices[mesh.tets[t_id][k]].pos.transpose()<<endl;
//                    }
//                    cout<<"intersecting points"<<endl;
//                    for (int i = 0; i < t_es.size(); i++) {
//                        const auto &le = t_es[i];
//                        std::array<int, 2> e = {{mesh.tets[t_id][le[0]], mesh.tets[t_id][le[1]]}};
//                        if (e[0] > e[1])
//                            std::swap(e[0], e[1]);
//                        if (map_edge_to_intersecting_point.find(e) != map_edge_to_intersecting_point.end()) {
//                            cout << e[0] << " " << e[1] << ": p" << map_edge_to_intersecting_point[e] << " "
//                                 << points[map_edge_to_intersecting_point[e]].transpose() << endl;
//                            cout<<"is_mark_sf = "<<is_mark_sf<<endl;
//                            if(cut_mesh.map_v_ids.find(e[0])!=cut_mesh.map_v_ids.end()) {
//                                cout << "e[0] to_plane_dists = " << cut_mesh.to_plane_dists[cut_mesh.map_v_ids[e[0]]]
//                                     << endl;
//                                cout << "e[0] is_snapped = " << cut_mesh.is_snapped[cut_mesh.map_v_ids[e[0]]]
//                                     << endl;
//                            }
//                            if(cut_mesh.map_v_ids.find(e[1])!=cut_mesh.map_v_ids.end()) {
//                                cout << "e[1] to_plane_dists = " << cut_mesh.to_plane_dists[cut_mesh.map_v_ids[e[1]]]
//                                     << endl;
//                                cout << "e[1] is_snapped = " << cut_mesh.is_snapped[cut_mesh.map_v_ids[e[1]]]
//                                     << endl;
//                            }
//
////                            Vector3 p;
////                            Scalar _;
////                            seg_plane_intersection(mesh.tet_vertices[e[0]].pos, mesh.tet_vertices[e[1]].pos,
////                                                   cut_mesh.p_vs[0], cut_mesh.p_n, p, _);
////                            cout << "new p:" << p.transpose() << endl;
////                            cout << (p - mesh.tet_vertices[e[1]].pos).squaredNorm() << " " << SCALAR_ZERO_2 << endl;
//                        }
//                    }
//
//                    pausee();
//                }
//                //fortest

                if (cnt == 0)
                    min_q = volume;
                else if (volume < min_q)
                    min_q = volume;
                cnt++;
            }

            return min_q;
        };

        int diag_config_id = 0;
        std::vector<std::pair<int, Vector3>> centroids;
        if (!my_diags.empty()) {
            const auto &all_diags = CutTable::get_diag_confs(config_id);
            std::vector<std::pair<int, Scalar>> min_qualities(all_diags.size());
            std::vector<std::vector<std::pair<int, Vector3>>> all_centroids(all_diags.size());
            for (int i = 0; i < all_diags.size(); i++) {
                if (my_diags != all_diags[i]) {
                    min_qualities[i] = std::make_pair(i, -1);
                    continue;
                }

                std::vector<std::pair<int, Vector3>> tmp_centroids;
                Scalar min_q = check_config(i, tmp_centroids);
//                if (min_q < SCALAR_ZERO_3)
//                    continue;
                min_qualities[i] = std::make_pair(i, min_q);
                all_centroids[i] = tmp_centroids;
            }
            std::sort(min_qualities.begin(), min_qualities.end(),
                      [](const std::pair<int, Scalar> &a, const std::pair<int, Scalar> &b) {
                          return a.second < b.second;
                      });

            if (min_qualities.back().second < SCALAR_ZERO_3) { // if tet quality is too bad
//                cout<<std::setprecision(16)<<"return 1 "<<min_qualities.back().second<<endl;
                return false;
            }

            diag_config_id = min_qualities.back().first;
            centroids = all_centroids[diag_config_id];
        } else {
            Scalar min_q = check_config(diag_config_id, centroids);
            if (min_q < SCALAR_ZERO_3) {
//                cout<<std::setprecision(16)<<"return 2 "<<min_q<<endl;
                return false;
            }
        }

        for (int i = 0; i < centroids.size(); i++) {
            map_lv_to_v_id[centroids[i].first] = vp_size + i;
            points.push_back(centroids[i].second);
        }

        //add new tets
        const auto &config = CutTable::get_tet_conf(config_id, diag_config_id);
        const auto &new_is_surface_fs = CutTable::get_surface_conf(config_id, diag_config_id);
        const auto &new_local_f_ids = CutTable::get_face_id_conf(config_id, diag_config_id);
        for (int i = 0; i < config.size(); i++) {
            const auto &t = config[i];
            new_tets.push_back(MeshTet(map_lv_to_v_id[t[0]], map_lv_to_v_id[t[1]],
                                       map_lv_to_v_id[t[2]], map_lv_to_v_id[t[3]]));

            //fortest
//            cout << map_lv_to_v_id[t[0]] << " " << map_lv_to_v_id[t[1]] << " " << map_lv_to_v_id[t[2]] << " "
//                 << map_lv_to_v_id[t[3]] << endl;

            new_track_surface_fs.emplace_back();
            for (int j = 0; j < 4; j++) {
                if (new_is_surface_fs[i][j] && is_mark_sf) {
                    (new_track_surface_fs.back())[j].push_back(insert_f_id);

                    covered_tet_fs.push_back({{new_tets.back()[(j + 1) % 4], new_tets.back()[(j + 3) % 4],
                                                      new_tets.back()[(j + 2) % 4]}});
//                    //fortest
//                    cout << "sf " << i << " " << j << endl;
//                    for (int k = 0; k < 3; k++) {
//                        int v_id = new_tets.back()[(j + k + 1) % 4];
//                        Vector3 p;
//                        if (v_id < v_size)
//                            p = mesh.tet_vertices[v_id].pos;
//                        else
//                            p = points[v_id - v_size];
//                        double dist = cut_mesh.get_to_plane_dist(p);
//                        cout << "v_id = " << v_id << endl;
//                        if (dist > mesh.params.eps_2_coplanar) {
//                            cout << "dist > mesh.params.eps_2_coplanar " << dist << endl;
//                            pausee();
//                        }
//                    }
//                    //fortest
                }

                int old_local_f_id = new_local_f_ids[i][j];
                if (old_local_f_id < 0)
                    continue;
                (new_track_surface_fs.back())[j].insert((new_track_surface_fs.back())[j].end(),
                                                        track_surface_fs[t_id][old_local_f_id].begin(),
                                                        track_surface_fs[t_id][old_local_f_id].end());
                (new_tets.back()).is_bbox_fs[j] = mesh.tets[t_id].is_bbox_fs[old_local_f_id];
                (new_tets.back()).is_surface_fs[j] = mesh.tets[t_id].is_surface_fs[old_local_f_id];
                (new_tets.back()).surface_tags[j] = mesh.tets[t_id].surface_tags[old_local_f_id];
            }
        }
        modified_t_ids.push_back(t_id);
    }

    return true;
}

void floatTetWild::pair_track_surface_fs(Mesh &mesh, std::vector<std::array<std::vector<int>, 4>> &track_surface_fs) {
    std::vector<std::array<bool, 4>> is_visited(track_surface_fs.size(), {{false, false, false, false}});
    for (int t_id = 0; t_id < track_surface_fs.size(); t_id++) {
        for (int j = 0; j < 4; j++) {
            //
            if (is_visited[t_id][j])
                continue;
            is_visited[t_id][j] = true;
            if (track_surface_fs[t_id][j].empty())
                continue;
            //
            int opp_t_id = get_opp_t_id(t_id, j, mesh);
            if (opp_t_id < 0)
                continue;
            int k = get_local_f_id(opp_t_id, mesh.tets[t_id][(j + 1) % 4], mesh.tets[t_id][(j + 2) % 4],
                                   mesh.tets[t_id][(j + 3) % 4], mesh);
            is_visited[opp_t_id][k] = true;
            //
            std::sort(track_surface_fs[t_id][j].begin(), track_surface_fs[t_id][j].end());
            std::sort(track_surface_fs[opp_t_id][k].begin(), track_surface_fs[opp_t_id][k].end());
            std::vector<int> f_ids;
            if (track_surface_fs[t_id][j] != track_surface_fs[opp_t_id][k]) {
                std::set_union(track_surface_fs[t_id][j].begin(), track_surface_fs[t_id][j].end(),
                               track_surface_fs[opp_t_id][k].begin(), track_surface_fs[opp_t_id][k].end(),
                               std::back_inserter(f_ids));
                track_surface_fs[t_id][j] = f_ids;
                track_surface_fs[opp_t_id][k] = f_ids;
            }
        }
    }
}

void floatTetWild::find_boundary_edges(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
                                       const std::vector<bool> &is_face_inserted, const std::vector<bool>& old_is_face_inserted,
                                       std::vector<std::pair<std::array<int, 2>, std::vector<int>>> &b_edge_infos,
                                       std::vector<bool>& is_on_cut_edges,
                                       std::vector<std::array<int, 2>>& b_edges) {
    std::vector<std::array<int, 2>> edges;
    std::vector<std::vector<int>> conn_tris(input_vertices.size());
    std::vector<std::vector<int>> uninserted_conn_tris(input_vertices.size());
    for (int i = 0; i < input_faces.size(); i++) {
        if(!is_face_inserted[i]) {///use currently inserted faces as mesh
            for (int j = 0; j < 3; j++)
                uninserted_conn_tris[input_faces[i][j]].push_back(i);
            continue;
        }
        const auto &f = input_faces[i];
        for (int j = 0; j < 3; j++) {
            //edges
            std::array<int, 2> e = {{f[j], f[(j + 1) % 3]}};
            if (e[0] > e[1])
                std::swap(e[0], e[1]);
            edges.push_back(e);
            //conn_tris
            conn_tris[input_faces[i][j]].push_back(i);
        }
    }
    vector_unique(edges);

    int cnt1 = 0;
    int cnt2 = 0;
    for (const auto &e: edges) {
        std::vector<int> n12_f_ids;
        std::set_intersection(conn_tris[e[0]].begin(), conn_tris[e[0]].end(),
                              conn_tris[e[1]].begin(), conn_tris[e[1]].end(), std::back_inserter(n12_f_ids));
        std::vector<int> uninserted_n12_f_ids;
        std::set_intersection(uninserted_conn_tris[e[0]].begin(), uninserted_conn_tris[e[0]].end(),
                              uninserted_conn_tris[e[1]].begin(), uninserted_conn_tris[e[1]].end(),
                              std::back_inserter(uninserted_n12_f_ids));

        bool needs_preserve = false;
        for(int f_id: n12_f_ids){
            if(!old_is_face_inserted[f_id]) {
                needs_preserve = true;
                break;
            }
        }

        if (n12_f_ids.size() == 1) {//open boundary
            b_edges.push_back(e);
            if(needs_preserve) {
                b_edge_infos.push_back(std::make_pair(e, n12_f_ids));
                if(!uninserted_n12_f_ids.empty())
                    is_on_cut_edges.push_back(true);
                else
                    is_on_cut_edges.push_back(false);
            }
            cnt1++;
        } else {
            int f_id = n12_f_ids[0];
            int j = 0;
            for (; j < 3; j++) {
                if ((input_faces[f_id][j] == e[0] && input_faces[f_id][mod3(j + 1)] == e[1])
                    || (input_faces[f_id][j] == e[1] && input_faces[f_id][mod3(j + 1)] == e[0]))
                    break;
            }

            Vector3 n = get_normal(input_vertices[input_faces[f_id][0]], input_vertices[input_faces[f_id][1]],
                                   input_vertices[input_faces[f_id][2]]);
            int t = get_t(input_vertices[input_faces[f_id][0]], input_vertices[input_faces[f_id][1]],
                          input_vertices[input_faces[f_id][2]]);

            bool is_fine = false;
            for (int k = 0; k < n12_f_ids.size(); k++) {
                if (n12_f_ids[k] == f_id)
                    continue;
                Vector3 n1 = get_normal(input_vertices[input_faces[n12_f_ids[k]][0]],
                                        input_vertices[input_faces[n12_f_ids[k]][1]],
                                        input_vertices[input_faces[n12_f_ids[k]][2]]);
                if (abs(n1.dot(n)) < 1 - SCALAR_ZERO) {
                    is_fine = true;
                    break;
                }
            }
            if (is_fine)
                continue;

            is_fine = false;
            int ori;
            for (int k = 0; k < n12_f_ids.size(); k++) {
                for (int r = 0; r < 3; r++) {
                    if (input_faces[n12_f_ids[k]][r] != input_faces[f_id][j] &&
                        input_faces[n12_f_ids[k]][r] != input_faces[f_id][mod3(j + 1)]) {
                        if (k == 0) {
                            ori = Predicates::orient_2d(to_2d(input_vertices[input_faces[f_id][j]], t),
                                                        to_2d(input_vertices[input_faces[f_id][mod3(j + 1)]], t),
                                                        to_2d(input_vertices[input_faces[n12_f_ids[k]][r]], t));
                            break;
                        }
                        int new_ori = Predicates::orient_2d(to_2d(input_vertices[input_faces[f_id][j]], t),
                                                            to_2d(input_vertices[input_faces[f_id][mod3(j + 1)]],
                                                                  t),
                                                            to_2d(input_vertices[input_faces[n12_f_ids[k]][r]], t));
                        if (new_ori != ori)
                            is_fine = true;
                        break;
                    }
                }
                if (is_fine)
                    break;
            }
            if (is_fine)
                continue;

            cnt2++;
            b_edges.push_back(e);
            if(needs_preserve) {
                b_edge_infos.push_back(std::make_pair(e, n12_f_ids));
                if(!uninserted_n12_f_ids.empty())
                    is_on_cut_edges.push_back(true);
                else
                    is_on_cut_edges.push_back(false);
            }
        }
    }

    cout << "#boundary_e1 = " << cnt1 << endl;
    cout << "#boundary_e2 = " << cnt2 << endl;
}
//double time_e1 = 0;
//double time_e2 = 0;
//double time_e3 = 0;
//double time_e4 = 0;

bool floatTetWild::insert_boundary_edges(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
                                         std::vector<std::pair<std::array<int, 2>, std::vector<int>>>& b_edge_infos,
                                         std::vector<bool>& is_on_cut_edges,
                                         std::vector<std::array<std::vector<int>, 4>> &track_surface_fs, Mesh& mesh,
                                         AABBWrapper &tree,
                                         std::vector<bool> &is_face_inserted, bool is_again,
                                         std::vector<std::array<int, 3>>& known_surface_fs,
                                         std::vector<std::array<int, 3>>& known_not_surface_fs) {
//    time_e1 = 0;
//    time_e2 = 0;
//    time_e3 = 0;
//    time_e4 = 0;
    igl::Timer timer;

    //fortest
    auto check_corvered_area = [&](int I, const std::vector<std::array<int, 3>>& cut_fs){
        int f_id = b_edge_infos[I].second[0];
        {
            auto &infos = cut_fs;
            Eigen::MatrixXd V(infos.size() * 3, 3), C(infos.size() * 3, 3);
            Eigen::MatrixXi F(infos.size(), 3);
            for (int i = 0; i < infos.size(); i++) {
                for (int j = 0; j < 3; j++) {
                    V.row(i * 3 + j) = mesh.tet_vertices[infos[i][j]].pos;
                    C.row(i * 3 + j) << 0, 0, 255;
                }
                F.row(i) << i * 3, i * 3 + 1, i * 3 + 2;
            }
            igl::writeOFF("_covered_tet_fs_" + std::to_string(f_id) + ".off", V, F, C);
        }
        {
            Eigen::MatrixXd V(3, 3), C(3, 3);
            Eigen::MatrixXi F(1, 3);
            for (int j = 0; j < 3; j++) {
                V.row(j) = input_vertices[input_faces[f_id][j]];
                C.row(j) << 255, 0, 0;
            }
            F.row(0) << 0, 1, 2;
            igl::writeOFF("_input_f_" + std::to_string(f_id) + ".off", V, F, C);
        }
        pausee();
    };
    //fortest

    auto mark_known_surface_fs = [&](const std::array<int, 3> &f, int tag) {
        std::vector<int> n_t_ids;
        set_intersection(mesh.tet_vertices[f[0]].conn_tets, mesh.tet_vertices[f[1]].conn_tets,
                         mesh.tet_vertices[f[2]].conn_tets, n_t_ids);
        if (n_t_ids.size() != 2)//todo:?????
            return;

        for (int t_id:n_t_ids) {
            int j = get_local_f_id(t_id, f[0], f[1], f[2], mesh);
            mesh.tets[t_id].is_surface_fs[j] = tag;
        }
    };


    auto record_boundary_info = [&](const std::vector<Vector3> &points, const std::vector<int> &snapped_v_ids,
                                    const std::array<int, 2> &e, bool is_on_cut) {
        const int tet_vertices_size = mesh.tet_vertices.size();
        for (int i = points.size(); i > 0; i--) {
            mesh.tet_vertices[tet_vertices_size - i].is_on_boundary = true;
            mesh.tet_vertices[tet_vertices_size - i].is_on_cut = is_on_cut;
        }

        for (int v_id: snapped_v_ids) {
            Scalar dis_2 = p_seg_squared_dist_3d(mesh.tet_vertices[v_id].pos, input_vertices[e[0]],
                                                 input_vertices[e[1]]);
            if (dis_2 <= mesh.params.eps_2) {
                mesh.tet_vertices[v_id].is_on_boundary = true;
                mesh.tet_vertices[v_id].is_on_cut = is_on_cut;
            }
        }

//        b_edges.push_back(e);

//        //fortest
//        for (int v_id: snapped_v_ids) {
//            Scalar dis_2 = p_seg_squared_dist_3d(mesh.tet_vertices[v_id].pos, input_vertices[e[0]],
//                                                 input_vertices[e[1]]);
//            if (dis_2 <= mesh.params.eps_2)
//                fout_v << mesh.tet_vertices[v_id].pos[0] << " "
//                     << mesh.tet_vertices[v_id].pos[1] << " "
//                     << mesh.tet_vertices[v_id].pos[2] << endl;
//        }
//        //fortest
    };


    timer.start();
    std::vector<std::vector<std::pair<int, int>>> covered_fs_infos(input_faces.size());
    for (int i = 0; i < track_surface_fs.size(); i++) {
        if (mesh.tets[i].is_removed)
            continue;
        for (int j = 0; j < 4; j++) {
            for (int f_id: track_surface_fs[i][j])
                covered_fs_infos[f_id].push_back(std::make_pair(i, j));
        }
    }
    logger().info("time1 = {}", timer.getElapsedTime());

    bool is_all_inserted = true;
    int cnt = 0;
    double time2 = 0;//fortest
    double time3 = 0;
    double time4 = 0;
    double time5 = 0;
    double time6 = 0;
    for (int I = 0; I < b_edge_infos.size(); I++) {
        const auto &e = b_edge_infos[I].first;
        auto &n_f_ids = b_edge_infos[I].second;///it is sorted
        bool is_on_cut = is_on_cut_edges[I];
        if(!is_again) {
            mesh.tet_vertices[e[0]].is_on_boundary = true;
            mesh.tet_vertices[e[1]].is_on_boundary = true;
            mesh.tet_vertices[e[0]].is_on_cut = is_on_cut;
            mesh.tet_vertices[e[1]].is_on_cut = is_on_cut;
        }
//        b_edges.push_back(e);

        timer.start();
        ///double check neighbors
        for (int i = 0; i < n_f_ids.size(); i++) {
            if (!is_face_inserted[n_f_ids[i]]) {
                n_f_ids.erase(n_f_ids.begin() + i);
                i--;
                break;
            }
        }
        time2+=timer.getElapsedTime();
        if (n_f_ids.empty()) {
            cout << "FAIL n_f_ids.empty()" << endl;
            continue;
        }

        timer.start();
        ///compute intersection
        std::vector<Vector3> points;
        std::map<std::array<int, 2>, int> map_edge_to_intersecting_point;
        std::vector<int> snapped_v_ids;
        std::vector<std::array<int, 3>> cut_fs;
        if (!insert_boundary_edges_get_intersecting_edges_and_points(covered_fs_infos,
                                                                     input_vertices, input_faces, e, n_f_ids,
                                                                     track_surface_fs,
                                                                     mesh, points, map_edge_to_intersecting_point,
                                                                     snapped_v_ids, cut_fs,
                                                                     is_again)) {
            for (int f_id: n_f_ids)
                is_face_inserted[f_id] = false;
            is_all_inserted = false;

            cout << "FAIL insert_boundary_edges_get_intersecting_edges_and_points" << endl;
            time3+=timer.getElapsedTime();
            continue;
        }
        time3+=timer.getElapsedTime();
        if (points.empty()) { ///if all snapped
            record_boundary_info(points, snapped_v_ids, e, is_on_cut);
            cnt++;
            continue;
        }

        timer.start();
        ///subdivision
        std::vector<int> cut_t_ids;
        for (const auto &m: map_edge_to_intersecting_point) {
            const auto &tet_e = m.first;
            std::vector<int> tmp;
            set_intersection(mesh.tet_vertices[tet_e[0]].conn_tets, mesh.tet_vertices[tet_e[1]].conn_tets, tmp);
            cut_t_ids.insert(cut_t_ids.end(), tmp.begin(), tmp.end());
        }
        vector_unique(cut_t_ids);
        std::vector<bool> is_mark_surface(cut_t_ids.size(), false);
        CutMesh empty_cut_mesh(mesh, Vector3(0, 0, 0), std::array<Vector3, 3>());
        //
        std::vector<MeshTet> new_tets;
        std::vector<std::array<std::vector<int>, 4>> new_track_surface_fs;
        std::vector<int> modified_t_ids;
        if (!subdivide_tets(-1, mesh, empty_cut_mesh, points, map_edge_to_intersecting_point, track_surface_fs,
                            cut_t_ids, is_mark_surface,
                            new_tets, new_track_surface_fs, modified_t_ids)) {
            bool is_inside_envelope = true;
            for (auto &f: cut_fs) {
#ifdef NEW_ENVELOPE
                if(tree.is_out_sf_envelope_exact({{mesh.tet_vertices[f[0]].pos, mesh.tet_vertices[f[1]].pos,
                                                    mesh.tet_vertices[f[2]].pos}})){
                    is_inside_envelope = false;
                    break;
                }
#else
    #ifdef STORE_SAMPLE_POINTS
                    std::vector<GEO::vec3> ps;
                    sample_triangle({{mesh.tet_vertices[f[0]].pos, mesh.tet_vertices[f[1]].pos,
                                             mesh.tet_vertices[f[2]].pos}}, ps, mesh.params.dd);
                    if (tree.is_out_sf_envelope(ps, mesh.params.eps_2)) {
    #else
                    GEO::index_t prev_facet = GEO::NO_FACET;
                    if(sample_triangle_and_check_is_out({{mesh.tet_vertices[f[0]].pos, mesh.tet_vertices[f[1]].pos,
                                             mesh.tet_vertices[f[2]].pos}}, mesh.params.dd, mesh.params.eps_2, tree, prev_facet)){
    #endif
                        is_inside_envelope = false;
                        break;
                    }
#endif
            }
            if (!is_inside_envelope) {
                for (int f_id: n_f_ids)
                    is_face_inserted[f_id] = false;
                for (auto &f: cut_fs) {
                    mark_known_surface_fs(f, KNOWN_NOT_SURFACE);
                }
                known_not_surface_fs.insert(known_not_surface_fs.end(), cut_fs.begin(), cut_fs.end());
                cout << "FAIL subdivide_tets" << endl;
            } else {
                for (auto &f: cut_fs)
                    mark_known_surface_fs(f, KNOWN_SURFACE);
                known_surface_fs.insert(known_surface_fs.end(), cut_fs.begin(), cut_fs.end());
                cout << "SEMI-FAIL subdivide_tets" << endl;
            }

            is_all_inserted = false;//unless now
            time4+=timer.getElapsedTime();
            continue;
        }
        time4+=timer.getElapsedTime();

        //
        timer.start();
        push_new_tets(mesh, track_surface_fs, points, new_tets, new_track_surface_fs, modified_t_ids, is_again);
        time5+=timer.getElapsedTime();

        //
        ///mark boundary vertices
        ///notice, here we assume points list is inserted in the end of mesh.tet_vertices
        timer.start();
        record_boundary_info(points, snapped_v_ids, e, is_on_cut);
        time6+=timer.getElapsedTime();
        cnt++;
    }

    logger().info("uninsert boundary #e = {}/{}", b_edge_infos.size() - cnt, b_edge_infos.size());
    logger().info("time2 = {}", time2);
    logger().info("time3 = {}", time3);
    logger().info("time4 = {}", time4);
    logger().info("time5 = {}", time5);
    logger().info("time6 = {}", time6);

//    logger().info("time_e1 = {}", time_e1);
//    logger().info("time_e2 = {}", time_e2);
//    logger().info("time_e3 = {}", time_e3);
//    //fortest
//    std::ofstream fout_e("b_edges1.obj");
////    std::ofstream fout_v("b_edges1_points.xyz");
//    for (auto &e: b_edges) {
//        fout_e << "v " << input_vertices[e[0]][0] << " "
//               << input_vertices[e[0]][1] << " "
//               << input_vertices[e[0]][2] << endl;
//        fout_e << "v " << input_vertices[e[1]][0] << " "
//               << input_vertices[e[1]][1] << " "
//               << input_vertices[e[1]][2] << endl;
//    }
//    for (int i = 0; i < b_edges.size(); i++) {
//        fout_e << "l " << i * 2 + 1 << " " << i * 2 + 2 << endl;
//    }
//    fout_e.close();
////    fout_v.close();
//    //fortest

//    cout << "b_edges.size = " << b_edges.size() << endl;

    return is_all_inserted;
}

bool floatTetWild::insert_boundary_edges_get_intersecting_edges_and_points(
        const std::vector<std::vector<std::pair<int, int>>>& covered_fs_infos,
        const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
        const std::array<int, 2> &e, const std::vector<int> &n_f_ids,
        std::vector<std::array<std::vector<int>, 4>> &track_surface_fs, Mesh &mesh,
        std::vector<Vector3>& points, std::map<std::array<int, 2>, int>& map_edge_to_intersecting_point,
        std::vector<int>& snapped_v_ids, std::vector<std::array<int, 3>>& cut_fs,
        bool is_again) {

//    igl::Timer timer;

    auto is_cross = [](int a, int b) {
        if ((a == Predicates::ORI_POSITIVE && b == Predicates::ORI_NEGATIVE)
            || (a == Predicates::ORI_NEGATIVE && b == Predicates::ORI_POSITIVE))
            return true;
        return false;
    };

    int t = get_t(input_vertices[input_faces[n_f_ids.front()][0]],
                  input_vertices[input_faces[n_f_ids.front()][1]],
                  input_vertices[input_faces[n_f_ids.front()][2]]);
    std::array<Vector2, 2> evs_2d = {{to_2d(input_vertices[e[0]], t), to_2d(input_vertices[e[1]], t)}};
    Vector3 n = (input_vertices[input_faces[n_f_ids.front()][1]] -
                 input_vertices[input_faces[n_f_ids.front()][0]]).cross(
            input_vertices[input_faces[n_f_ids.front()][2]] - input_vertices[input_faces[n_f_ids.front()][0]]);
    n.normalize();
    const Vector3 &pp = input_vertices[input_faces[n_f_ids.front()][0]];

//    timer.start();
    std::vector<bool> is_visited(mesh.tets.size(), false);
    std::queue<int> t_ids_queue;
    ///find seed t_ids
    if (!is_again) {
        std::vector<int> t_ids;
        for (int f_id: n_f_ids) {
            for (const auto &info: covered_fs_infos[f_id])
                t_ids.push_back(info.first);
        }
        for (int v_id: e)
            t_ids.insert(t_ids.end(), mesh.tet_vertices[v_id].conn_tets.begin(),
                         mesh.tet_vertices[v_id].conn_tets.end());
        vector_unique(t_ids);
        for (int t_id: t_ids) {
            t_ids_queue.push(t_id);
            is_visited[t_id] = true;
        }
    } else {
        Vector3 min_e, max_e;
        get_bbox_face(input_vertices[e[0]], input_vertices[e[0]], input_vertices[e[1]], min_e, max_e);

        std::vector<int> t_ids;
        for (int f_id: n_f_ids) {
            for (const auto &info: covered_fs_infos[f_id])
                t_ids.push_back(info.first);
        }
#ifdef FLOAT_TETWILD_USE_TBB
        tbb::concurrent_vector<int> tbb_t_ids;
        tbb::parallel_for(size_t(0), mesh.tets.size(), [&](size_t t_id){
            if (mesh.tets[t_id].is_removed)
                return;
            if (track_surface_fs[t_id][0].empty() && track_surface_fs[t_id][1].empty()
                && track_surface_fs[t_id][2].empty() && track_surface_fs[t_id][3].empty())
                return;

            Vector3 min_t, max_t;
            get_bbox_tet(mesh.tet_vertices[mesh.tets[t_id][0]].pos, mesh.tet_vertices[mesh.tets[t_id][1]].pos,
                         mesh.tet_vertices[mesh.tets[t_id][2]].pos, mesh.tet_vertices[mesh.tets[t_id][3]].pos,
                         min_t, max_t);
            if (!is_bbox_intersected(min_e, max_e, min_t, max_t))
                return;
            tbb_t_ids.push_back(t_id);
        });
        t_ids.insert(t_ids.end(), tbb_t_ids.begin(), tbb_t_ids.end());
#else
        for (int t_id = 0; t_id < mesh.tets.size(); t_id++) {
            if (mesh.tets[t_id].is_removed)
                continue;
            if (track_surface_fs[t_id][0].empty() && track_surface_fs[t_id][1].empty()
                && track_surface_fs[t_id][2].empty() && track_surface_fs[t_id][3].empty())
                continue;

            Vector3 min_t, max_t;
            get_bbox_tet(mesh.tet_vertices[mesh.tets[t_id][0]].pos, mesh.tet_vertices[mesh.tets[t_id][1]].pos,
                         mesh.tet_vertices[mesh.tets[t_id][2]].pos, mesh.tet_vertices[mesh.tets[t_id][3]].pos,
                         min_t, max_t);
            if (!is_bbox_intersected(min_e, max_e, min_t, max_t))
                continue;

            t_ids.push_back(t_id);
        }
#endif

        vector_unique(t_ids);
        for (int t_id: t_ids) {
            t_ids_queue.push(t_id);
            is_visited[t_id] = true;
        }
    }
//    time_e1+=timer.getElapsedTimeInSec();

//    timer.start();
    std::vector<int> v_oris(mesh.tet_vertices.size(), Predicates::ORI_UNKNOWN);
    while (!t_ids_queue.empty()) {
        int t_id = t_ids_queue.front();
        t_ids_queue.pop();

        std::array<bool, 4> is_cut_vs = {{false, false, false, false}};
        for (int j = 0; j < 4; j++) {
            ///check if contains
            if (track_surface_fs[t_id][j].empty())
                continue;
            std::sort(track_surface_fs[t_id][j].begin(), track_surface_fs[t_id][j].end());
            std::vector<int> tmp;
            std::set_intersection(track_surface_fs[t_id][j].begin(), track_surface_fs[t_id][j].end(),
                                  n_f_ids.begin(), n_f_ids.end(), std::back_inserter(tmp));
            if (tmp.empty())
                continue;

            ///check if cut through
            //check tri side of seg
            std::array<int, 3> f_v_ids = {{mesh.tets[t_id][(j + 1) % 4], mesh.tets[t_id][(j + 2) % 4],
                                                  mesh.tets[t_id][(j + 3) % 4]}};
//            std::array<Vector2, 3> fvs_2d = {{to_2d(mesh.tet_vertices[f_v_ids[0]].pos, t),
//                                                     to_2d(mesh.tet_vertices[f_v_ids[1]].pos, t),
//                                                     to_2d(mesh.tet_vertices[f_v_ids[2]].pos, t)}};
            std::array<Vector2, 3> fvs_2d = {{to_2d(mesh.tet_vertices[f_v_ids[0]].pos, n, pp, t),
                                                     to_2d(mesh.tet_vertices[f_v_ids[1]].pos, n, pp, t),
                                                     to_2d(mesh.tet_vertices[f_v_ids[2]].pos, n, pp, t)}};
            int cnt_pos = 0;
            int cnt_neg = 0;
            int cnt_on = 0;
            for (int k = 0; k < 3; k++) {
                int &ori = v_oris[f_v_ids[k]];
                if (ori == Predicates::ORI_UNKNOWN)
                    ori = Predicates::orient_2d(evs_2d[0], evs_2d[1], fvs_2d[k]);
                if (ori == Predicates::ORI_ZERO) {
                    cnt_on++;
                } else {
                    Scalar dis_2 = p_seg_squared_dist_3d(mesh.tet_vertices[f_v_ids[k]].pos, input_vertices[e[0]],
                                                         input_vertices[e[1]]);
                    if (dis_2 < mesh.params.eps_2_coplanar) {
                        ori = Predicates::ORI_ZERO;
                        cnt_on++;
                        continue;
                    }
                    if (ori == Predicates::ORI_POSITIVE)
                        cnt_pos++;
                    else
                        cnt_neg++;
                }
            }
            if (cnt_on >= 2) {
                cut_fs.push_back(f_v_ids);
                std::sort(cut_fs.back().begin(), cut_fs.back().end());
                continue;
            }
            if (cnt_neg == 0 || cnt_pos == 0)
                continue;

            //check tri edge - seg intersection
            bool is_intersected = false;
            for (auto &p: evs_2d) { ///first check if endpoints are contained inside the triangle
                if (is_p_inside_tri_2d(p, fvs_2d)) {
                    is_intersected = true;
                    break;
                }
            }
            if (!is_intersected) { ///then check if there's intersection
                for (int k = 0; k < 3; k++) {
                    //if cross
                    if (!is_cross(v_oris[f_v_ids[k]], v_oris[f_v_ids[(k + 1) % 3]]))
                        continue;
                    //if already know intersect
                    std::array<int, 2> tri_e = {{f_v_ids[k], f_v_ids[(k + 1) % 3]}};
                    if (tri_e[0] > tri_e[1])
                        std::swap(tri_e[0], tri_e[1]);
                    if (map_edge_to_intersecting_point.find(tri_e) != map_edge_to_intersecting_point.end()) {
                        is_intersected = true;
                        break;
                    }
                    //if intersect
                    double t2 = -1;
                    if (seg_seg_intersection_2d(evs_2d, {{fvs_2d[k], fvs_2d[(k + 1) % 3]}}, t2)) {
                        Vector3 p = (1 - t2) * mesh.tet_vertices[f_v_ids[k]].pos
                                    + t2 * mesh.tet_vertices[f_v_ids[(k + 1) % 3]].pos;
                        double dis1 = (p - mesh.tet_vertices[f_v_ids[k]].pos).squaredNorm();
                        double dis2 = (p - mesh.tet_vertices[f_v_ids[(k + 1) % 3]].pos).squaredNorm();
//                        if (dis1 < SCALAR_ZERO_2) {
                        if (dis1 < mesh.params.eps_2_coplanar) {
                            v_oris[f_v_ids[k]] = Predicates::ORI_ZERO;
                            is_intersected = true;
                            break;
                        }
//                        if (dis2 < SCALAR_ZERO_2) {
                        if (dis2 < mesh.params.eps_2_coplanar) {
                            v_oris[f_v_ids[k]] = Predicates::ORI_ZERO;
                            is_intersected = true;
                            break;
                        }
                        points.push_back(p);
                        map_edge_to_intersecting_point[tri_e] = points.size() - 1;
                        is_intersected = true;
                        break;
                    }
                }
                if (!is_intersected) /// no need to return false here
                    continue;
            }

            std::sort(f_v_ids.begin(), f_v_ids.end());
            cut_fs.push_back(f_v_ids);
            is_cut_vs[(j + 1) % 4] = true;
            is_cut_vs[(j + 2) % 4] = true;
            is_cut_vs[(j + 3) % 4] = true;
        }
        for (int j = 0; j < 4; j++) {
            if (!is_cut_vs[j])
                continue;
            for (int n_t_id: mesh.tet_vertices[mesh.tets[t_id][j]].conn_tets) {
                if (!is_visited[n_t_id]) {
                    t_ids_queue.push(n_t_id);
                    is_visited[n_t_id] = true;
                }
            }
        }
    }
    vector_unique(cut_fs);
//    time_e2+=timer.getElapsedTimeInSec();

//    timer.start();
    std::vector<std::array<int, 2>> tet_edges;
    for (const auto &f:cut_fs) {
        for (int j = 0; j < 3; j++) {
            if (f[j] < f[(j + 1) % 3])
                tet_edges.push_back({{f[j], f[(j + 1) % 3]}});
            else
                tet_edges.push_back({{f[(j + 1) % 3], f[j]}});
        }
    }
    vector_unique(tet_edges);
    //
    for (const auto &tet_e:tet_edges) {
        bool is_snapped = false;
        for (int j = 0; j < 2; j++) {
            if (v_oris[tet_e[j]] == Predicates::ORI_ZERO) {
                snapped_v_ids.push_back(tet_e[j]);
                is_snapped = true;
            }
        }
        if (is_snapped)
            continue;
        if (map_edge_to_intersecting_point.find(tet_e) != map_edge_to_intersecting_point.end())
            continue;
        std::array<Vector2, 2> tri_evs_2d = {{to_2d(mesh.tet_vertices[tet_e[0]].pos, n, pp, t),
                                                     to_2d(mesh.tet_vertices[tet_e[1]].pos, n, pp, t)}};
        Scalar t_seg = -1;
        if (seg_line_intersection_2d(tri_evs_2d, evs_2d, t_seg)) {
            points.push_back((1 - t_seg) * mesh.tet_vertices[tet_e[0]].pos
                             + t_seg * mesh.tet_vertices[tet_e[1]].pos);
            map_edge_to_intersecting_point[tet_e] = points.size() - 1;
        }
    }
    vector_unique(snapped_v_ids);
//    time_e3+=timer.getElapsedTimeInSec();

//    //fortest
//    Eigen::MatrixXd V(cut_fs.size() * 3, 3), C(cut_fs.size() * 3, 3);
//    Eigen::MatrixXi F(cut_fs.size(), 3);
//    for (int i = 0; i < cut_fs.size(); i++) {
//        for (int j = 0; j < 3; j++) {
//            V.row(i * 3 + j) = mesh.tet_vertices[cut_fs[i][j]].pos;
//            C.row(i * 3 + j) << 0, 0, 1;
//        }
//        F.row(i) << i * 3, i * 3 + 1, i * 3 + 2;
//    }
//    igl::writeOFF("test_f.off", V, F, C);
////            Eigen::MatrixXd V(3, 3), C(3, 3);
////            Eigen::MatrixXi F(1, 3);
////            for (int j = 0; j < 3; j++) {
////                V.row(j) = mesh.tet_vertices[cut_fs[i][j]].pos;
////                C.row(j) << 0, 0, 1;
////            }
////            F.row(i) << 0, 1, 2;
////            igl::writeOFF("test_f.off", V, F, C);
//    //
//    std::ofstream fout("test_e.obj");
//    fout << "v " << input_vertices[e[0]].transpose() << endl;
//    fout << "v " << input_vertices[e[1]].transpose() << endl;
//    fout << "l 1 2" << endl;
//    fout.close();
//    //
//    pausee();
//    //fortest

    return true;
}

void floatTetWild::mark_surface_fs(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
                                   const std::vector<int> &input_tags,
                                   std::vector<std::array<std::vector<int>, 4>> &track_surface_fs,
                                   const std::vector<bool> &is_face_inserted,
                                   const std::vector<std::array<int, 3>>& known_surface_fs,
                                   const std::vector<std::array<int, 3>>& known_not_surface_fs,
                                   std::vector<std::array<int, 2>>& b_edges,
                                   Mesh &mesh, AABBWrapper &tree) {

    auto is_on_bounded_side = [&](const std::array<Vector2, 3> &ps_2d, const Vector2 &c) {
        int cnt_pos = 0;
        int cnt_neg = 0;
        for (int j = 0; j < 3; j++) {
            int ori = Predicates::orient_2d(ps_2d[j], ps_2d[(j + 1) % 3], c);
            if (ori == Predicates::ORI_POSITIVE)
                cnt_pos++;
            else if (ori == Predicates::ORI_NEGATIVE)
                cnt_neg++;
        }

        if (cnt_neg > 0 && cnt_pos > 0)
            return false;
        return true;
    };

    std::vector<std::array<bool, 4>> is_visited(track_surface_fs.size(), {{false, false, false, false}});
    for (int t_id = 0; t_id < track_surface_fs.size(); t_id++) {
        for (int j = 0; j < 4; j++) {
//            //fortest
//            if (track_surface_fs[t_id][j].size() > 0)
////            if (std::find(track_surface_fs[t_id][j].begin(), track_surface_fs[t_id][j].end(), III) != track_surface_fs[t_id][j].end())
//                mesh.tets[t_id].is_surface_fs[j] = -1;
//            continue;
//            //fortest

            //
//            if (mesh.tets[t_id].is_surface_fs[j] != NOT_SURFACE || is_visited[t_id][j])
//                continue;
//            is_visited[t_id][j] = true;
//            if (track_surface_fs[t_id][j].empty())
//                continue;
//            //
//            int opp_t_id = get_opp_t_id(t_id, j, mesh);
//            if (opp_t_id < 0)
//                continue;
//            int k = get_local_f_id(opp_t_id, mesh.tets[t_id][(j + 1) % 4], mesh.tets[t_id][(j + 2) % 4],
//                                   mesh.tets[t_id][(j + 3) % 4], mesh);
//            is_visited[opp_t_id][k] = true;
//            //
//            std::sort(track_surface_fs[t_id][j].begin(), track_surface_fs[t_id][j].end());
//            std::sort(track_surface_fs[opp_t_id][k].begin(), track_surface_fs[opp_t_id][k].end());
//            std::vector<int> f_ids;
//            if (track_surface_fs[t_id][j] != track_surface_fs[opp_t_id][k])
//                std::set_union(track_surface_fs[t_id][j].begin(), track_surface_fs[t_id][j].end(),
//                               track_surface_fs[opp_t_id][k].begin(), track_surface_fs[opp_t_id][k].end(),
//                               std::back_inserter(f_ids));
//            else
//                f_ids = track_surface_fs[t_id][j];

            int ff_id = -1;
            int opp_t_id = -1;
            int k = -1;
            if (mesh.tets[t_id].is_surface_fs[j] == KNOWN_SURFACE
                || mesh.tets[t_id].is_surface_fs[j] == KNOWN_NOT_SURFACE) {
                opp_t_id = get_opp_t_id(t_id, j, mesh);
                if (opp_t_id < 0) {
                    mesh.tets[t_id].is_surface_fs[j] = NOT_SURFACE;
                    continue;
                }
                k = get_local_f_id(opp_t_id, mesh.tets[t_id][(j + 1) % 4], mesh.tets[t_id][(j + 2) % 4],
                                   mesh.tets[t_id][(j + 3) % 4], mesh);
                is_visited[t_id][j] = true;
                is_visited[opp_t_id][k] = true;
                if (mesh.tets[t_id].is_surface_fs[j] == KNOWN_NOT_SURFACE || track_surface_fs[t_id][j].empty()) {
                    mesh.tets[t_id].is_surface_fs[j] = NOT_SURFACE;
                    mesh.tets[opp_t_id].is_surface_fs[k] = NOT_SURFACE;
                    continue;
                } else
                    ff_id = track_surface_fs[t_id][j].front();
            } else {
                if (mesh.tets[t_id].is_surface_fs[j] != NOT_SURFACE || is_visited[t_id][j])
                    continue;
                is_visited[t_id][j] = true;
                if (track_surface_fs[t_id][j].empty())
                    continue;

                auto &f_ids = track_surface_fs[t_id][j];

                auto &tp1_3d = mesh.tet_vertices[mesh.tets[t_id][(j + 1) % 4]].pos;
                auto &tp2_3d = mesh.tet_vertices[mesh.tets[t_id][(j + 2) % 4]].pos;
                auto &tp3_3d = mesh.tet_vertices[mesh.tets[t_id][(j + 3) % 4]].pos;

                //
//                static const Scalar dist_2_max = mesh.params.bbox_diag_length * mesh.params.bbox_diag_length;
//                std::array<Scalar, 3> min_dist_2 = {{dist_2_max, dist_2_max, dist_2_max}};
//                std::array<GEO::vec3, 3> t_vs_geo = {{GEO::vec3(tp1_3d[0], tp1_3d[1], tp1_3d[2]),
//                                                             GEO::vec3(tp2_3d[0], tp2_3d[1], tp2_3d[2]),
//                                                             GEO::vec3(tp3_3d[0], tp3_3d[1], tp3_3d[2])}};
//                for (int f_id: f_ids) {
//                    std::array<GEO::vec3, 3> f_vs_geo;
//                    for (int k = 0; k < 3; k++) {
//                        f_vs_geo[k] = GEO::vec3(input_vertices[input_faces[f_id][k]][0],
//                                                input_vertices[input_faces[f_id][k]][1],
//                                                input_vertices[input_faces[f_id][k]][2]);
//                    }
//                    for (int k = 0; k < 3; k++) {
//                        double dist_2 = GEO::Geom::point_triangle_squared_distance(t_vs_geo[k], f_vs_geo[0],
//                                                                                   f_vs_geo[1], f_vs_geo[2]);
//                        if (dist_2 < min_dist_2[k])
//                            min_dist_2[k] = dist_2;
//                    }
//                }
//                double eps = mesh.params.eps * mesh.params.eps;
//                if (min_dist_2[0] <= eps && min_dist_2[1] <= eps && min_dist_2[2] <= eps)
//                    ff_id = track_surface_fs[t_id][j].front();
//                else {
//                    //fortest
////                    int t = get_t(tp1_3d, tp2_3d, tp3_3d);
////                    std::array<Vector2, 3> tps_2d = {{to_2d(tp1_3d, t), to_2d(tp2_3d, t), to_2d(tp3_3d, t)}};
////                    std::array<Vector2, 4> cs = {{(tps_2d[0] + tps_2d[1] + tps_2d[2]) / 3}};
////                    for (int k = 0; k < 3; k++)
////                        cs[k + 1] = (tps_2d[k] + cs[0]) / 3;
////
////                    std::array<Vector2, 3> ps_2d;
////                    for (int f_id: f_ids) {
////                        if (!is_face_inserted[f_id])
////                            continue;
////
////                        ps_2d = {{to_2d(input_vertices[input_faces[f_id][0]], t),
////                                         to_2d(input_vertices[input_faces[f_id][1]], t),
////                                         to_2d(input_vertices[input_faces[f_id][2]], t)}};
////
////                        bool is_in = false;
////                        for (auto &c:cs) {
////                            if (is_on_bounded_side(ps_2d, c)) {
////                                is_in = true;
////                                break;
////                            }
////                        }
////                        if (is_in) {
////                            ff_id = f_id;
////                            break;
////                        }
////                    }
////                    if(ff_id>0) {
////                        cout<<t_id<<endl;
////                        cout << min_dist_2[0] << " " << min_dist_2[1] << " " << min_dist_2[2] << endl;
////                        cout << eps << endl;
////                        {
////                            Eigen::MatrixXd V(f_ids.size() * 3, 3), C(f_ids.size() * 3, 3);
////                            Eigen::MatrixXi F(f_ids.size(), 3);
////                            for (int i = 0; i < f_ids.size(); i++) {
////                                for (int k = 0; k < 3; k++) {
////                                    V.row(i * 3 + k) = input_vertices[input_faces[f_ids[i]][k]];
////                                    C.row(i * 3 + k) << 0, 0, 255;
////                                }
////                                F.row(i) << i * 3, i * 3 + 1, i * 3 + 2;
////                            }
////                            igl::writeOFF("_covered_input_fs_" + std::to_string(t_id) + ".off", V, F, C);
////                        }
////                        {
////                            Eigen::MatrixXd V(3, 3), C(3, 3);
////                            Eigen::MatrixXi F(1, 3);
////                            for (int k = 0; k < 3; k++) {
////                                V.row(k) = mesh.tet_vertices[mesh.tets[t_id][(j + 1 + k) % 4]].pos;
////                                C.row(k) << 255, 0, 0;
////                            }
////                            F.row(0) << 0, 1, 2;
////                            igl::writeOFF("_covered_tet_f_" + std::to_string(t_id) + ".off", V, F, C);
////                        }
////                        pausee();
////                    }
//                    //fortest
//                    continue;
//                }
                //


//                int t = get_t(tp1_3d, tp2_3d, tp3_3d);
//                std::array<Vector2, 3> tps_2d = {{to_2d(tp1_3d, t), to_2d(tp2_3d, t), to_2d(tp3_3d, t)}};
//                std::array<Vector2, 4> cs = {{(tps_2d[0] + tps_2d[1] + tps_2d[2]) / 3}};
//                for (int k = 0; k < 3; k++)
//                    cs[k + 1] = (tps_2d[k] + cs[0]) / 3;
//                std::array<Vector2, 3> ps_2d;
//                bool is_all_out = true;
//                for (int f_id: f_ids) {
//                    if (!is_face_inserted[f_id])
//                        continue;
//
//                    ps_2d = {{to_2d(input_vertices[input_faces[f_id][0]], t),
//                                     to_2d(input_vertices[input_faces[f_id][1]], t),
//                                     to_2d(input_vertices[input_faces[f_id][2]], t)}};
//
//                    bool is_in = false;
//                    for (auto &c:cs) {
//                        if (is_on_bounded_side(ps_2d, c)) {
//                            is_in = true;
//                            break;
//                        }
//                    }
//                    if (is_in) {
//                        ff_id = f_id;
//                        is_all_out = false;
//                        break;
//                    }
//                }
//                if(is_all_out)
//                    continue;
//                //

#ifdef NEW_ENVELOPE
                if (tree.is_out_sf_envelope_exact({{tp1_3d, tp2_3d, tp3_3d}}))
                    continue;
                else
                    ff_id = track_surface_fs[t_id][j].front();
#else
                    double eps_2 = (mesh.params.eps + mesh.params.eps_simplification) / 2;
                    double dd = (mesh.params.dd + mesh.params.dd_simplification) / 2;
                    eps_2 *= eps_2;
    #ifdef STORE_SAMPLE_POINTS
                    std::vector<GEO::vec3> ps;
                    sample_triangle({{tp1_3d, tp2_3d, tp3_3d}}, ps, dd);
                    if (tree.is_out_sf_envelope(ps, eps_2))
    #else
                    GEO::index_t prev_facet = GEO::NO_FACET;
                    if(sample_triangle_and_check_is_out({{tp1_3d, tp2_3d, tp3_3d}}, dd, eps_2, tree, prev_facet))
    #endif
                        continue;
                    else
                        ff_id = track_surface_fs[t_id][j].front();
#endif

//                int t = get_t(tp1_3d, tp2_3d, tp3_3d);
//                std::array<Vector2, 3> tps_2d = {{to_2d(tp1_3d, t), to_2d(tp2_3d, t), to_2d(tp3_3d, t)}};
//                std::array<Vector2, 4> cs = {{(tps_2d[0] + tps_2d[1] + tps_2d[2]) / 3}};
//                for(int k=0;k<3;k++)
//                    cs[k+1] = (tps_2d[k] + cs[0]) / 3;
//
//                std::array<Vector2, 3> ps_2d;
//                for (int f_id: f_ids) {
//                    if (!is_face_inserted[f_id])
//                        continue;
//
//                    ps_2d = {{to_2d(input_vertices[input_faces[f_id][0]], t),
//                                     to_2d(input_vertices[input_faces[f_id][1]], t),
//                                     to_2d(input_vertices[input_faces[f_id][2]], t)}};
//
//                    bool is_in = false;
//                    for(auto& c:cs) {
//                        if (is_on_bounded_side(ps_2d, c)){
//                            is_in = true;
//                            break;
//                        }
//                    }
//                    if(is_in) {
//                        ff_id = f_id;
//                        break;
//                    }
//                }
//                if (ff_id < 0) {
//                    continue;
////                    std::vector<GEO::vec3> ps;
////                    sample_triangle({{tp1_3d, tp2_3d, tp3_3d}}, ps, mesh.params.dd);
////                    if (tree.is_out_sf_envelope(ps, mesh.params.eps_2))
////                        continue;
////                    else
////                        ff_id = track_surface_fs[t_id][j].front();
//                }
//                else {
//                    std::vector<GEO::vec3> ps;
//                    sample_triangle({{tp1_3d, tp2_3d, tp3_3d}}, ps, mesh.params.dd);
//                    if (tree.is_out_sf_envelope(ps, mesh.params.eps_2))
//                        continue;
//                }

                opp_t_id = get_opp_t_id(t_id, j, mesh);
                if (opp_t_id < 0) {
                    mesh.tets[t_id].is_surface_fs[j] = NOT_SURFACE;
                    continue;
                }
                k = get_local_f_id(opp_t_id, mesh.tets[t_id][(j + 1) % 4], mesh.tets[t_id][(j + 2) % 4],
                                   mesh.tets[t_id][(j + 3) % 4], mesh);
                is_visited[opp_t_id][k] = true;
            }

            myassert(ff_id>=0, "ff_id<0!!!");//fortest

            mesh.tets[t_id].surface_tags[j] = input_tags[ff_id];
            mesh.tets[opp_t_id].surface_tags[k] = input_tags[ff_id];

            auto &fv1 = input_vertices[input_faces[ff_id][0]];
            auto &fv2 = input_vertices[input_faces[ff_id][1]];
            auto &fv3 = input_vertices[input_faces[ff_id][2]];
            //
            int ori = Predicates::orient_3d(fv1, fv2, fv3, mesh.tet_vertices[mesh.tets[t_id][j]].pos);
            int opp_ori = Predicates::orient_3d(fv1, fv2, fv3, mesh.tet_vertices[mesh.tets[opp_t_id][k]].pos);
            //
            if (ori == Predicates::ORI_POSITIVE && opp_ori == Predicates::ORI_NEGATIVE
                || ori == Predicates::ORI_NEGATIVE && opp_ori == Predicates::ORI_POSITIVE) {
                mesh.tets[t_id].is_surface_fs[j] = ori;
                mesh.tets[opp_t_id].is_surface_fs[k] = opp_ori;
                continue;
            }
            //
            if (ori == Predicates::ORI_ZERO && opp_ori != Predicates::ORI_ZERO) {
                mesh.tets[t_id].is_surface_fs[j] = -opp_ori;
                mesh.tets[opp_t_id].is_surface_fs[k] = opp_ori;
                continue;
            }
            if (opp_ori == Predicates::ORI_ZERO && ori != Predicates::ORI_ZERO) {
                mesh.tets[t_id].is_surface_fs[j] = ori;
                mesh.tets[opp_t_id].is_surface_fs[k] = -ori;
                continue;
            }
            //
            if (ori == opp_ori) {
                Vector3 n = (fv2 - fv1).cross(fv3 - fv1);
                n.normalize();
                Scalar dist = n.dot(mesh.tet_vertices[mesh.tets[t_id][j]].pos - fv1);
                Scalar opp_dist = n.dot(mesh.tet_vertices[mesh.tets[opp_t_id][k]].pos - fv1);
                if (ori == Predicates::ORI_ZERO) {
//                    cout << "impossible!! " << dist << " " << opp_dist << endl;
//                    cout << (mesh.tet_vertices[mesh.tets[t_id][j]].pos == mesh.tet_vertices[mesh.tets[opp_t_id][k]].pos)
//                         << endl;
//                    cout << t_id << " " << j << " " << mesh.tets[t_id][j] << endl;
//                    cout << opp_t_id << " " << k << " " << mesh.tets[opp_t_id][k] << endl;
//                    cout << "n = " << n << endl;

                    auto &tv1 = mesh.tet_vertices[mesh.tets[t_id][(j + 1) % 4]].pos;
                    auto &tv2 = mesh.tet_vertices[mesh.tets[t_id][(j + 2) % 4]].pos;
                    auto &tv3 = mesh.tet_vertices[mesh.tets[t_id][(j + 3) % 4]].pos;
                    Vector3 nt;
                    if (Predicates::orient_3d(tv1, tv2, tv3, mesh.tet_vertices[mesh.tets[t_id][j]].pos)
                        == Predicates::ORI_POSITIVE)
                        nt = (tv2 - tv1).cross(tv3 - tv1);
                    else
                        nt = (tv3 - tv1).cross(tv2 - tv1);
                    nt.normalize();
                    if (n.dot(nt) > 0) {
                        mesh.tets[opp_t_id].is_surface_fs[k] = 1;
                        mesh.tets[t_id].is_surface_fs[j] = -1;
                    } else {
                        mesh.tets[opp_t_id].is_surface_fs[k] = -1;
                        mesh.tets[t_id].is_surface_fs[j] = 1;
                    }
                } else {
                    if (dist < opp_dist) {
                        mesh.tets[opp_t_id].is_surface_fs[k] = opp_ori;
                        mesh.tets[t_id].is_surface_fs[j] = -ori;
                    } else {
                        mesh.tets[opp_t_id].is_surface_fs[k] = -opp_ori;
                        mesh.tets[t_id].is_surface_fs[j] = ori;
                    }
                }
            }
            //
        }
    }

    cout<<"known_surface_fs.size = "<<known_surface_fs.size()<<endl;
    cout<<"known_not_surface_fs.size = "<<known_not_surface_fs.size()<<endl;
    if(known_surface_fs.empty() && known_not_surface_fs.empty())
        return;

//    Eigen::MatrixXd V(known_surface_fs.size()*3, 3);
//    Eigen::MatrixXi F(known_surface_fs.size(), 3);
//    for(int i=0;i<known_surface_fs.size();i++) {
//        for (int j = 0; j < 3; j++)
//            V.row(i * 3 + j) = mesh.tet_vertices[known_surface_fs[i][j]].pos;
//        F.row(i) << i * 3, i * 3 + 1, i * 3 + 2;
//    }
//    igl::writeOFF("known_surface_fs.off",V, F);
//    pausee("writing known_surface_fs.off");


    //b_edges
    std::vector<std::array<int, 2>> tmp_edges;
    for (const auto &f: known_surface_fs) {
        for (int j = 0; j < 3; j++) {
            if (f[j] < f[(j + 1) % 3])
                tmp_edges.push_back({{f[j], f[(j + 1) % 3]}});
            else
                tmp_edges.push_back({{f[(j + 1) % 3], f[j]}});
        }
    }
    for (const auto &f: known_not_surface_fs) {
        for (int j = 0; j < 3; j++) {
            if (f[j] < f[(j + 1) % 3])
                tmp_edges.push_back({{f[j], f[(j + 1) % 3]}});
            else
                tmp_edges.push_back({{f[(j + 1) % 3], f[j]}});
        }
    }
    vector_unique(tmp_edges);

    for (auto &e: tmp_edges) {
        int cnt = 0;
        std::vector<int> n_t_ids;
        set_intersection(mesh.tet_vertices[e[0]].conn_tets, mesh.tet_vertices[e[1]].conn_tets, n_t_ids);
        for (int t_id: n_t_ids) {
            for (int j = 0; j < 4; j++) {
                if (mesh.tets[t_id][j] != e[0] && mesh.tets[t_id][j] != e[1]) {
                    if (mesh.tets[t_id].is_surface_fs[j] != NOT_SURFACE) {
                        cnt++;
                        if (cnt > 2)
                            break;
                    }
                }
                if (cnt > 2)
                    break;
            }
        }
//        cout<<"cnt = "<<cnt<<endl;
        if (cnt == 2) {
            b_edges.push_back(e);
            mesh.tet_vertices[e[0]].is_on_boundary = true;
//            cout<<"b_edges.push_back(e);"<<endl;
        }
    }
}

bool floatTetWild::is_uninserted_face_covered(int uninserted_f_id, const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
        const std::vector<int>& cut_t_ids, Mesh &mesh){

    std::array<Vector3, 3> vs = {{input_vertices[input_faces[uninserted_f_id][0]],
                                         input_vertices[input_faces[uninserted_f_id][1]],
                                         input_vertices[input_faces[uninserted_f_id][2]]}};
    std::vector<GEO::vec3> ps;
    sample_triangle(vs, ps, mesh.params.dd);

    std::vector<int> n_t_ids;
    for(int t_id: cut_t_ids) {
        for (int j = 0; j < 4; j++)
            n_t_ids.insert(n_t_ids.end(), mesh.tet_vertices[mesh.tets[t_id][j]].conn_tets.begin(),
                           mesh.tet_vertices[mesh.tets[t_id][j]].conn_tets.end());
    }
    vector_unique(n_t_ids);

    std::vector<std::array<int, 3>> faces;
    for(int t_id: n_t_ids) {
        for (int j = 0; j < 4; j++){
            if(mesh.tets[t_id].is_surface_fs[j] != NOT_SURFACE) {
                faces.push_back(
                        {{mesh.tets[t_id][(j + 1) % 4], mesh.tets[t_id][(j + 2) % 4], mesh.tets[t_id][(j + 3) % 4]}});
                std::sort(faces.back().begin(), faces.back().end());
            }
        }
    }
    vector_unique(faces);

    for(auto& p: ps){
        bool is_valid = false;
        for(auto& f: faces) {
            double dis_2 = GEO::Geom::point_triangle_squared_distance(p, to_geo_p(mesh.tet_vertices[f[0]].pos),
                                                                      to_geo_p(mesh.tet_vertices[f[1]].pos),
                                                                      to_geo_p(mesh.tet_vertices[f[2]].pos));
            if (dis_2 < mesh.params.eps_2) {
                is_valid = true;
                break;
            }
        }
        if(!is_valid)
            return false;
    }

    cout<<"covered!!!!!!!"<<endl;
    return true;
}

int floatTetWild::get_opp_t_id(int t_id, int j, const Mesh &mesh){
    std::vector<int> tmp;
    set_intersection(mesh.tet_vertices[mesh.tets[t_id][(j + 1) % 4]].conn_tets,
                     mesh.tet_vertices[mesh.tets[t_id][(j + 2) % 4]].conn_tets,
                     mesh.tet_vertices[mesh.tets[t_id][(j + 3) % 4]].conn_tets,
                     tmp);
//    //fortest
//    if(tmp.size() != 1 && tmp.size() != 2) {
//        cout << "tmp.size() = " << tmp.size() << endl;
//        cout << "t_id = " << t_id << ", j = " << j << endl;
//        vector_print(mesh.tet_vertices[mesh.tets[t_id][(j + 1) % 4]].conn_tets);
//        vector_print(mesh.tet_vertices[mesh.tets[t_id][(j + 2) % 4]].conn_tets);
//        vector_print(mesh.tet_vertices[mesh.tets[t_id][(j + 3) % 4]].conn_tets);
//        vector_print(tmp);
//        for(int i: tmp){
//            mesh.tets[i].print();
//        }
//        pausee();
//    }
//    //fortest
    if (tmp.size() == 2)
        return tmp[0] == t_id ? tmp[1] : tmp[0];
    else
        return -1;
}

void floatTetWild::myassert(bool b, const std::string& s) {
    if (b == false) {
        cout << "myassert fail: " << s << endl;
        pausee();
    }
}

void floatTetWild::check_track_surface_fs(Mesh &mesh, std::vector<std::array<std::vector<int>, 4>> &track_surface_fs,
                                          const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
                                          const std::vector<int> &sorted_f_ids) {
    return;

//    //check connection
//    std::vector<std::vector<int>> conn_tets(mesh.tet_vertices.size());
//    for (int i = 0; i < mesh.tets.size(); i++) {
//        for (int j = 0; j < 4; j++)
//            conn_tets[mesh.tets[i][j]].push_back(i);
//    }
//    for (int i = 0; i < mesh.tet_vertices.size(); i++) {
//        std::sort(mesh.tet_vertices[i].conn_tets.begin(), mesh.tet_vertices[i].conn_tets.end());
//        if (mesh.tet_vertices[i].conn_tets != conn_tets[i]) {
//            cout << "mesh.tet_vertices[i].conn_tets!=conn_tets[i]" << endl;
//            pausee();
//        }
//    }
//    cout<<"check 1 done"<<endl;
//
//    for (int i = 0; i < mesh.tets.size(); i++) {
//        for (int j = 0; j < 4; j++) {
//            int opp_t_id = get_opp_t_id(i, j, mesh);
//        }
//    }
//    cout<<"check 2 done"<<endl;

//    for (int t_id = 0; t_id < track_surface_fs.size(); t_id++) {
//        for (int j = 0; j < 4; j++) {
//            std::sort(track_surface_fs[t_id][j].begin(), track_surface_fs[t_id][j].end());
//            int opp_t_id = get_opp_t_id(t_id, j, mesh);
//            if (opp_t_id < 0) {
//                if (!track_surface_fs[t_id][j].empty()) {
//                    cout << "bbox face but !track_surface_fs[t_id][j].empty()" << endl;
//                    pausee();
//                }
//                continue;
//            }
//            int k = 0;
//            for (k = 0; k < 4; k++) {
//                if (mesh.tets[opp_t_id][k] != mesh.tets[t_id][(j + 1) % 4]
//                    && mesh.tets[opp_t_id][k] != mesh.tets[t_id][(j + 2) % 4]
//                    && mesh.tets[opp_t_id][k] != mesh.tets[t_id][(j + 3) % 4])
//                    break;
//            }
//            std::sort(track_surface_fs[opp_t_id][k].begin(), track_surface_fs[opp_t_id][k].end());
//            if (track_surface_fs[t_id][j] != track_surface_fs[opp_t_id][k]) {
//                cout << "track_surface_fs[t_id][j]!=track_surface_fs[opp_t_id][k]" << endl;
//                cout << "t " << t_id << ": ";
//                vector_print(track_surface_fs[t_id][j]);
//                cout << "opp_t " << opp_t_id << ": ";
//                vector_print(track_surface_fs[opp_t_id][k]);
//                pausee();
//            }
//        }
//    }
//    cout<<"check 3 done"<<endl;

    std::vector<std::vector<std::pair<int, int>>> covered_fs_infos(input_faces.size());
    for (int i = 0; i < track_surface_fs.size(); i++) {
        for (int j = 0; j < 4; j++) {
            for (int f_id: track_surface_fs[i][j])
                covered_fs_infos[f_id].push_back(std::make_pair(i, j));
        }
    }

    //check is covered
    for (int i = 0; i < sorted_f_ids.size(); i++) {
        int f_id = sorted_f_ids[i];
        if (covered_fs_infos[f_id].empty())//not inserted
            continue;
        std::array<Vector3, 3> f_vs = {{input_vertices[input_faces[f_id][0]],
                                               input_vertices[input_faces[f_id][1]],
                                               input_vertices[input_faces[f_id][2]]}};
        int t = get_t(f_vs[0], f_vs[1], f_vs[2]);
        std::array<Vector2, 3> f_vs_2d = {{to_2d(input_vertices[input_faces[f_id][0]], t),
                                                  to_2d(input_vertices[input_faces[f_id][1]], t),
                                                  to_2d(input_vertices[input_faces[f_id][2]], t)}};
        //
//        cout<<"covered_fs_infos[f_id].size = "<<covered_fs_infos[f_id].size()<<"->";
//        cout<<"f_id = "<<f_id<<endl;
        for (int n = 0; n < covered_fs_infos[f_id].size(); n++) {
            int t_id = covered_fs_infos[f_id][n].first;
            int j = covered_fs_infos[f_id][n].second;
            std::array<Vector3, 3> t_vs = {{mesh.tet_vertices[mesh.tets[t_id][(j + 1) % 4]].pos,
                                                   mesh.tet_vertices[mesh.tets[t_id][(j + 2) % 4]].pos,
                                                   mesh.tet_vertices[mesh.tets[t_id][(j + 3) % 4]].pos}};
            std::array<Vector2, 3> t_vs_2d = {{to_2d(mesh.tet_vertices[mesh.tets[t_id][(j + 1) % 4]].pos, t),
                                                      to_2d(mesh.tet_vertices[mesh.tets[t_id][(j + 2) % 4]].pos, t),
                                                      to_2d(mesh.tet_vertices[mesh.tets[t_id][(j + 3) % 4]].pos, t)}};
//            int k = -1;
//            for (int r = 0; r < 3; r++) {
//                if(is_p_inside_tri_2d(t_vs_2d[r], f_vs_2d))
//                    continue;
//                double dist_2 = GEO::Geom::point_triangle_squared_distance(
//                        GEO::vec3(t_vs[r][0], t_vs[r][1], t_vs[r][2]),
//                        GEO::vec3(f_vs[0][0], f_vs[0][1], f_vs[0][2]),
//                        GEO::vec3(f_vs[1][0], f_vs[1][1], f_vs[1][2]),
//                        GEO::vec3(f_vs[2][0], f_vs[2][1], f_vs[2][2]));
//                if (dist_2 > mesh.params.eps_2) {
////                    cout << "r = " << r << endl;
////                    cout << "dist_2 = " << dist_2 << endl;
//                    k = r;
//                    break;
//                }
//            }
//            if (k < 0)
//                continue;
//            bool is_out = true;
//            for (int j = 0; j < 3; j++) {
//                int ori_out = Predicates::orient_2d(f_vs_2d[j], f_vs_2d[(j + 1) % 3], t_vs_2d[k]);
////                cout << "ori_out = " << ori_out << endl;
//                is_out = true;
//                for (int r = 0; r < 3; r++) {
//                    int ori_r = Predicates::orient_2d(f_vs_2d[j], f_vs_2d[(j + 1) % 3], t_vs_2d[r]);
////                    cout << "ori_r = " << ori_r << endl;
//                    if (ori_r == ori_out || ori_r == Predicates::ORI_ZERO)
//                        continue;
//                    double dist_2 = p_line_squared_dist_3d(t_vs[r], f_vs[j], f_vs[(j + 1) % 3]);
////                    cout << "dist_2 = " << dist_2 << "; " << mesh.params.eps_2_coplanar << endl;
//                    if (dist_2 <= mesh.params.eps_2_coplanar)
//                        continue;
//                    is_out = false;
//                    break;
//                }
//                if (is_out) {
//                    break;
//                }
//            }
            Vector2 c = (t_vs_2d[0]+t_vs_2d[1]+t_vs_2d[2])/3;
            if(!is_p_inside_tri_2d(c, f_vs_2d)){
                bool is_crossed = false;
                for (int k = 0; k < 3; k++) {//edges of input face
                    std::array<int, 3> oris;
                    for (int r = 0; r < 3; r++) {//vertices of tet face
                        oris[r] = Predicates::orient_2d(f_vs_2d[k], f_vs_2d[(k + 1) % 3], t_vs_2d[r]);
                        if (oris[r] == Predicates::ORI_ZERO)
                            continue;
                        double dist_2 = p_line_squared_dist_3d(t_vs[r], f_vs[k], f_vs[(k + 1) % 3]);
                        if (dist_2 <= mesh.params.eps_2_coplanar)
                            oris[r] = Predicates::ORI_ZERO;
                    }
                    for (int r = 0; r < 3; r++) {
                        if (is_crossing(oris[r], oris[(r + 1) % 3])) {
                            is_crossed = true;
                            break;
                        }
                    }
                    if (is_crossed)
                        break;
                }
                if(!is_crossed){
//                    cout << "is out!!" << endl;
                    vector_erase(track_surface_fs[t_id][j], f_id);
                    covered_fs_infos[f_id].erase(covered_fs_infos[f_id].begin() + n);
                    n--;
                }
            }
        }
        if (covered_fs_infos[f_id].empty()) {
            cout << "covered_fs_infos[f_id].empty()!!!" << endl;
            pausee();
            continue;
        }
//        cout<<covered_fs_infos[f_id].size()<<endl;
        //
        std::vector<GEO::vec3> ps;
        sample_triangle(f_vs, ps, mesh.params.dd);
        //
        double eps = mesh.params.eps_coplanar * mesh.params.eps_coplanar;
        //
        if (ps.size() < 4) {
            cout << "ps.size() < 4!!" << endl;
            ps.clear();
            sample_triangle(f_vs, ps, mesh.params.dd / 4);
        }
        //
        bool is_inside = false;
        for (const auto &p: ps) {
            is_inside = false;
            for (const auto &info: covered_fs_infos[f_id]) {
                int t_id = info.first;
                int j = info.second;
                std::array<int, 3> f = {{mesh.tets[t_id][(j + 1) % 4], mesh.tets[t_id][(j + 2) % 4],
                                                mesh.tets[t_id][(j + 3) % 4]}};
                std::array<GEO::vec3, 3> vs;
                for (int k = 0; k < 3; k++) {
                    vs[k] = GEO::vec3(mesh.tet_vertices[f[k]].pos[0], mesh.tet_vertices[f[k]].pos[1],
                                      mesh.tet_vertices[f[k]].pos[2]);
                }
                double dis_2 = GEO::Geom::point_triangle_squared_distance(p, vs[0], vs[1], vs[2]);
                if (dis_2 <= eps) {
                    is_inside = true;
                    break;
                }
            }
            if (!is_inside) {
                cout << "check is_input_face_covered fail!!" << endl;
                cout << "f_id = " << f_id << endl;
                cout << "i = " << i << endl;
                pausee();
                break;
            }
        }
        //
//        bool is_inside = true;//fortest
        if (is_inside/* && i < sorted_f_ids.size() * 0.98*/)
//        if (f_id != 3083)
            continue;
        cout << i << " f_id = " << f_id << endl;
        //output input/tet triangles in different colors
        {
            auto &infos = covered_fs_infos[f_id];
            Eigen::MatrixXd V(infos.size() * 3, 3), C(infos.size() * 3, 3);
            Eigen::MatrixXi F(infos.size(), 3);
            for (int i = 0; i < infos.size(); i++) {
                int t_id = infos[i].first;
                int k = infos[i].second;
                for (int j = 0; j < 3; j++) {
                    V.row(i * 3 + j) = mesh.tet_vertices[mesh.tets[t_id][(k + j + 1) % 4]].pos;
                    C.row(i * 3 + j) << 0, 0, 255;
                }
                F.row(i) << i * 3, i * 3 + 1, i * 3 + 2;
            }
            igl::writeOFF(std::to_string(i) + "_covered_tet_fs_" + std::to_string(f_id) + "_"
                          + std::to_string(is_inside) + ".off", V, F, C);
        }
        {
            Eigen::MatrixXd V(3, 3), C(3, 3);
            Eigen::MatrixXi F(1, 3);
            for (int j = 0; j < 3; j++) {
                V.row(j) = f_vs[j];
                C.row(j) << 255, 0, 0;
            }
            F.row(0) << 0, 1, 2;
            igl::writeOFF(std::to_string(i) + "_covered_input_f_" + std::to_string(f_id) + "_"
                          + std::to_string(is_inside) + ".off", V, F, C);
        }
//        pausee();
    }
}

int floatTetWild::orient_rational(const Vector3_r& p1, const Vector3_r& p2, const Vector3_r& p3, const Vector3_r& p){
    auto nv = (p2-p1).cross(p3-p1);
    triwild::Rational res = nv.dot(p-p1);
    if(res == 0)
        return Predicates::ORI_ZERO;
    if(res < 0)
        return Predicates::ORI_POSITIVE;
    else
        return Predicates::ORI_NEGATIVE;
}
