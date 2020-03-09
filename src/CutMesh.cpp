// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

//
// Created by Yixin Hu on 9/12/19.
//

#include <floattetwild/CutMesh.h>
#include <floattetwild/TriangleInsertion.h>
#include <floattetwild/LocalOperations.h>
#include <floattetwild/Predicates.hpp>
#include <floattetwild/intersections.h>
#include <floattetwild/Logger.hpp>

#include <igl/Timer.h>

double time_cut_mesh11 = 0;
double time_cut_mesh12 = 0;
double time_cut_mesh13 = 0;
double time_cut_mesh14 = 0;
double time_get_intersecting_edges_and_points1 = 0;
double time_get_intersecting_edges_and_points2 = 0;
double time_get_intersecting_edges_and_points3 = 0;
double time_get_intersecting_edges_and_points4 = 0;
//igl::Timer timer;

void floatTetWild::print_times1(){
//    logger().debug("\t\t\t- time_cut_mesh11 = {}s", time_cut_mesh11);
//    logger().debug("\t\t\t- time_cut_mesh12 = {}s", time_cut_mesh12);
//    logger().debug("\t\t\t- time_cut_mesh13 = {}s", time_cut_mesh13);
//    logger().debug("\t\t\t- time_cut_mesh14 = {}s", time_cut_mesh14);
    logger().debug("\t\t- time_get_intersecting_edges_and_points1 = {}s", time_get_intersecting_edges_and_points1);
    logger().debug("\t\t- time_get_intersecting_edges_and_points2 = {}s", time_get_intersecting_edges_and_points2);
    logger().debug("\t\t- time_get_intersecting_edges_and_points3 = {}s", time_get_intersecting_edges_and_points3);
    logger().debug("\t\t- time_get_intersecting_edges_and_points4 = {}s", time_get_intersecting_edges_and_points4);
}

void floatTetWild::CutMesh::construct(const std::vector<int>& cut_t_ids) {
    v_ids.reserve(cut_t_ids.size() * 4);
    for (int t_id:cut_t_ids) {
        for (int j = 0; j < 4; j++)
            v_ids.push_back(mesh.tets[t_id][j]);
    }
    vector_unique(v_ids);

    for (int i = 0; i < v_ids.size(); i++)
        map_v_ids[v_ids[i]] = i;

    tets.resize(cut_t_ids.size());
    for (int i = 0; i < cut_t_ids.size(); i++) {
        tets[i] = {{map_v_ids[mesh.tets[cut_t_ids[i]][0]],
                           map_v_ids[mesh.tets[cut_t_ids[i]][1]],
                           map_v_ids[mesh.tets[cut_t_ids[i]][2]],
                           map_v_ids[mesh.tets[cut_t_ids[i]][3]]}};
    }

//    std::vector<std::vector<int>> conn_tets(v_ids.size());
//    for (int i = 0; i < tets.size(); i++) {
//        for (int j = 0; j < 4; j++)
//            conn_tets[tets[i][j]].push_back(i);
//    }
//
//    timer.start();
//    opp_t_ids.resize(tets.size(), {{-1, -1, -1, -1}});
//    time_cut_mesh12 += timer.getElapsedTime();
//
//    for (int i = 0; i < tets.size(); i++) {//todo: construct conn_tets/opp_t_ids only when expension is required.
//        for (int j = 0; j < 4; j++) {
//            if (opp_t_ids[i][j] >= 0)
//                continue;
//
//            timer.start();
//            std::vector<int> n_t_ids;
//            set_intersection_sorted(conn_tets[tets[i][(j + 1) % 4]], conn_tets[tets[i][(j + 2) % 4]],
//                                    conn_tets[tets[i][(j + 3) % 4]], n_t_ids);
//
//            time_cut_mesh13 += timer.getElapsedTime();
//            assert(!n_t_ids.empty());
//            if (n_t_ids.size() < 2)
//                continue;
//
//            timer.start();
//            int n_t_id = n_t_ids[0] == i ? n_t_ids[1] : n_t_ids[0];
//            opp_t_ids[i][j] = n_t_id;
//            for (int k = 0; k < 4; k++) {
//                if (tets[n_t_id][k] != tets[i][(j + 1) % 4] && tets[n_t_id][k] != tets[i][(j + 2) % 4]
//                    && tets[n_t_id][k] != tets[i][(j + 3) % 4]) {
//                    opp_t_ids[n_t_id][k] = i;
//                    break;
//                }
//            }
//            time_cut_mesh14 += timer.getElapsedTime();
//        }
//    }
}

bool floatTetWild::CutMesh::snap_to_plane() {
    bool snapped = false;
    to_plane_dists.resize(map_v_ids.size());
    is_snapped.resize(map_v_ids.size());
    for (auto &v:map_v_ids) {
        int v_id = v.first;
        int lv_id = v.second;

        int ori = Predicates::orient_3d(p_vs[0], p_vs[1], p_vs[2], mesh.tet_vertices[v_id].pos);
        if (ori == Predicates::ORI_ZERO) {
            to_plane_dists[lv_id] = 0;
            continue;
        }
        to_plane_dists[lv_id] = get_to_plane_dist(mesh.tet_vertices[v_id].pos);
        if (ori == Predicates::ORI_POSITIVE && to_plane_dists[lv_id] > 0
            || ori == Predicates::ORI_NEGATIVE && to_plane_dists[lv_id] < 0){
//            cout<<"reverted!!! "<<to_plane_dists[lv_id]<<endl;
            to_plane_dists[lv_id] = -to_plane_dists[lv_id];
        }

        if (std::fabs(to_plane_dists[lv_id]) < mesh.params.eps_coplanar) {
            is_snapped[lv_id] = true;
            snapped = true;
        }
    }

    revert_totally_snapped_tets(0, tets.size());

    return snapped;
}

void floatTetWild::CutMesh::expand(std::vector<int>& cut_t_ids) {
//    timer.start();
    std::vector<std::vector<int>> conn_tets(v_ids.size());
    for (int i = 0; i < tets.size(); i++) {
        for (int j = 0; j < 4; j++)
            conn_tets[tets[i][j]].push_back(i);
    }

    std::vector<std::array<int, 4>> opp_t_ids(tets.size(), {{-1, -1, -1, -1}});
    for (int i = 0; i < tets.size(); i++) {
        for (int j = 0; j < 4; j++) {
            if (opp_t_ids[i][j] >= 0)
                continue;

            std::vector<int> n_t_ids;
            set_intersection_sorted(conn_tets[tets[i][(j + 1) % 4]], conn_tets[tets[i][(j + 2) % 4]],
                                    conn_tets[tets[i][(j + 3) % 4]], n_t_ids);

            assert(!n_t_ids.empty());
            if (n_t_ids.size() < 2)
                continue;

            int n_t_id = n_t_ids[0] == i ? n_t_ids[1] : n_t_ids[0];
            opp_t_ids[i][j] = n_t_id;
            for (int k = 0; k < 4; k++) {
                if (tets[n_t_id][k] != tets[i][(j + 1) % 4] && tets[n_t_id][k] != tets[i][(j + 2) % 4]
                    && tets[n_t_id][k] != tets[i][(j + 3) % 4]) {
                    opp_t_ids[n_t_id][k] = i;
                    break;
                }
            }
        }
    }
//    time_cut_mesh12 += timer.getElapsedTime();

    while (true) {
//        timer.start();
        std::vector<std::array<int, 5>> new_opp_t_ids;
        for (int i = 0; i < opp_t_ids.size(); i++) {
            for (int j = 0; j < 4; j++) {
                if (opp_t_ids[i][j] >= 0)
                    continue;
                if (!is_snapped[tets[i][(j + 1) % 4]] && !is_snapped[tets[i][(j + 2) % 4]]
                    && !is_snapped[tets[i][(j + 3) % 4]])
                    continue;

                int n_gt_id = get_opp_t_id(cut_t_ids[i], j, mesh);
                if (n_gt_id < 0)
                    continue;
                new_opp_t_ids.push_back({{-1, -1, -1, -1, n_gt_id}});
                for (int k = 0; k < 4; k++) {
                    if (mesh.tets[n_gt_id][k] != v_ids[tets[i][(j + 1) % 4]]
                        && mesh.tets[n_gt_id][k] != v_ids[tets[i][(j + 2) % 4]]
                        && mesh.tets[n_gt_id][k] != v_ids[tets[i][(j + 3) % 4]]) {
                        new_opp_t_ids.back()[k] = i;
                        break;
                    }
                }
            }
        }
//        time_cut_mesh13 += timer.getElapsedTime();
        if(new_opp_t_ids.empty())
            return;

//        timer.start();
        std::sort(new_opp_t_ids.begin(), new_opp_t_ids.end(),
                  [&](const std::array<int, 5> &a, const std::array<int, 5> &b) {
                      return a.back() < b.back();
                  });
        for (int i = 0; i < new_opp_t_ids.size() - 1; i++) {
            if (new_opp_t_ids[i].back() == new_opp_t_ids[i + 1].back()) {
                for (int j = 0; j < 4; j++) {
                    if (new_opp_t_ids[i][j] >= 0) {
                        new_opp_t_ids[i + 1][j] = new_opp_t_ids[i][j];
                        break;
                    }
                }
                new_opp_t_ids.erase(new_opp_t_ids.begin() + i);
                i--;
            }
        }

        const int old_tets_size = tets.size();
        for (int i = 0; i < new_opp_t_ids.size(); i++) {
            ///
            int cnt_on = 0;
            for (int j = 0; j < 4; j++) {
                int v_id = mesh.tets[new_opp_t_ids[i].back()][j];
                if (map_v_ids.find(v_id) != map_v_ids.end() && is_v_on_plane(map_v_ids[v_id]))
                    cnt_on++;
            }
            if (cnt_on != 3) {
                int cnt_pos = 0;
                int cnt_neg = 0;
                for (int j = 0; j < 4; j++) {
                    int ori = Predicates::orient_3d(p_vs[0], p_vs[1], p_vs[2],
                                                    mesh.tet_vertices[mesh.tets[new_opp_t_ids[i].back()][j]].pos);
                    if (ori == Predicates::ORI_POSITIVE)
                        cnt_pos++;
                    else if (ori == Predicates::ORI_NEGATIVE)
                        cnt_neg++;
                }
                if (cnt_neg == 0 || cnt_pos == 0)
                    continue;
            }

            ///
            cut_t_ids.push_back(new_opp_t_ids[i].back());

            int t_id = tets.size();
            tets.emplace_back();
            auto &t = tets.back();
            const int old_v_ids_size = v_ids.size();
            for (int j = 0; j < 4; j++) {
                int v_id = mesh.tets[new_opp_t_ids[i].back()][j];
                int lv_id;
                if (map_v_ids.find(v_id) == map_v_ids.end()) {
                    v_ids.push_back(v_id);
                    lv_id = v_ids.size() - 1;
                    map_v_ids[v_id] = lv_id;
//                    to_plane_dists.push_back(get_to_plane_dist(mesh.tet_vertices[v_id].pos));
                    double dist = get_to_plane_dist(mesh.tet_vertices[v_id].pos);
                    int ori = Predicates::orient_3d(p_vs[0], p_vs[1], p_vs[2], mesh.tet_vertices[v_id].pos);
                    if ((ori == Predicates::ORI_NEGATIVE && dist < 0)//todo: change get_to_plane_dist return value sign
                        || (ori == Predicates::ORI_POSITIVE && dist > 0))
                        dist = -dist;
                    else if (ori == Predicates::ORI_ZERO)
                        dist = 0;
                    to_plane_dists.push_back(dist);

                    if (ori != Predicates::ORI_ZERO && std::fabs(to_plane_dists[lv_id]) < mesh.params.eps_coplanar)
                        is_snapped.push_back(true);
                    else
                        is_snapped.push_back(false);
                    conn_tets.emplace_back();
                } else
                    lv_id = map_v_ids[v_id];
                t[j] = map_v_ids[v_id];
                conn_tets[lv_id].push_back(t_id);
            }

            opp_t_ids.push_back({{new_opp_t_ids[i][0], new_opp_t_ids[i][1], new_opp_t_ids[i][2], new_opp_t_ids[i][3]}});
            for (int j = 0; j < 4; j++) {
                if (opp_t_ids.back()[j] < 0) {
                    if (t[(j + 1) % 4] < old_v_ids_size && t[(j + 2) % 4] < old_v_ids_size
                        && t[(j + 3) % 4] < old_v_ids_size) {
                        std::vector<int> tmp;
                        set_intersection(conn_tets[t[(j + 1) % 4]], conn_tets[t[(j + 2) % 4]],
                                         conn_tets[t[(j + 3) % 4]], tmp);
                        if(tmp.size() == 1)//一个三角形刚好填了一个v字的缺口，所有点都在v_ids中，但三角形不在cut_mesh中
                            continue;
                        opp_t_ids.back()[j] = tmp[0] == t_id ? tmp[1] : tmp[0];
                    } else
                        continue;
                }
                int opp_t_id = opp_t_ids.back()[j];
                for (int k = 0; k < 4; k++) {
                    if (tets[opp_t_id][k] != t[(j + 1) % 4] && tets[opp_t_id][k] != t[(j + 2) % 4]
                        && tets[opp_t_id][k] != t[(j + 3) % 4]) {
                        opp_t_ids[opp_t_id][k] = t_id;
                        break;
                    }
                }
            }
        }
//        time_cut_mesh14 += timer.getElapsedTime();
        if (old_tets_size == tets.size())
            break;

        revert_totally_snapped_tets(old_tets_size, tets.size());
    }

//    //fortest
//    for(int i=0;i<v_ids.size();i++){
//        int v_id = v_ids[i];
//        double dis = get_to_plane_dist(mesh.tet_vertices[v_id].pos);
//        if(to_plane_dists[map_v_ids[v_id]] != dis){
//            cout<<"expand.. to_plane_dists[map_v_ids[v_id]] != dis"<<endl;
//            cout<<to_plane_dists[map_v_ids[v_id]] <<" "<< dis<<endl;
//            cout<<i<<endl;
//            pausee();
//        }
//    }
//    //fortest
}

void floatTetWild::CutMesh::expand_new(std::vector<int> &cut_t_ids) {
    const int t = get_t(p_vs[0], p_vs[1], p_vs[2]);
    const std::array<Vector2, 3> tri_2d = {{to_2d(p_vs[0], t), to_2d(p_vs[1], t), to_2d(p_vs[2], t)}};

    std::vector<bool> is_in_cutmesh(mesh.tets.size(), false);
    for (int t_id:cut_t_ids)
        is_in_cutmesh[t_id] = true;

    int cnt_loop=0;
    std::vector<bool> is_interior(v_ids.size(), false);
    while (true) {
        cnt_loop++;
        /////
        std::vector<bool> is_visited(mesh.tets.size(), false);
        for (int t_id:cut_t_ids)
            is_visited[t_id] = true;

        /////
        int old_cut_t_ids = cut_t_ids.size();
        for (auto m: map_v_ids) {
            int gv_id = m.first;
            int lv_id = m.second;

            if (is_interior[lv_id])
                continue;
            if (!is_snapped[lv_id])
                continue;
//            if(is_projected[lv_id])
//                continue;

            bool is_in = true;
            for (int gt_id: mesh.tet_vertices[gv_id].conn_tets) {
                if (is_in_cutmesh[gt_id])
                    continue;
                is_in = false;

                if (is_visited[gt_id])
                    continue;
                is_visited[gt_id] = true;

                ///
                int cnt = 0;
                int cnt_on = 0;
                for (int j = 0; j < 4; j++) {
                    int tmp_gv_id = mesh.tets[gt_id][j];
                    if (map_v_ids.find(tmp_gv_id) != map_v_ids.end()) {
                        cnt++;
                        if (is_v_on_plane(map_v_ids[tmp_gv_id]))
                            cnt_on++;
                    }
                }
                if (cnt < 3)
                    continue;
//                if (cnt_on < 3) {
                    int cnt_pos = 0;
                    int cnt_neg = 0;
                    for (int j = 0; j < 4; j++) {
                        int ori = Predicates::orient_3d(p_vs[0], p_vs[1], p_vs[2],
                                                        mesh.tet_vertices[mesh.tets[gt_id][j]].pos);
                        if (ori == Predicates::ORI_POSITIVE)
                            cnt_pos++;
                        else if (ori == Predicates::ORI_NEGATIVE)
                            cnt_neg++;
                    }
                    if (cnt_neg == 0 || cnt_pos == 0)
                        continue;

                    bool is_overlapped = false;
                    std::array<Vector2, 4> tet_2d;
                    for (int j = 0; j < 4; j++) {
                        Scalar dist = get_to_plane_dist(mesh.tet_vertices[mesh.tets[gt_id][j]].pos);
                        Vector3 proj_p = mesh.tet_vertices[mesh.tets[gt_id][j]].pos - dist * p_n;
                        tet_2d[j] = to_2d(proj_p, t);
                    }
                    for(int j=0;j<4;j++) {
                        if (is_tri_tri_cutted_2d({{tet_2d[(j + 1) % 4], tet_2d[(j + 2) % 4], tet_2d[(j + 3) % 4]}},
                                                 tri_2d)) {
                            is_overlapped = true;
                            break;
                        }
                    }
                    if(!is_overlapped)
                        continue;
//                }

                ///
                cut_t_ids.push_back(gt_id);
                is_in_cutmesh[gt_id] = true;

                ///
                tets.emplace_back();
                auto &t = tets.back();
                for (int j = 0; j < 4; j++) {
                    int new_gv_id = mesh.tets[gt_id][j];
                    int new_lv_id;
                    if (map_v_ids.find(new_gv_id) == map_v_ids.end()) {
                        //
                        v_ids.push_back(new_gv_id);
                        is_interior.push_back(false);
                        new_lv_id = v_ids.size() - 1;
                        map_v_ids[new_gv_id] = new_lv_id;
                        //
                        double dist = get_to_plane_dist(mesh.tet_vertices[new_gv_id].pos);
                        int ori = Predicates::orient_3d(p_vs[0], p_vs[1], p_vs[2], mesh.tet_vertices[new_gv_id].pos);
                        if ((ori == Predicates::ORI_NEGATIVE && dist < 0)
                            || (ori == Predicates::ORI_POSITIVE && dist > 0))
                            dist = -dist;
                        else if (ori == Predicates::ORI_ZERO)
                            dist = 0;
                        to_plane_dists.push_back(dist);
                        //
                        if (ori != Predicates::ORI_ZERO &&
                            std::fabs(to_plane_dists[new_lv_id]) < mesh.params.eps_coplanar)
                            is_snapped.push_back(true);
                        else
                            is_snapped.push_back(false);
                        is_projected.push_back(false);
                    } else
                        new_lv_id = map_v_ids[new_gv_id];
                    t[j] = new_lv_id;
                }
            }
            if (is_in)
                is_interior[lv_id] = true;
        }
        if (cut_t_ids.size() == old_cut_t_ids)
            break;
    }
    revert_totally_snapped_tets(0, tets.size());

//    cout<<"("<<cnt_loop<<")";

//    //fortest
//    std::vector <std::vector<int>> conn_tets(v_ids.size());
//    for (int i = 0; i < tets.size(); i++) {
//        for (int j = 0; j < 4; j++)
//            conn_tets[tets[i][j]].push_back(i);
//    }
//
//    std::vector <std::array<int, 4>> opp_t_ids(tets.size(), {{-1, -1, -1, -1}});
//    for (int i = 0; i < tets.size(); i++) {
//        for (int j = 0; j < 4; j++) {
//            if (opp_t_ids[i][j] >= 0)
//                continue;
//
//            std::vector<int> n_t_ids;
//            set_intersection_sorted(conn_tets[tets[i][(j + 1) % 4]], conn_tets[tets[i][(j + 2) % 4]],
//                                    conn_tets[tets[i][(j + 3) % 4]], n_t_ids);
//
//            assert(!n_t_ids.empty());
//            if (n_t_ids.size() < 2)
//                continue;
//
//            int n_t_id = n_t_ids[0] == i ? n_t_ids[1] : n_t_ids[0];
//            opp_t_ids[i][j] = n_t_id;
//            for (int k = 0; k < 4; k++) {
//                if (tets[n_t_id][k] != tets[i][(j + 1) % 4] && tets[n_t_id][k] != tets[i][(j + 2) % 4]
//                    && tets[n_t_id][k] != tets[i][(j + 3) % 4]) {
//                    opp_t_ids[n_t_id][k] = i;
//                    break;
//                }
//            }
//        }
//    }
//    int tmp_cnt = 0;
//    for (int i = 0; i < cut_t_ids.size(); i++) {
//        for (int j = 0; j < 4; j++) {
//            if (opp_t_ids[i][j] < 0
//                && is_v_on_plane(tets[i][(j + 1) % 4])
//                && is_v_on_plane(tets[i][(j + 2) % 4])
//                && is_v_on_plane(tets[i][(j + 3) % 4])) {
//                tmp_cnt++;
//                cout << i << " gt_id " << cut_t_ids[i] << " should include the neighbor!!!" << endl;
//                cout << "opp_t_id = " << get_opp_t_id(cut_t_ids[i], j, mesh) << endl;
//                cout << is_snapped[tets[i][(j + 1) % 4]] << " "
//                     << is_snapped[tets[i][(j + 2) % 4]] << " "
//                     << is_snapped[tets[i][(j + 3) % 4]] << endl;
//                pausee();
//            }
//        }
//    }
//    if(tmp_cnt>0) {
//        cout << tmp_cnt << " tets missing!!" << endl;
//    }
//    //fortest
}

int floatTetWild::CutMesh::project_to_plane(int input_vertices_size) {
    is_projected.resize(v_ids.size(), false);

    int cnt = 0;
    for (int i = 0; i < is_snapped.size(); i++) {
        if (!is_snapped[i] || is_projected[i])
            continue;
        if (v_ids[i] < input_vertices_size)
            continue;
        Scalar dist = get_to_plane_dist(mesh.tet_vertices[v_ids[i]].pos);
//        if (std::abs(dist) <= 1e-9) {
//            cnt++;
//            continue;
//        }
        Vector3 proj_p = mesh.tet_vertices[v_ids[i]].pos - p_n * dist;
//        cout << get_to_plane_dist(proj_p) << endl;
        bool is_snappable = true;
        for (int t_id: mesh.tet_vertices[v_ids[i]].conn_tets) {
            int j = mesh.tets[t_id].find(v_ids[i]);
            if (is_inverted(mesh, t_id, j, proj_p)) {
                is_snappable = false;
                break;
            }
        }
        if (is_snappable) {
            mesh.tet_vertices[v_ids[i]].pos = proj_p;
            is_projected[i] = true;
            to_plane_dists[i] = get_to_plane_dist(proj_p);
            cnt++;
        }
    }
//    cout<<std::count(is_snapped.begin(), is_snapped.end(), true)<<"/"<<cnt<<endl;
    return cnt;
}

void floatTetWild::CutMesh::revert_totally_snapped_tets(int a, int b) {
//    return;

    int cnt = 0;
    for (int i = a; i < b; i++) {
        const auto &t = tets[i];
        if (is_v_on_plane(t[0]) && is_v_on_plane(t[1]) && is_v_on_plane(t[2]) && is_v_on_plane(t[3])) {
            auto tmp_t = t;
            std::sort(tmp_t.begin(), tmp_t.end(), [&](int a, int b) {
                return fabs(to_plane_dists[a]) > fabs(to_plane_dists[b]);
            });
            for (int j = 0; j < 3; j++) {
                if (is_snapped[tmp_t[j]] == true) {
                    is_snapped[tmp_t[j]] = false;
//                    //fortest
//                    cnt++;
//                    if(j!=0){
//                        cout<<"snapping j!=0"<<endl;
//                        pausee();
//                    }
//                    //fortest
                    break;
                }
            }
        }
    }
}

bool floatTetWild::CutMesh::get_intersecting_edges_and_points(std::vector<Vector3> &points,
                                                              std::map<std::array<int, 2>, int>& map_edge_to_intersecting_point,
                                                              std::vector<int>& subdivide_t_ids) {
    std::vector<std::array<int, 2>> edges;
    for (auto &t: tets) {
        for (int j = 0; j < 3; j++) {
            std::array<int, 2> e;
            if (t[0] < t[j + 1])
                e = {{t[0], t[j + 1]}};
            else
                e = {{t[j + 1], t[0]}};
            edges.push_back(e);
            e = {{t[j + 1], t[mod3(j + 1) + 1]}};
            if (e[0] > e[1])
                std::swap(e[0], e[1]);
            edges.push_back(e);
        }
    }
    vector_unique(edges);

    std::vector<int> e_v_ids;
    for (int i = 0; i < edges.size(); i++) {
        const auto &e = edges[i];
        if (is_v_on_plane(e[0]) || is_v_on_plane(e[1]))
            continue;
        if (to_plane_dists[e[0]] * to_plane_dists[e[1]] >= 0)
            continue;

//        if (to_plane_dists[e[0]] > 0 && to_plane_dists[e[1]] > 0
//            || to_plane_dists[e[0]] < 0 && to_plane_dists[e[1]] < 0)
//            continue;

//        timer.start();
        int v1_id = v_ids[e[0]];
        int v2_id = v_ids[e[1]];
        Vector3 p;
        Scalar _;
        bool is_result = seg_plane_intersection(mesh.tet_vertices[v1_id].pos, mesh.tet_vertices[v2_id].pos,
                                                p_vs[0], p_n, p, _);
//        time_get_intersecting_edges_and_points2 += timer.getElapsedTime();
        if (!is_result) {
            //fortest
            cout << "seg_plane_intersection no result!" << endl;
            cout << to_plane_dists[e[0]] << ", " << to_plane_dists[e[1]] << endl;
            cout << get_to_plane_dist(mesh.tet_vertices[v1_id].pos) << ", "
                 << get_to_plane_dist(mesh.tet_vertices[v2_id].pos) << endl;
            cout << "e[0] = " << e[0] << endl;
            //fortest
            return false;
        }

//        Scalar dist1 = (p-mesh.tet_vertices[v1_id].pos).squaredNorm();
//        Scalar dist2 = (p-mesh.tet_vertices[v2_id].pos).squaredNorm();
//        if(dist1 <= SCALAR_ZERO_2){
//            cout<<"snapped e[0]"<<endl;
//            is_snapped[e[0]] = true;
//            continue;
//        }
//        if(dist2 <= SCALAR_ZERO_2){
//            cout<<"snapped e[1]"<<endl;
//            is_snapped[e[1]] = true;
//            continue;
//        }

//        timer.start();
        points.push_back(p);
        if (v1_id < v2_id)
            map_edge_to_intersecting_point[{{v1_id, v2_id}}] = points.size() - 1;
        else
            map_edge_to_intersecting_point[{{v2_id, v1_id}}] = points.size() - 1;
//        time_get_intersecting_edges_and_points3 += timer.getElapsedTime();

//        timer.start();
//        std::vector<int> tmp;
//        set_intersection(mesh.tet_vertices[v1_id].conn_tets, mesh.tet_vertices[v2_id].conn_tets, tmp);
//        subdivide_t_ids.insert(subdivide_t_ids.end(), tmp.begin(), tmp.end());
        e_v_ids.push_back(v1_id);
        e_v_ids.push_back(v2_id);
//        time_get_intersecting_edges_and_points4 += timer.getElapsedTime();
    }
//    timer.start();
    vector_unique(e_v_ids);
    for (int v_id: e_v_ids)
        subdivide_t_ids.insert(subdivide_t_ids.end(), mesh.tet_vertices[v_id].conn_tets.begin(),
                               mesh.tet_vertices[v_id].conn_tets.end());
    vector_unique(subdivide_t_ids);
//    time_get_intersecting_edges_and_points1 += timer.getElapsedTime();

    return true;
}

void floatTetWild::CutMesh::get_one_ring_t_ids(std::vector<int> &old_t_ids, std::vector<int> &neighbor_t_ids) {
//    std::vector<int> tmp_lv_ids;
//    for (int i = 0; i < tets.size(); i++) {
//        for (int j = 0; j < 4; j++) {
//            if (opp_t_ids[i][j] < 0) {
//                tmp_lv_ids.push_back(tets[i][(j + 1) % 4]);
//                tmp_lv_ids.push_back(tets[i][(j + 2) % 4]);
//                tmp_lv_ids.push_back(tets[i][(j + 3) % 4]);
//            }
//        }
//    }
//    vector_unique(tmp_lv_ids);
//
//    std::vector<int> tmp_t_ids;
//    for (int lv_id:tmp_lv_ids) {
//        tmp_t_ids.insert(neighbor_t_ids.end(), mesh.tet_vertices[v_ids[lv_id]].conn_tets.begin(),
//                         mesh.tet_vertices[v_ids[lv_id]].conn_tets.end());
//    }
//    vector_unique(tmp_t_ids);
//
//    std::set_difference(tmp_t_ids.begin(), tmp_t_ids.end(), old_t_ids.begin(), old_t_ids.end(),
//                        std::back_inserter(neighbor_t_ids));
}

bool floatTetWild::CutMesh::check() {
    return true;

    bool is_good = true;
    for (auto &m:map_v_ids) {
        int gv_id = m.first;
        int lv_id = m.second;

        Scalar dist = get_to_plane_dist(mesh.tet_vertices[gv_id].pos);
        if (std::fabs(dist) < mesh.params.eps_coplanar && to_plane_dists[lv_id] != 0 && !is_snapped[lv_id]) {
            cout << "wrong vertex in cut mesh" << endl;
            cout << dist << endl;
            cout << to_plane_dists[lv_id] << endl;
            cout << is_snapped[lv_id] << endl;
            is_good = false;
        }
    }
    return is_good;
}