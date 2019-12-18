// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include <floattetwild/FloatTetCuttingParallel.h>
#include <floattetwild/intersections.h>
#include <floattetwild/FloatTetCutting.h>
#include <floattetwild/LocalOperations.h>
#include <floattetwild/MeshImprovement.h>

#ifdef FLOAT_TETWILD_USE_TBB
#include <tbb/parallel_for.h>
#endif

#include <igl/Timer.h>

void floatTetWild::generate_coloring_graph(const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces,
                             const Mesh& mesh, std::vector<std::array<int, 2>>& graph, bool is_again) {
    igl::Timer igl_timer;

    igl_timer.start();
    auto &tets = mesh.tets;
    auto &tet_vertices = mesh.tet_vertices;

    std::vector<std::array<Vector3, 2>> graph_nodes(input_faces.size());
#ifdef FLOAT_TETWILD_USE_TBB
    tbb::parallel_for( size_t(0), size_t(input_faces.size()), [&]( size_t f_id ) {
#else
        for (int f_id = 0; f_id < input_faces.size(); f_id++) {
#endif
        Vector3 min_f, max_f;
        get_bbox_face(input_vertices[input_faces[f_id][0]], input_vertices[input_faces[f_id][0]],
                      input_vertices[input_faces[f_id][0]], min_f, max_f, mesh.params.eps);

        std::queue<int> t_ids_queue;
        std::vector<bool> is_visited(tets.size(), false);
        if (!is_again) {
            for (int t_id:tet_vertices[input_faces[f_id][0]].conn_tets) {
                t_ids_queue.push(t_id);
                is_visited[t_id] = true;
            }
        } else {
            for (int t_id = 0; t_id < tets.size(); t_id++) {
                if (tets[t_id].is_removed)
                    continue;
                is_visited[t_id] = true;
                Vector3 min_t, max_t;
                get_bbox_tet(tet_vertices[tets[t_id][0]].pos, tet_vertices[tets[t_id][1]].pos,
                             tet_vertices[tets[t_id][2]].pos, tet_vertices[tets[t_id][2]].pos,
                             min_t, max_t, mesh.params.eps);
                if (!is_bbox_intersected(min_f, max_f, min_t, max_t))
                    continue;
                t_ids_queue.push(t_id);
                break;
            }
        }

        bool is_init = true;
        Vector3 &min_b = graph_nodes[f_id][0];
        Vector3 &max_b = graph_nodes[f_id][1];
        while (!t_ids_queue.empty()) {
            int t_id = t_ids_queue.front();
            t_ids_queue.pop();

            Vector3 min_t, max_t;
            get_bbox_tet(tet_vertices[tets[t_id][0]].pos, tet_vertices[tets[t_id][1]].pos,
                         tet_vertices[tets[t_id][2]].pos, tet_vertices[tets[t_id][2]].pos,
                         min_t, max_t, mesh.params.eps);
            if (!is_bbox_intersected(min_f, max_f, min_t, max_t))
                continue;

            if (is_init) {
                min_b = min_t;
                max_b = max_t;
                is_init = false;
            } else {
                for (int j = 0; j < 3; j++) {
                    if (min_t[j] < min_b[j])
                        min_b[j] = min_t[j];
                    if (max_t[j] > max_b[j])
                        max_b[j] = max_t[j];
                }
            }

            for (int j = 0; j < 4; j++) {
                for (int n_t_id:tet_vertices[tets[t_id][j]].conn_tets) {
                    if (!is_visited[n_t_id]) {
                        t_ids_queue.push(n_t_id);
                        is_visited[n_t_id] = true;
                    }
                }
            }

        }

#ifdef FLOAT_TETWILD_USE_TBB
    });
#else
    }
#endif
    cout<<igl_timer.getElapsedTime()<<endl;

    igl_timer.start();
    for (int n = 0; n < graph_nodes.size(); n++) {
        for (int m = n + 1; m < graph_nodes.size(); m++) {
            if (is_bbox_intersected(graph_nodes[n][0], graph_nodes[n][1], graph_nodes[m][0], graph_nodes[m][1]))
                graph.push_back({{n, m}});
        }
    }
    cout<<igl_timer.getElapsedTime()<<endl;
}

void floatTetWild::generate_coloring(const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces,
                             const Mesh& mesh, std::vector<int> &colors, bool is_again) {
    std::vector<std::array<int, 2>> graph;
    generate_coloring_graph(input_vertices, input_faces, mesh, graph, is_again);
    std::vector<std::unordered_set<int>> conn_box_ids(input_faces.size());
    for (const auto &edge : graph) {
        conn_box_ids[edge[0]].insert(edge[1]);
        conn_box_ids[edge[1]].insert(edge[0]);
    }

    colors.resize(input_faces.size());
    std::fill(colors.begin(), colors.end(), -1);
    colors[0] = 0;

    std::vector<bool> available(input_faces.size(), true);
    for (int i = 1; i < input_faces.size(); ++i) {
        const auto &f = input_faces[i];
        const auto &ring = conn_box_ids[i];

        for (const auto n : ring) {
            if (colors[n] != -1)
                available[colors[n]] = false;
        }

        int first_available_col;
        for (first_available_col = 0; first_available_col < available.size(); first_available_col++) {
            if (available[first_available_col])
                break;
        }

        assert(available[first_available_col]);

        colors[i] = first_available_col;

        for (const auto n : ring) {
            if (colors[n] != -1)
                available[colors[n]] = true;
        }
    }
}

void floatTetWild::box_sets(const int threshold, const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces,
                             const Mesh& mesh, std::vector<std::vector<int>> &concurrent_sets, std::vector<int> &serial_set, bool is_again) {
    std::vector<int> colors;
    generate_coloring(input_vertices, input_faces, mesh, colors, is_again);
    int max_c = -1;
    for (const auto c : colors)
        max_c = std::max(max_c, int(c));

    concurrent_sets.clear();
    concurrent_sets.resize(max_c + 1);
    serial_set.clear();

    for (size_t i = 0; i < colors.size(); ++i) {
        const int col = colors[i];
        //removed vertex
        if (col < 0)
            serial_set.push_back(i);
        else
            concurrent_sets[col].push_back(i);
    }


    for (int i = concurrent_sets.size() - 1; i >= 0; --i) {
        if (concurrent_sets[i].size() < threshold) {
            serial_set.insert(serial_set.end(), concurrent_sets[i].begin(), concurrent_sets[i].end());
            concurrent_sets.erase(concurrent_sets.begin() + i);
        }
    }
}

#include <floattetwild/MeshIO.hpp>
void floatTetWild::partition_mesh(const std::vector<std::vector<int>>& partition_t_ids,
        const Mesh& meshin, std::vector<Mesh>& sub_meshes,
//        tbb::concurrent_vector<Mesh>& sub_meshes,
                                  const std::vector<std::array<std::vector<int>, 4>>& cut_f_ids,
                                  std::vector<std::vector<std::array<std::vector<int>, 4>>>& sub_cut_f_ids,
//                                  tbb::concurrent_vector<std::vector<std::array<std::vector<int>, 4>>>& sub_cut_f_ids,
                                  std::vector<std::array<int, 2>>& map_input_to_partition,
                                  std::vector<std::vector<std::array<int, 2>>>& v_p_ids, bool is_again) {
//    std::vector<std::vector<int>> partition_t_ids;
//    meshin.partition(N, partition_t_ids);

    std::vector<std::array<int, 2>> map_v_partition(meshin.tet_vertices.size(), {{-1, -1}});
//    for(int n = 0; n < partition_t_ids.size(); n++){
//        for(int t_id: partition_t_ids[n])
//            for(int j=0;j<4;j++)
//                map_v_partition[meshin.tets[t_id][j]] = n;
//    }

    sub_meshes.resize(partition_t_ids.size());
    sub_cut_f_ids.resize(partition_t_ids.size());
    v_p_ids.resize(partition_t_ids.size());
    for (int n = 0; n < partition_t_ids.size(); n++) {
        cout<<partition_t_ids[n].size()<<endl;
        ////partition topology
        sub_meshes[n].params = meshin.params;//copy params
        sub_meshes[n].tets.resize(partition_t_ids[n].size());
        sub_cut_f_ids[n].resize(partition_t_ids[n].size());

        sub_meshes[n].tet_vertices.reserve(sub_meshes[n].tets.size() / 6);
        std::vector<int> map_v_ids(meshin.tet_vertices.size(), -1);
        int cnt = 0;
        for (int i = 0; i < partition_t_ids[n].size(); i++) {
            int t_id = partition_t_ids[n][i];
            auto &t = meshin.tets[t_id];
            sub_cut_f_ids[n][i] = cut_f_ids[t_id];
            sub_meshes[n].tets[i] = t;//copy tags
//            sub_meshes[n].tets[i].opp_t_ids = {{-1,-1,-1,-1}};
            sub_meshes[n].tets[i].opp_t_ids = {{OPP_T_ID_UNKNOWN, OPP_T_ID_UNKNOWN, OPP_T_ID_UNKNOWN, OPP_T_ID_UNKNOWN}};
            for (int j = 0; j < 4; j++) {
                if (map_v_ids[t[j]] < 0) {
                    map_v_ids[t[j]] = cnt++;
                    sub_meshes[n].tet_vertices.push_back(meshin.tet_vertices[t[j]]);
                    sub_meshes[n].tet_vertices.back().conn_tets.clear();
                    bool is_visited = true;
                    int &p_id = map_v_partition[t[j]][0];
                    int &lv_id = map_v_partition[t[j]][1];
                    if (p_id < 0) {
                        p_id = n;
                        lv_id = map_v_ids[t[j]];
                        is_visited = false;
                    }
                    v_p_ids[n].push_back({{p_id, lv_id}});

                    if (!is_again) {
                        if (t[j] < map_input_to_partition.size()) {
                            if(is_visited)
                                map_input_to_partition[t[j]] = {{-p_id, lv_id}};
                            else
                                map_input_to_partition[t[j]] = {{p_id, lv_id}};
//                            if (map_input_to_partition[t[j]][0] >= 0 && map_input_to_partition[t[j]][0] != n)
////                                map_input_to_partition[t[j]] = {{-1, -1}};//input vertex is on the interface
//                                map_input_to_partition[t[j]] = {{-n, map_v_ids[t[j]]}};//input vertex is on the interface
//                            else
//                                map_input_to_partition[t[j]] = {{n, map_v_ids[t[j]]}};
                        }
                    }
                }
                sub_meshes[n].tets[i][j] = map_v_ids[t[j]];//map v_ids
            }
        }


        ////update connectivity
        //update conn_tets
        for(int i=0;i<sub_meshes[n].tets.size();i++) {
            for (int j = 0; j < 4; j++)
//                sub_meshes[n].tet_vertices[sub_meshes[n].tets[i][j]].conn_tets.insert(i);
                sub_meshes[n].tet_vertices[sub_meshes[n].tets[i][j]].conn_tets.push_back(i);
        }

        //update opp_t_ids
//        for(int t_id = 0;t_id<sub_meshes[n].tets.size();t_id++) {
//            auto &t = sub_meshes[n].tets[t_id];
//            for (int j = 0; j < 4; j++) {
//                if (t.opp_t_ids[j] >= 0)
//                    continue;
//                std::unordered_set<int> tmp;
//                set_intersection(sub_meshes[n].tet_vertices[t[(j + 1) % 4]].conn_tets,
//                                 sub_meshes[n].tet_vertices[t[(j + 2) % 4]].conn_tets, tmp);
//                std::vector<int> pair;
//                set_intersection(sub_meshes[n].tet_vertices[t[(j + 3) % 4]].conn_tets, tmp, pair);
//                if (pair.size() == 2) {
//                    int opp_t_id = pair[0] == t_id ? pair[1] : pair[0];
//                    t.opp_t_ids[j] = opp_t_id;
//                    auto &opp_t = sub_meshes[n].tets[opp_t_id];
//                    for (int k = 0; k < 4; k++) {
//                        if (opp_t[k] != t[(j + 1) % 4] && opp_t[k] != t[(j + 2) % 4] && opp_t[k] != t[(j + 3) % 4])
//                            opp_t.opp_t_ids[k] = t_id;
//                    }
//                }
//            }
//        }
//        MeshIO::write_mesh("hehe"+std::to_string(n)+".msh", sub_meshes[n]);
    }
}

void floatTetWild::merge_meshes(const std::vector<Mesh>& sub_meshes, //const tbb::concurrent_vector<Mesh>& sub_meshes,
        Mesh& mesh,
        const std::vector<std::vector<std::array<std::vector<int>, 4>>>& sub_cut_f_ids,
//        const tbb::concurrent_vector<std::vector<std::array<std::vector<int>, 4>>>& sub_cut_f_ids,
        std::vector<std::array<std::vector<int>, 4>>& cut_f_ids,
        const std::vector<std::array<int, 2>>& map_input_to_partition, std::vector<int>& map_input_to_merged_mesh,
        const std::vector<std::vector<std::array<int, 2>>>& v_p_ids) {
    ////merge topology
    mesh = sub_meshes.front();
    cut_f_ids = sub_cut_f_ids.front();

    if (sub_meshes.size() == 1)
        return;

    std::vector<int> incr(sub_meshes.size() + 1);
    incr[0] = 0;
    for (int n = 1; n < sub_meshes.size(); n++) {
        int old_v_num = mesh.tet_vertices.size();
        int old_t_num = mesh.tets.size();
        incr[n] = old_v_num;
        mesh.tet_vertices.insert(mesh.tet_vertices.end(),
                                 sub_meshes[n].tet_vertices.begin(), sub_meshes[n].tet_vertices.end());
        mesh.tets.insert(mesh.tets.end(), sub_meshes[n].tets.begin(), sub_meshes[n].tets.end());
        cut_f_ids.insert(cut_f_ids.end(), sub_cut_f_ids[n].begin(), sub_cut_f_ids[n].end());
        for (int i = old_t_num; i < mesh.tets.size(); i++) {
            for (int j = 0; j < 4; j++)
                mesh.tets[i][j] += old_v_num;
        }
    }
    incr.back() = mesh.tets.size();

    for (int iv_id = 0; iv_id < map_input_to_partition.size(); iv_id++) {
        int p_id = map_input_to_partition[iv_id][0];
        int lv_id = map_input_to_partition[iv_id][1];
        if (p_id < 0)
            map_input_to_merged_mesh[iv_id] = incr[-p_id] + lv_id;
        else
            map_input_to_merged_mesh[iv_id] = incr[p_id] + lv_id;
    }

    ////update connectivity
    for (auto &v: mesh.tet_vertices)
        v.conn_tets.clear();
    for (int i = 0; i < mesh.tets.size(); i++) {
        for (int j = 0; j < 4; j++)
//            mesh.tet_vertices[mesh.tets[i][j]].conn_tets.insert(i);
            mesh.tet_vertices[mesh.tets[i][j]].conn_tets.push_back(i);
    }

    for (int n = 0; n < v_p_ids.size(); n++) {
        for (int i = 0; i < v_p_ids[n].size(); i++) {//only happens on tet_vertices existing before cutting
            int m = v_p_ids[n][i][0];
            int lv_id = v_p_ids[n][i][1];
            if (m != n) {
                int old_v_id = incr[n] + i;
                int new_v_id = incr[m] + lv_id;
                mesh.tet_vertices[old_v_id].is_removed = true;
                for (int t_id: mesh.tet_vertices[old_v_id].conn_tets) {
                    int j = mesh.tets[t_id].find(old_v_id);
                    mesh.tets[t_id][j] = new_v_id;
//                    mesh.tet_vertices[new_v_id].conn_tets.insert(t_id);
                    mesh.tet_vertices[new_v_id].conn_tets.push_back(t_id);
                }
            }
        }
    }

//    MeshIO::write_mesh("hehehe.msh", mesh);
//    output_surface(mesh, "hehe_bbox");
//    //pausee();

    //update opp_t_ids
    for (auto &t: mesh.tets)
        t.opp_t_ids = {{OPP_T_ID_UNKNOWN, OPP_T_ID_UNKNOWN, OPP_T_ID_UNKNOWN, OPP_T_ID_UNKNOWN}};
//    for (int t_id = 0; t_id < mesh.tets.size(); t_id++) {
//        auto &t = mesh.tets[t_id];
//        for (int j = 0; j < 4; j++) {
//            if (t.opp_t_ids[j] >= 0)
//                continue;
//            std::vector<int> pair;
//            set_intersection(mesh.tet_vertices[t[(j + 1) % 4]].conn_tets,
//                             mesh.tet_vertices[t[(j + 2) % 4]].conn_tets,
//                             mesh.tet_vertices[t[(j + 3) % 4]].conn_tets, pair);
//            if (pair.size() == 2) {
//                int opp_t_id = pair[0] == t_id ? pair[1] : pair[0];
//                t.opp_t_ids[j] = opp_t_id;
//                auto &opp_t = mesh.tets[opp_t_id];
//                for (int k = 0; k < 4; k++) {
//                    if (opp_t[k] != t[(j + 1) % 4] && opp_t[k] != t[(j + 2) % 4] && opp_t[k] != t[(j + 3) % 4]) {
//                        opp_t.opp_t_ids[k] = t_id;
//                        break;
//                    }
//                }
//            }
//        }
//    }
}

bool floatTetWild::is_cutting_cross_partitions(Mesh& sub_mesh, const std::vector<int>& intersection_results_wn) {
#ifdef FLOAT_TETWILD_USE_TBB
    for (int t_id:intersection_results_wn) {
        for (int j = 0; j < 4; j++) {
            if (sub_mesh.tets[t_id].opp_t_ids[j] == OPP_T_ID_UNKNOWN) {
                sub_mesh.tets[t_id].opp_t_ids[j] = get_opp_t_id(sub_mesh, t_id, j);
                int opp_t_id = sub_mesh.tets[t_id].opp_t_ids[j];
                if (sub_mesh.tets[t_id].opp_t_ids[j] != OPP_T_ID_BOUNDARY) {
                    int k = sub_mesh.tets[opp_t_id].find_opp(sub_mesh.tets[t_id][mod4(j + 1)],
                                                             sub_mesh.tets[t_id][mod4(j + 2)],
                                                             sub_mesh.tets[t_id][mod4(j + 3)]);
                    sub_mesh.tets[opp_t_id].opp_t_ids[k] = t_id;
                }
            }

            if (sub_mesh.tets[t_id].opp_t_ids[j] == OPP_T_ID_BOUNDARY && sub_mesh.tets[t_id].is_bbox_fs[j] == NOT_BBOX)
                return true;
        }
    }
#endif
    return false;
}