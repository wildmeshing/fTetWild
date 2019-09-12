//
// Created by Yixin Hu on 9/12/19.
//

#include <floattetwild/CutMesh.h>
#include <floattetwild/TriangleInsertion.h>
#include <floattetwild/LocalOperations.h>
#include <floattetwild/Predicates.hpp>
#include <floattetwild/intersections.h>

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

    std::vector<std::vector<int>> conn_tets(v_ids.size());
    for (int i = 0; i < tets.size(); i++) {
        for (int j = 0; j < 4; j++)
            conn_tets[tets[i][j]].push_back(i);
    }

    opp_t_ids.resize(tets.size(), {{-1, -1, -1, -1}});
    for (int i = 0; i < tets.size(); i++) {
        for (int j = 0; j < 4; j++) {
            if (opp_t_ids[i][j] >= 0)
                continue;
            std::vector<int> n_t_ids;
            set_intersection(conn_tets[tets[i][(j + 1) % 4]], conn_tets[tets[i][(j + 2) % 4]],
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
        if (ori == Predicates::ORI_POSITIVE && to_plane_dists[lv_id] < 0
            || ori == Predicates::ORI_NEGATIVE && to_plane_dists[lv_id] > 0)
            to_plane_dists[lv_id] = -to_plane_dists[lv_id];

        if (std::abs(to_plane_dists[lv_id]) < mesh.params.eps_2_coplanar) {
            is_snapped[lv_id] = true;
            snapped = true;
        }
    }
    for (auto &t: tets) {
        if (is_v_on_plane(t[0]) && is_v_on_plane(t[1]) && is_v_on_plane(t[2]) && is_v_on_plane(t[3])) {
            auto tmp_t = t;
            std::sort(tmp_t.begin(), tmp_t.end(), [&](int a, int b) {
                return to_plane_dists[a] < to_plane_dists[b];
            });
            for (int j = 3; j >= 0; j--) {
                if (is_snapped[tmp_t.back()] == true)
                    is_snapped[tmp_t.back()] = false;
            }
        }
    }

    return snapped;
}

void floatTetWild::CutMesh::expand(std::vector<int>& cut_t_ids) {
//    cout<<"expanding"<<endl;

    std::vector<std::vector<int>> conn_tets(v_ids.size());
    for (int i = 0; i < tets.size(); i++) {
        for (int j = 0; j < 4; j++)
            conn_tets[tets[i][j]].push_back(i);
    }

    while (true) {
//        //fortest
//        cout<<"loop"<<endl;
//        cout<<"cut_t_ids.size() = "<<cut_t_ids.size()<<endl;

        std::vector<std::array<int, 5>> new_opp_t_ids;
        for (int i = 0; i < opp_t_ids.size(); i++) {
            for (int j = 0; j < 4; j++) {
                if (opp_t_ids[i][j] >= 0)
                    continue;
                if (!is_snapped[tets[i][(j + 1) % 4]] && !is_snapped[tets[i][(j + 2) % 4]]
                    && !is_snapped[tets[i][(j + 3) % 4]])
                    continue;

//                std::vector<int> n_t_ids;
//                set_intersection(mesh.tet_vertices[v_ids[tets[i][(j + 1) % 4]]].conn_tets,
//                                 mesh.tet_vertices[v_ids[tets[i][(j + 2) % 4]]].conn_tets,
//                                 mesh.tet_vertices[v_ids[tets[i][(j + 3) % 4]]].conn_tets,
//                                 n_t_ids);
//                if (n_t_ids.size() == 1)
//                    continue;
//                int new_t_id = n_t_ids[0] == cut_t_ids[i] ? n_t_ids[1] : n_t_ids[0];
//
//                //fortest
//                if(get_opp_t_id(cut_t_ids[i], j, mesh)!=new_t_id){
//                    cout<<"get_opp_t_id(cut_t_ids[i], j, mesh)!=new_t_id"<<endl;
//                    pausee();
//                }
//                //fortest

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

//                //fortest
//                if(n_gt_id==5539) {
//                    cout << "n_gt_id==5539" << endl;
//                    std::vector<int> n_t_ids;
//                    set_intersection(mesh.tet_vertices[v_ids[tets[i][(j + 1) % 4]]].conn_tets,
//                                     mesh.tet_vertices[v_ids[tets[i][(j + 2) % 4]]].conn_tets,
//                                     mesh.tet_vertices[v_ids[tets[i][(j + 3) % 4]]].conn_tets,
//                                     n_t_ids);
//                    vector_print(n_t_ids);
//                    cout << "n_gt_id " << n_gt_id << ": ";
//                    mesh.tets[n_gt_id].print();
//                    cout << (std::find(cut_t_ids.begin(), cut_t_ids.end(), n_gt_id) == cut_t_ids.end()) << endl;
//                    cout << "cut_t_ids[i] " << cut_t_ids[i] << ": ";
//                    mesh.tets[cut_t_ids[i]].print();
//                    for (int k = 0; k < 5; k++)
//                        cout << new_opp_t_ids.back()[k] << " ";
//                    cout << endl;
//                    pausee();
//                }
//                //fortest
            }
        }
        if(new_opp_t_ids.empty())
            return;

        std::sort(new_opp_t_ids.begin(), new_opp_t_ids.end(),
                  [&](const std::array<int, 5> &a, const std::array<int, 5> &b) {
                      return a.back() < b.back();
                  });
//        //fortest
//        for(auto& m: new_opp_t_ids) {
//            for (int k = 0; k < 5; k++)
//                cout << m[k] << " ";
//            cout<<endl;
//        }
//        pausee();
//        //fortest
        for (int i = 0; i < new_opp_t_ids.size() - 1; i++) {
            if (new_opp_t_ids[i].back() == new_opp_t_ids[i + 1].back()) {
                for (int j = 0; j < 4; j++) {
                    if (new_opp_t_ids[i][j] >= 0) {
//                        //fortest
//                        if(new_opp_t_ids[i + 1][j] >= 0){
//                            cout<<"new_opp_t_ids[i + 1][j] >= 0"<<endl;
//                            cout<<cut_t_ids[new_opp_t_ids[i + 1][j]]<<endl;
//                            mesh.tets[cut_t_ids[new_opp_t_ids[i + 1][j]]].print();
//                            cout<<cut_t_ids[new_opp_t_ids[i][j]]<<endl;
//                            mesh.tets[cut_t_ids[new_opp_t_ids[i][j]]].print();
//                            pausee();
//                        }
//                        //fortest
                        new_opp_t_ids[i + 1][j] = new_opp_t_ids[i][j];
                        break;
                    }
                }
                new_opp_t_ids.erase(new_opp_t_ids.begin() + i);
                i--;
            }
        }
//        //fortest
//        cout<<"new_opp_t_ids.size() = "<<new_opp_t_ids.size()<<endl;
//        for(auto& m: new_opp_t_ids) {
//            for (int k = 0; k < 5; k++)
//                cout << m[k] << " ";
//            cout<<endl;
//        }
//        pausee();
//        //fortest

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
                if (cnt_neg == 0 || cnt_pos == 0) {
//                    cout << new_opp_t_ids[i].back() << " is skipped" << endl;//fortest
                    continue;
                }
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
                    to_plane_dists.push_back(get_to_plane_dist(mesh.tet_vertices[lv_id].pos));
                    if (std::abs(to_plane_dists[lv_id]) < mesh.params.eps_2_coplanar) {
                        is_snapped.push_back(true);
//                        snapped = true;
                    } else
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
//                        //fortest
//                        if (opp_t_ids[opp_t_id][k] >= 0) {
//                            cout << "old opp_t_ids[opp_t_id][k]>=0!!" << endl;
//                            cout << "j = " << j << endl;
//                            cout << "k = " << k << endl;
//                            cout << opp_t_ids[opp_t_id][0] << " " << opp_t_ids[opp_t_id][1] << " "
//                                 << opp_t_ids[opp_t_id][2] << " " << opp_t_ids[opp_t_id][3] << endl;
//                            cout << opp_t_ids.back()[0] << " " << opp_t_ids.back()[1] << " "
//                                 << opp_t_ids.back()[2] << " " << opp_t_ids.back()[3] << endl;
//                            cout << "opp_t_id " << opp_t_id << ": "
//                                 << tets[opp_t_id][0] << " " << tets[opp_t_id][1] << " "
//                                 << tets[opp_t_id][2] << " " << tets[opp_t_id][3] << endl;
//                            cout << "opp_t_ids[opp_t_id][k] " << opp_t_ids[opp_t_id][k] << ": "
//                                 << tets[opp_t_ids[opp_t_id][k]][0] << " "
//                                 << tets[opp_t_ids[opp_t_id][k]][1] << " "
//                                 << tets[opp_t_ids[opp_t_id][k]][2] << " "
//                                 << tets[opp_t_ids[opp_t_id][k]][3] << endl;
//                            cout << "t_id " << t_id << ": "
//                                 << tets[t_id][0] << " " << tets[t_id][1] << " " << tets[t_id][2] << " "
//                                 << tets[t_id][3] << endl;
//                            pausee();
//                        }
//                        //fortest
                        opp_t_ids[opp_t_id][k] = t_id;
                        break;
                    }
                }
            }
        }
        if (old_tets_size == tets.size())
            break;

        for (int i = old_tets_size; i < tets.size(); i++) {
            const auto &t = tets[i];
            if (is_v_on_plane(t[0]) && is_v_on_plane(t[1]) && is_v_on_plane(t[2]) && is_v_on_plane(t[3])) {
                auto tmp_t = t;
                std::sort(tmp_t.begin(), tmp_t.end(), [&](int a, int b) {
                    return to_plane_dists[a] < to_plane_dists[b];
                });
                for (int j = 3; j >= 0; j--) {
                    if (is_snapped[tmp_t.back()] == true)
                        is_snapped[tmp_t.back()] = false;
                }
            }
        }

//        if (!snapped)
//            break;
    }
}

bool floatTetWild::CutMesh::get_intersecting_edges_and_points(std::vector<Vector3> &points,
                                                              std::map<std::array<int, 2>, int>& map_edge_to_intersecting_point,
                                                              std::vector<int>& subdivide_t_ids) {
    std::vector<std::array<int, 2>> edges;
    for (auto &t: tets) {
        for (int j = 0; j < 3; j++) {
            std::array<int, 2> e;
            if (e[0] > e[1])
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

    for (int i = 0; i < edges.size(); i++) {
        auto &e = edges[i];
        if (is_v_on_plane(e[0]) || is_v_on_plane(e[1]))
            continue;
        if (to_plane_dists[e[0]] > 0 && to_plane_dists[e[1]] > 0
            || to_plane_dists[e[0]] < 0 && to_plane_dists[e[1]] < 0)
            continue;

        int v1_id = v_ids[e[0]];
        int v2_id = v_ids[e[1]];
        Vector3 p;
        Scalar _;
        bool is_result = seg_plane_intersection(mesh.tet_vertices[v1_id].pos, mesh.tet_vertices[v2_id].pos,
                                                p_vs[0], p_n, p, _);
        if (!is_result)
            return false;

        points.push_back(p);
        if (v1_id < v2_id)
            map_edge_to_intersecting_point[{{v1_id, v2_id}}] = points.size() - 1;
        else
            map_edge_to_intersecting_point[{{v2_id, v1_id}}] = points.size() - 1;

        std::vector<int> tmp;
        set_intersection(mesh.tet_vertices[v1_id].conn_tets, mesh.tet_vertices[v2_id].conn_tets, tmp);
        subdivide_t_ids.insert(subdivide_t_ids.end(), tmp.begin(), tmp.end());
    }
    vector_unique(subdivide_t_ids);

    return true;
}

void floatTetWild::CutMesh::get_one_ring_t_ids(std::vector<int> &old_t_ids, std::vector<int> &neighbor_t_ids) {
    std::vector<int> tmp_lv_ids;
    for (int i = 0; i < tets.size(); i++) {
        for (int j = 0; j < 4; j++) {
            if (opp_t_ids[i][j] < 0) {
                tmp_lv_ids.push_back(tets[i][(j + 1) % 4]);
                tmp_lv_ids.push_back(tets[i][(j + 2) % 4]);
                tmp_lv_ids.push_back(tets[i][(j + 3) % 4]);
            }
        }
    }
    vector_unique(tmp_lv_ids);

    std::vector<int> tmp_t_ids;
    for (int lv_id:tmp_lv_ids) {
        tmp_t_ids.insert(neighbor_t_ids.end(), mesh.tet_vertices[v_ids[lv_id]].conn_tets.begin(),
                         mesh.tet_vertices[v_ids[lv_id]].conn_tets.end());
    }
    vector_unique(tmp_t_ids);

    std::set_difference(tmp_t_ids.begin(), tmp_t_ids.end(), old_t_ids.begin(), old_t_ids.end(),
                        std::back_inserter(neighbor_t_ids));
}