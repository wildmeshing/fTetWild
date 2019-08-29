//
// Created by Yixin Hu on 2019-08-27.
//

#include <floattetwild/TriangleInsertion.h>

#include <floattetwild/FloatTetCutting.h>
#include <floattetwild/FloatTetCuttingCheck.h>
#include <floattetwild/FloatTetCuttingParallel.h>
#include <floattetwild/LocalOperations.h>
#include <floattetwild/Predicates.hpp>
#include <floattetwild/MeshIO.hpp>
#include <floattetwild/auto_table.hpp>
#include <floattetwild/Logger.hpp>
#include <floattetwild/intersections.h>


#include <igl/writeSTL.h>
#include <igl/Timer.h>

#ifdef FLOAT_TETWILD_USE_TBB
#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include <tbb/atomic.h>
#include <floattetwild/FloatTetCuttingParallel.h>
#include <tbb/concurrent_queue.h>
#endif

#include <bitset>
#include <numeric>
#include <unordered_map>

void floatTetWild::insert_triangles(const std::vector<Vector3> &input_vertices,
        const std::vector<Vector3i> &input_faces, const std::vector<int> &input_tags,
        Mesh &mesh, std::vector<bool> &is_face_inserted, AABBWrapper &tree, bool is_again) {
    //todo: mark inserted faces on mesh before calling this function!!!

    std::vector<Vector3> new_vertices;
    std::vector<std::array<int, 4>> new_tets;
    for (int i = 0; i < input_faces.size(); i++) {
        if (is_face_inserted[i])
            continue;
        if(insert_one_triangle(i, input_vertices, input_faces, input_tags, mesh, tree, is_again))
            is_face_inserted[i] = true;
    }

    //todo: preserve boundary

}

bool floatTetWild::insert_one_triangle(int insert_f_id, const std::vector<Vector3> &input_vertices,
        const std::vector<Vector3i> &input_faces, const std::vector<int> &input_tags,
        Mesh &mesh, AABBWrapper &tree, bool is_again) {

    std::array<Vector3, 3> vs = {{input_vertices[input_faces[insert_f_id][0]],
                                         input_vertices[input_faces[insert_f_id][1]],
                                         input_vertices[input_faces[insert_f_id][2]]}};
    Vector3 n = (vs[1] - vs[0]).cross(vs[2] - vs[0]);
    int t = get_t(vs[0], vs[1], vs[2]);

    /////
    std::vector<int> cut_t_ids;
    //todo: find tet for cutting, BFS!!


    /////
    CutMesh cut_mesh(mesh, n, vs);
    cut_mesh.construct(cut_t_ids);
    if (cut_mesh.snap_to_plane())
        cut_mesh.expand(cut_t_ids);

    /////
    std::vector<Vector3> points;
    std::map<std::array<int, 2>, int> map_edge_to_intersecting_point;
    std::vector<int> subdivide_t_ids;
    if (!cut_mesh.get_intersecting_edges_and_points(points, map_edge_to_intersecting_point, subdivide_t_ids))
        return false;

    /////
    std::vector<MeshTet> new_tets;
    std::vector<int> modified_t_ids;
    if (!subdivide_tets(mesh, points, map_edge_to_intersecting_point, subdivide_t_ids, new_tets, modified_t_ids))
        return false;

    if (!is_again) {
        ///vs
        const int old_v_size = mesh.tet_vertices.size();
        mesh.tet_vertices.resize(mesh.tet_vertices.size() + points.size());
        for (int i = 0; i < points.size(); i++) {
            mesh.tet_vertices[old_v_size + i].pos = points[i];
            //todo: tags???
        }

        ///tets
        mesh.tets.reserve(mesh.tets.size() + new_tets.size());
        for (int i = 0; i < new_tets.size(); i++) {
            if (i < modified_t_ids.size())
                mesh.tets[modified_t_ids[i]] = new_tets[i];
            else
                mesh.tets.push_back(new_tets[i]);

            //todo: tags???
        }
    } else {
        //todo
    }

    return true;
}

bool floatTetWild::subdivide_tets(Mesh& mesh, std::vector<Vector3>& points,
                    std::map<std::array<int, 2>, int>& map_edge_to_intersecting_point,
                    const std::vector<int>& subdivide_t_ids,
                    std::vector<MeshTet>& new_tets, std::vector<int>& modified_t_ids) {
    static const std::array<std::array<int, 2>, 6> t_es = {{{{0, 1}}, {{1, 2}}, {{2, 0}}, {{0, 3}}, {{1, 3}}, {{2, 3}}}};
    static const std::array<std::array<int, 3>, 4> t_f_es = {{{{1, 5, 4}}, {{5, 3, 2}}, {{3, 0, 4}}, {{0, 1, 2}}}};
    static const std::array<std::array<int, 3>, 4> t_f_vs = {{{{3, 1, 2}}, {{0, 2, 3}}, {{1, 3, 0}}, {{2, 0, 1}}}};

    for (int t_id: subdivide_t_ids) {
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
        if (config_id == 0) //no intersection
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
        auto check_config = [&](int diag_config_id, std::vector<std::pair<int, Vector3>> &centroids) {
            const std::vector<Vector4i> &config = CutTable::get_tet_conf(config_id, diag_config_id);
            Scalar min_q;
            int cnt = 0;
            for (const auto &tet: config) {
                std::array<Vector3, 4> vs;
                for (int j = 0; j < 4; j++) {
                    if (map_lv_to_v_id.find(tet[j]) == map_lv_to_v_id.end()) {
//                        Vector3 c;//todo: compute centroid -> vs[j]
                        centroids.push_back(std::make_pair(tet[j], vs[j]));
                    } else {
                        int v_id = map_lv_to_v_id[tet[j]];
                        if (v_id < v_size)
                            vs[j] = mesh.tet_vertices[v_id].pos;
                        else
                            vs[j] = points[v_id - v_size];
                    }
                }

                Scalar volumn = Predicates::orient_3d_volume(vs[0], vs[1], vs[2], vs[3]);
                if (cnt == 0)
                    min_q = volumn;
                else if (volumn < min_q)
                    min_q = volumn;
            }

            return min_q;
        };

        int diag_config_id = 0;
        std::vector<std::pair<int, Vector3>> centroids;
        if (!my_diags.empty()) {
            const auto &all_diags = CutTable::get_diag_confs(config_id);
            std::vector<std::pair<int, Scalar>> min_qualities;
            std::vector<std::vector<std::pair<int, Vector3>>> all_centroids;

            for (int i = 0; i < all_diags.size(); i++) {
                if (my_diags != all_diags[i])
                    continue;

                std::map<int, int> tmp_map;
                std::vector<std::pair<int, Vector3>> tmp_centroids;
                Scalar min_q = check_config(i, tmp_centroids);
                if (min_q < SCALAR_ZERO_3)
                    continue;
                min_qualities.push_back(std::make_pair(i, min_q));
                all_centroids.push_back(tmp_centroids);
            }
            std::sort(min_qualities.begin(), min_qualities.end(),
                      [](const std::pair<int, Scalar> &a, const std::pair<int, Scalar> &b) {
                          return a.second < b.second;
                      });

            if (min_qualities.back().second < SCALAR_ZERO_3) // if tet quality is too bad
                return false;

            int diag_config_id = min_qualities.back().first;
            centroids = all_centroids[diag_config_id];
        }

        std::sort(centroids.begin(), centroids.end(),
                  [](const std::pair<int, Vector3> &a, const std::pair<int, Vector3> &b) {
                      return a.first < b.first;
                  });
        for (int i = 0; i < centroids.size(); i++) {
            map_lv_to_v_id[centroids[i].first] = vp_size + i;
            points.push_back(centroids[i].second);
        }

        //add new tets
        const auto &config = CutTable::get_tet_conf(config_id, diag_config_id);
        const auto &new_is_surface_fs = CutTable::get_surface_conf(config_id, diag_config_id);
        const auto &new_local_f_ids = CutTable::get_face_id_conf(config_id, diag_config_id);
        for (const auto &t:config) {
            new_tets.push_back(MeshTet(map_lv_to_v_id[t[0]], map_lv_to_v_id[t[1]],
                                       map_lv_to_v_id[t[2]], map_lv_to_v_id[t[3]]));
            //todo: mark is_surface_fs and is_bbox_fs

        }
    }

    return true;
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
                           map_v_ids[mesh.tets[cut_t_ids[i]][2]]}};
    }

    std::vector<std::vector<int>> conn_tets(v_ids.size());
    for (int i = 0; i < tets.size(); i++) {
        for (int j = 0; j < 4; j++)
            conn_tets[tets[i][j]].push_back(i);
    }

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
    for (auto &v:map_v_ids) {
        int v_id = v.first;
        int lv_id = v.second;

        int ori = Predicates::orient_3d(p_vs[0], p_vs[1], p_vs[2], mesh.tet_vertices[v_id].pos);
        if (ori == Predicates::ORI_ZERO) {
            to_plane_dists[lv_id] = 0;
            continue;
        }
        to_plane_dists[lv_id] = get_to_plane_dist(mesh.tet_vertices[lv_id].pos);
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
    while (true) {
        std::vector<std::array<int, 5>> new_opp_t_ids;
        for (int i = 0; i < opp_t_ids.size(); i++) {
            for (int j = 0; j < 4; j++) {
                if (opp_t_ids[i][j] >= 0)
                    continue;
                if (!is_snapped[tets[i][(j + 1) % 4]] && !is_snapped[tets[i][(j + 2) % 4]]
                    && !is_snapped[tets[i][(j + 3) % 4]])
                    continue;
                std::vector<int> n_t_ids;
                set_intersection(mesh.tet_vertices[v_ids[tets[i][(j + 1) % 4]]].conn_tets,
                                 mesh.tet_vertices[v_ids[tets[i][(j + 2) % 4]]].conn_tets,
                                 mesh.tet_vertices[v_ids[tets[i][(j + 3) % 4]]].conn_tets,
                                 n_t_ids);
                if (n_t_ids.size() == 1)
                    continue;
                int new_t_id = n_t_ids[0] == cut_t_ids[i] ? n_t_ids[1] : n_t_ids[0];
                new_opp_t_ids.push_back({{-1, -1, -1, -1, new_t_id}});
                for (int j = 0; j < 4; j++) {
                    if (mesh.tets[new_t_id][j] != v_ids[tets[i][(j + 1) % 4]]
                        && mesh.tets[new_t_id][j] != v_ids[tets[i][(j + 2) % 4]]
                        && mesh.tets[new_t_id][j] != v_ids[tets[i][(j + 3) % 4]]) {
                        new_opp_t_ids.back()[j] = i;
                        break;
                    }
                }
            }
        }
        std::sort(new_opp_t_ids.begin(), new_opp_t_ids.end(),
                  [&](const std::array<int, 5> &a, const std::array<int, 5> &b) {
                      return a.back() < b.back();
                  });
        for (int i = 0; i < new_opp_t_ids.size() - 1; i++) {
            if (new_opp_t_ids[i].back() == new_opp_t_ids[i + 1].back()) {
                for (int j = 0; j < 4; j++) {
                    if (new_opp_t_ids[i][j] >= 0) {
                        assert(new_opp_t_ids[i + 1][j] < 0);
                        new_opp_t_ids[i + 1][j] = new_opp_t_ids[i][j];
                        break;
                    }
                }
                new_opp_t_ids.erase(new_opp_t_ids.begin() + i);
                i--;
            }
        }

        //todo: update CutMesh, topo & snapping
        bool snapped = false;
        const int old_tets_size = tets.size();
        for (int i = 0; i < new_opp_t_ids.size(); i++) {
            ///
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

            ///
            cut_t_ids.push_back(new_opp_t_ids[i].back());

            int t_id = tets.size();
            tets.emplace_back();
            auto &t = tets.back();
            for (int j = 0; j < 4; j++) {
                int v_id = mesh.tets[new_opp_t_ids[i].back()][j];
                if (map_v_ids.find(v_id) == map_v_ids.end()) {
                    v_ids.push_back(v_id);
                    int lv_id = v_ids.size() - 1;
                    map_v_ids[v_id] = lv_id;
                    to_plane_dists[lv_id] = get_to_plane_dist(mesh.tet_vertices[lv_id].pos);
                    if (std::abs(to_plane_dists[lv_id]) < mesh.params.eps_2_coplanar) {
                        is_snapped[lv_id] = true;
                        snapped = true;
                    }
                }
                t[j] = map_v_ids[v_id];
            }

            opp_t_ids.emplace_back();
            auto &opp = opp_t_ids.back();
            for (int j = 0; j < 4; j++) {
                opp[j] = new_opp_t_ids[i][j];
                if (new_opp_t_ids[i][j] < 0)
                    continue;
                int opp_t_id = new_opp_t_ids[i][j];
                for (int k = 0; k < 4; k++) {
                    if (tets[opp_t_id][k] != t[(j + 1) % 4] && tets[opp_t_id][k] != t[(j + 2) % 4]
                        && tets[opp_t_id][k] != t[(j + 3) % 4]) {
                        opp_t_ids[opp_t_id][k] = t_id;
                        break;
                    }
                }
            }
        }
        if (old_tets_size == tets.size())
            break;
        if (!snapped)
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

    std::vector<bool> is_intersected(edges.size(), false);
    for (int i = 0; i < edges.size(); i++) {
        auto &e = edges[i];
        if (is_v_on_plane(e[0]) || is_v_on_plane(e[1]))
            continue;
        if (to_plane_dists[e[0]] > 0 && to_plane_dists[e[1]] > 0
            || to_plane_dists[e[0]] < 0 && to_plane_dists[e[1]] < 0)
            continue;

        Vector3 p;
        Scalar _;
        bool is_result = seg_plane_intersection(mesh.tet_vertices[v_ids[e[0]]].pos, mesh.tet_vertices[v_ids[e[1]]].pos,
                                                p_vs[0], p_n, p, _);
        if (!is_result)
            return false;

        points.push_back(p);
        map_edge_to_intersecting_point[e] = points.size()-1;
    }


    for (int i = 0; i < edges.size(); i++) {
        if (!is_intersected[i])
            continue;
        auto &e = edges[i];
        std::vector<int> tmp;
        set_intersection(mesh.tet_vertices[v_ids[e[0]]].conn_tets, mesh.tet_vertices[v_ids[e[1]]].conn_tets, tmp);
        subdivide_t_ids.insert(subdivide_t_ids.end(), tmp.begin(), tmp.end());
    }
    vector_unique(subdivide_t_ids);

    return true;
}