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

#include <floattetwild/MeshImprovement.h>//fortest

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

#define III -1

void floatTetWild::insert_triangles(const std::vector<Vector3> &input_vertices,
        const std::vector<Vector3i> &input_faces, const std::vector<int> &input_tags,
        Mesh &mesh, std::vector<bool> &is_face_inserted, AABBWrapper &tree, bool is_again) {

    logger().info("triangle insertion start, #f = {}, #v = {}, #t = {}",
                  input_faces.size(), mesh.tet_vertices.size(), mesh.tets.size());

    /////
    std::vector < std::array < std::vector < int > , 4 >> track_surface_fs(mesh.tets.size());
    if (!is_again) {
        match_surface_fs(mesh, input_vertices, input_faces, is_face_inserted, track_surface_fs);
        //todo: ???
    } else {
        //todo: ??
    }
    logger().info("match_surface_fs done, matched #f = {}",
                  std::count(is_face_inserted.begin(), is_face_inserted.end(), true));

    /////
    std::vector <Vector3> new_vertices;
    std::vector <std::array<int, 4>> new_tets;
    for (int i = 0; i < input_faces.size(); i++) {
        if (is_face_inserted[i])
            continue;

        cout<<endl<<"fid "<<i<<endl;
        if (insert_one_triangle(i, input_vertices, input_faces, input_tags, mesh, track_surface_fs, tree, is_again))
            is_face_inserted[i] = true;

//        pausee();//fortest
        if(i == III)
            break;//fortest
    }
    logger().info("insert_one_triangle * n done, #v = {}, #t = {}", mesh.tet_vertices.size(), mesh.tets.size());
    logger().info("uninserted #f = {}", std::count(is_face_inserted.begin(), is_face_inserted.end(), false));

    //fortest
    check_track_surface_fs(mesh, track_surface_fs);

    /////
    //todo: preserve open boundary

    /////
    //todo: update mesh is_surface_fs
    mark_surface_fs(input_vertices, input_faces, track_surface_fs, mesh);
    logger().info("mark_surface_fs done");
}

bool floatTetWild::insert_one_triangle(int insert_f_id, const std::vector<Vector3> &input_vertices,
        const std::vector<Vector3i> &input_faces, const std::vector<int> &input_tags,
        Mesh &mesh, std::vector<std::array<std::vector<int>, 4>>& track_surface_fs,
        AABBWrapper &tree, bool is_again) {

    std::array<Vector3, 3> vs = {{input_vertices[input_faces[insert_f_id][0]],
                                         input_vertices[input_faces[insert_f_id][1]],
                                         input_vertices[input_faces[insert_f_id][2]]}};
    Vector3 n = (vs[1] - vs[0]).cross(vs[2] - vs[0]);
    n.normalize();
    int t = get_t(vs[0], vs[1], vs[2]);

    /////
    std::vector<int> cut_t_ids;
    find_cutting_tets(insert_f_id, input_faces, vs, mesh, cut_t_ids);

    //fortest
    myassert(!cut_t_ids.empty());
    //fortest

    /////
    CutMesh cut_mesh(mesh, n, vs);
    cut_mesh.construct(cut_t_ids);
    if (cut_mesh.snap_to_plane()) {
        cout<<"mesh.tets.size() = "<<mesh.tets.size()<<endl;
        cout << "cut_t_ids.size() " << cut_t_ids.size() << "->";
        cut_mesh.expand(cut_t_ids);
        cout << cut_t_ids.size() << " expanded" << endl;
        cout << "snapped #v = " << std::count(cut_mesh.is_snapped.begin(), cut_mesh.is_snapped.end(), true) << endl;
    }

    /////
    std::vector<Vector3> points;
    std::map<std::array<int, 2>, int> map_edge_to_intersecting_point;
    std::vector<int> subdivide_t_ids;
    if (!cut_mesh.get_intersecting_edges_and_points(points, map_edge_to_intersecting_point, subdivide_t_ids))
        return false;
    //have to add all cut_t_ids
    vector_unique(cut_t_ids);
    std::vector<int> tmp;
    std::set_difference(subdivide_t_ids.begin(), subdivide_t_ids.end(), cut_t_ids.begin(), cut_t_ids.end(),
                        std::back_inserter(tmp));
    std::vector<bool> is_mark_surface(cut_t_ids.size(), true);
    cut_t_ids.insert(cut_t_ids.end(), tmp.begin(), tmp.end());
    is_mark_surface.resize(is_mark_surface.size() + tmp.size(), false);
    cout << "cut_mesh.get_intersecting_edges_and_points OK" << endl;

    /////
    std::vector<MeshTet> new_tets;
    std::vector<std::array<std::vector<int>, 4>> new_track_surface_fs;
    std::vector<int> modified_t_ids;
    if (!subdivide_tets(insert_f_id, mesh, cut_mesh, points, map_edge_to_intersecting_point, track_surface_fs,
                        cut_t_ids, is_mark_surface,
                        new_tets, new_track_surface_fs, modified_t_ids))
        return false;

    //fortest
    cout << "subdivide_tets OK" << endl;
    cout << "cut_t_ids.size() = " << cut_t_ids.size() << endl;
    cout << "new_tets.size() = " << new_tets.size() << endl;
    myassert(new_tets.size() == new_track_surface_fs.size());
    cout << "points.size() = " << points.size() << endl;
    //fortest

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
                mesh.tets.push_back(new_tets[i]);
                track_surface_fs.push_back(new_track_surface_fs[i]);
                for (int j = 0; j < 4; j++) {
                    mesh.tet_vertices[mesh.tets.back()[j]].conn_tets.push_back(mesh.tets.size() - 1);
                }
            }
            //todo: tags???
        }
    } else {
        //todo
    }

    return true;
}

void floatTetWild::find_cutting_tets(int f_id, const std::vector<Vector3i> &input_faces,
                                     const std::array<Vector3, 3>& vs, Mesh &mesh, std::vector<int> &cut_t_ids) {
    for (int t_id = 0; t_id < mesh.tets.size(); t_id++) {
        std::array<int, 4> oris;
        for (int j = 0; j < 4; j++) {
            oris[j] = Predicates::orient_3d(vs[0], vs[1], vs[2], mesh.tet_vertices[mesh.tets[t_id][j]].pos);
        }

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
            if (cnt_on == 3) {
                result = is_tri_tri_cutted_hint(vs[0], vs[1], vs[2], tp1, tp2, tp3, CUT_COPLANAR);
            } else if (cnt_pos > 0 && cnt_neg > 0) {
                result = is_tri_tri_cutted_hint(vs[0], vs[1], vs[2], tp1, tp2, tp3, CUT_FACE);
            }
            if (result == CUT_EMPTY)
                continue;

            cut_t_ids.push_back(t_id);
            break;
        }
    }

    return;//fortest

    std::vector<int> n_t_ids;
    for (int j = 0; j < 3; j++) {
        n_t_ids.insert(n_t_ids.end(), mesh.tet_vertices[input_faces[f_id][j]].conn_tets.begin(),
                       mesh.tet_vertices[input_faces[f_id][j]].conn_tets.end());
    }
    vector_unique(n_t_ids);

    std::vector<bool> is_visited(mesh.tets.size(), false);
    std::queue<int> queue_t_ids;
    for (int t_id: n_t_ids)
        queue_t_ids.push(t_id);
    while (!queue_t_ids.empty()) {
        int t_id = queue_t_ids.front();
        queue_t_ids.pop();
        if (is_visited[t_id])
            continue;
        is_visited[t_id] = true;

        std::array<int, 4> oris;
        for (int j = 0; j < 4; j++) {
            oris[j] = Predicates::orient_3d(vs[0], vs[1], vs[2], mesh.tet_vertices[mesh.tets[t_id][j]].pos);
        }

        bool is_cutted = false;
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
            if (cnt_on == 3) {
                result = is_tri_tri_cutted_hint(vs[0], vs[1], vs[2], tp1, tp2, tp3, CUT_COPLANAR);
            } else if (cnt_pos > 0 && cnt_neg > 0) {
                result = is_tri_tri_cutted_hint(vs[0], vs[1], vs[2], tp1, tp2, tp3, CUT_FACE);
            }
            if (result == CUT_EMPTY)
                continue;

            is_cutted = true;
            int opp_t_id = get_opp_t_id(t_id, j, mesh);
            if (opp_t_id >= 0 && !is_visited[opp_t_id])
                queue_t_ids.push(opp_t_id);
        }
        if (is_cutted)
            cut_t_ids.push_back(t_id);
    }
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
                        myassert(cut_mesh.map_v_ids.find(mesh.tets[t_id][(j + k + 1) % 4]) != cut_mesh.map_v_ids.end());//fortest
                        if (cut_mesh.is_v_on_plane(cut_mesh.map_v_ids[mesh.tets[t_id][(j + k + 1) % 4]])) {
                            cnt_on++;
//                            cout << (j + k + 1) % 4 << endl;//fortest
                        }
                    }
//                    cout << "cnt_on = " << cnt_on << endl;//fortest
                    if (cnt_on == 3) {
                        new_tets.push_back(mesh.tets[t_id]);
                        new_track_surface_fs.push_back(track_surface_fs[t_id]);
                        (new_track_surface_fs.back())[j].push_back(insert_f_id);
                        modified_t_ids.push_back(t_id);

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
            Scalar min_q;
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
                        } else
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
                if (cnt == 0)
                    min_q = volume;
                else if (volume < min_q)
                    min_q = volume;
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
                if (min_q < SCALAR_ZERO_3)
                    continue;
                min_qualities[i] = std::make_pair(i, min_q);
                all_centroids[i] = tmp_centroids;
            }
            std::sort(min_qualities.begin(), min_qualities.end(),
                      [](const std::pair<int, Scalar> &a, const std::pair<int, Scalar> &b) {
                          return a.second < b.second;
                      });

            if (min_qualities.back().second < SCALAR_ZERO_3) // if tet quality is too bad
                return false;

            diag_config_id = min_qualities.back().first;
            centroids = all_centroids[diag_config_id];
        } else {
            Scalar min_q = check_config(diag_config_id, centroids);
            if (min_q < SCALAR_ZERO_3)
                return false;
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

void floatTetWild::mark_surface_fs(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
                     std::vector<std::array<std::vector<int>, 4>> &track_surface_fs, Mesh &mesh) {
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

    for (int t_id = 0; t_id < track_surface_fs.size(); t_id++) {
        for (int j = 0; j < 4; j++) {
//            //fortest
//            if (track_surface_fs[t_id][j].size() > 0)
////            if (std::find(track_surface_fs[t_id][j].begin(), track_surface_fs[t_id][j].end(), III) != track_surface_fs[t_id][j].end())
//                mesh.tets[t_id].is_surface_fs[j] = 1;
//            continue;
//            //fortest

            auto &f_ids = track_surface_fs[t_id][j];
            auto &tp1_3d = mesh.tet_vertices[mesh.tets[t_id][(j + 1) % 4]].pos;
            auto &tp2_3d = mesh.tet_vertices[mesh.tets[t_id][(j + 2) % 4]].pos;
            auto &tp3_3d = mesh.tet_vertices[mesh.tets[t_id][(j + 3) % 4]].pos;
            int t = get_t(tp1_3d, tp2_3d, tp3_3d);
            std::array<Vector2, 3> tps_2d = {{to_2d(tp1_3d, t), to_2d(tp2_3d, t), to_2d(tp3_3d, t)}};
            Vector2 c = (tps_2d[0] + tps_2d[1] + tps_2d[2]) / 3;
            for (int f_id: f_ids) {
                std::array<Vector2, 3> ps_2d = {{to_2d(input_vertices[input_faces[f_id][0]], t),
                                                        to_2d(input_vertices[input_faces[f_id][1]], t),
                                                        to_2d(input_vertices[input_faces[f_id][2]], t)}};
                if (!is_on_bounded_side(ps_2d, c))
                    continue;

                mesh.tets[t_id].is_surface_fs[j] = 1;//todo: track orientation!!
                break;
            }
        }
    }

    //fortest: output and check
    output_surface(mesh, "surface.stl");
}

int floatTetWild::get_opp_t_id(int t_id, int j, const Mesh &mesh){
    std::vector<int> tmp;
    set_intersection(mesh.tet_vertices[mesh.tets[t_id][(j + 1) % 4]].conn_tets,
                     mesh.tet_vertices[mesh.tets[t_id][(j + 2) % 4]].conn_tets,
                     mesh.tet_vertices[mesh.tets[t_id][(j + 3) % 4]].conn_tets,
                     tmp);
    //fortest
    if(tmp.size() != 1 && tmp.size() != 2) {
        cout << "tmp.size() = " << tmp.size() << endl;
        cout << "t_id = " << t_id << ", j = " << j << endl;
        vector_print(mesh.tet_vertices[mesh.tets[t_id][(j + 1) % 4]].conn_tets);
        vector_print(mesh.tet_vertices[mesh.tets[t_id][(j + 2) % 4]].conn_tets);
        vector_print(mesh.tet_vertices[mesh.tets[t_id][(j + 3) % 4]].conn_tets);
        vector_print(tmp);
        for(int i: tmp){
            mesh.tets[i].print();
        }
        pausee();
    }
    //fortest
    if (tmp.size() == 2)
        return tmp[0] == t_id ? tmp[1] : tmp[0];
    else
        return -1;
}

void floatTetWild::myassert(bool b) {
    if (b == false) {
        cout<<"myassert fail!"<<endl;
        pausee();
    }
}

void floatTetWild::check_track_surface_fs(Mesh &mesh, std::vector<std::array<std::vector<int>, 4>> &track_surface_fs){
    //check connection
    std::vector<std::vector<int>> conn_tets(mesh.tet_vertices.size());
    for (int i = 0; i < mesh.tets.size(); i++) {
        for (int j = 0; j < 4; j++)
            conn_tets[mesh.tets[i][j]].push_back(i);
    }
    for (int i = 0; i < mesh.tet_vertices.size(); i++) {
        std::sort(mesh.tet_vertices[i].conn_tets.begin(), mesh.tet_vertices[i].conn_tets.end());
        if (mesh.tet_vertices[i].conn_tets != conn_tets[i]) {
            cout << "mesh.tet_vertices[i].conn_tets!=conn_tets[i]" << endl;
            pausee();
        }
    }
    cout<<"check 1 done"<<endl;

    for (int i = 0; i < mesh.tets.size(); i++) {
        for (int j = 0; j < 4; j++) {
            int opp_t_id = get_opp_t_id(i, j, mesh);
        }
    }
    cout<<"check 2 done"<<endl;

    for (int t_id = 0; t_id < track_surface_fs.size(); t_id++) {
        for (int j = 0; j < 4; j++) {
            std::sort(track_surface_fs[t_id][j].begin(), track_surface_fs[t_id][j].end());
            int opp_t_id = get_opp_t_id(t_id, j, mesh);
            if (opp_t_id < 0) {
                if (!track_surface_fs[t_id][j].empty()) {
                    cout << "bbox face !track_surface_fs[t_id][j].empty()" << endl;
                    pausee();
                }
                continue;
            }
            int k = 0;
            for (k = 0; k < 4; k++) {
                if (mesh.tets[opp_t_id][k] != mesh.tets[t_id][(j + 1) % 4]
                    && mesh.tets[opp_t_id][k] != mesh.tets[t_id][(j + 2) % 4]
                    && mesh.tets[opp_t_id][k] != mesh.tets[t_id][(j + 3) % 4])
                    break;
            }
            std::sort(track_surface_fs[opp_t_id][k].begin(), track_surface_fs[opp_t_id][k].end());
            if (track_surface_fs[t_id][j] != track_surface_fs[opp_t_id][k]) {
                cout << "track_surface_fs[t_id][j]!=track_surface_fs[opp_t_id][k]" << endl;
                cout << "t " << t_id << ": ";
                vector_print(track_surface_fs[t_id][j]);
                cout << "opp_t " << opp_t_id << ": ";
                vector_print(track_surface_fs[opp_t_id][k]);
                pausee();
            }
        }
    }
    cout<<"check 3 done"<<endl;
}