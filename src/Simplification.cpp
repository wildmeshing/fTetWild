// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include <floattetwild/Simplification.h>
#include <floattetwild/Logger.hpp>
#include <floattetwild/LocalOperations.h>

#include <igl/remove_duplicate_vertices.h>
#include <igl/writeOFF.h>
#include <igl/Timer.h>
#include <igl/unique_rows.h>

#ifdef FLOAT_TETWILD_USE_TBB
#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include <tbb/atomic.h>
#include <tbb/concurrent_vector.h>
#include <tbb/parallel_sort.h>
#include <tbb/concurrent_unordered_set.h>
#endif

void floatTetWild::simplify(std::vector<Vector3>& input_vertices, std::vector<Vector3i>& input_faces, std::vector<int>& input_tags,
        const AABBWrapper& tree, const Parameters& params, bool skip_simplify) {

    remove_duplicates(input_vertices, input_faces, input_tags, params);
    if (skip_simplify)
        return;

    std::vector<bool> v_is_removed(input_vertices.size(), false);
    std::vector<bool> f_is_removed(input_faces.size(), false);
    std::vector<std::unordered_set<int>> conn_fs(input_vertices.size());
    for (int i = 0; i < input_faces.size(); i++) {
        for (int j = 0; j < 3; j++)
            conn_fs[input_faces[i][j]].insert(i);
    }

    igl::Timer timer;
    timer.start();
    collapsing(input_vertices, input_faces, tree, params, v_is_removed, f_is_removed, conn_fs);
    std::cout<<"collapsing "<<timer.getElapsedTime()<<std::endl;

    timer.start();
    swapping(input_vertices, input_faces, tree, params, v_is_removed, f_is_removed, conn_fs);
    std::cout<<"swapping "<<timer.getElapsedTime()<<std::endl;

    //clean up vs, fs
    //v
    std::vector<int> map_v_ids(input_vertices.size(), -1);
    int cnt = 0;
    for (int i = 0; i < input_vertices.size(); i++) {
        if (v_is_removed[i] || conn_fs[i].empty())
            continue;
        map_v_ids[i] = cnt;
        cnt++;
    }

    std::vector<Vector3> new_input_vertices(cnt);
    cnt = 0;
    for (int i = 0; i < input_vertices.size(); i++) {
        if (v_is_removed[i] || conn_fs[i].empty())
            continue;
        new_input_vertices[cnt++] = input_vertices[i];
    }
    input_vertices = new_input_vertices;

    //f
    cnt = 0;
    for (int i = 0; i < input_faces.size(); i++) {
        if (f_is_removed[i])
            continue;
        for (int j = 0; j < 3; j++)
            input_faces[i][j] = map_v_ids[input_faces[i][j]];
        cnt++;
    }

    std::vector<Vector3i> new_input_faces(cnt);
    std::vector<int> new_input_tags(cnt);
    cnt = 0;
    for (int i = 0; i < input_faces.size(); i++) {
        if (f_is_removed[i])
            continue;
        new_input_faces[cnt] = input_faces[i];
        new_input_tags[cnt] = input_tags[i];
        cnt++;
    }
    input_faces = new_input_faces;
    input_tags = new_input_tags;

//    flattening(input_vertices, input_faces, tree, params);

    remove_duplicates(input_vertices, input_faces, input_tags, params);

    logger().info("#v = {}", input_vertices.size());
    logger().info("#f = {}", input_faces.size());

    ////////////////////////
    //output
    if(params.log_level < 3) {
        Eigen::MatrixXd V(input_vertices.size(), 3);
        Eigen::MatrixXi F(input_faces.size(), 3);
        for (int i = 0; i < input_vertices.size(); i++) {
            V.row(i) = input_vertices[i];
        }
        for (int i = 0; i < input_faces.size(); i++) {
            F.row(i) = input_faces[i];
        }
        if (!params.output_path.empty()) {
            igl::writeOFF(params.output_path + "_" + params.postfix + "_simplify.off", V, F);
        }
    }

//    ////////////////////////
//    //check
//    auto tmp = input_vertices;
//    std::sort(tmp.begin(), tmp.end(), [](const Vector3& a, const Vector3& b){
//       return std::tuple<Scalar, Scalar, Scalar>(a[0], a[1], a[2]) < std::tuple<Scalar, Scalar, Scalar>(b[0], b[1], b[2]);
//    });
//    for(int i=0;i<tmp.size()-1;i++){
//        if(tmp[i]==tmp[i+1]){
//            cout<<"find duplicates! "<<i<<endl;
//        }
//    }
}

bool floatTetWild::remove_duplicates(std::vector<Vector3>& input_vertices, std::vector<Vector3i>& input_faces, std::vector<int>& input_tags, const Parameters& params) {
//    std::vector<size_t> indices(input_vertices.size());
//    for(size_t i=0;i<input_vertices.size();i++)
//        indices[i] = i;
//
//    std::sort(indices.begin(), indices.end(), [&input_vertices](size_t i1, size_t i2) {
//        return std::make_tuple(input_vertices[i1][0], input_vertices[i1][1], input_vertices[i1][2])
//               < std::make_tuple(input_vertices[i2][0], input_vertices[i2][1], input_vertices[i2][2]);
//    });
//    indices.erase(std::unique(indices.begin(), indices.end(), [&input_vertices](size_t i1, size_t i2) {
//        return std::make_tuple(input_vertices[i1][0], input_vertices[i1][1], input_vertices[i1][2])
//               == std::make_tuple(input_vertices[i2][0], input_vertices[i2][1], input_vertices[i2][2]);
//    }), indices.end());

    MatrixXs V_tmp(input_vertices.size(), 3), V_in;
    Eigen::MatrixXi F_tmp(input_faces.size(), 3), F_in;
    for (int i = 0; i < input_vertices.size(); i++)
        V_tmp.row(i) = input_vertices[i];
    for (int i = 0; i < input_faces.size(); i++)
        F_tmp.row(i) = input_faces[i];

    //
    Eigen::VectorXi IV, _;
    igl::remove_duplicate_vertices(V_tmp, F_tmp, SCALAR_ZERO * params.bbox_diag_length, V_in, IV, _, F_in);
    //
    for (int i = 0; i < F_in.rows(); i++) {
        int j_min = 0;
        for (int j = 1; j < 3; j++) {
            if (F_in(i, j) < F_in(i, j_min))
                j_min = j;
        }
        if (j_min == 0)
            continue;
        int v0_id = F_in(i, j_min);
        int v1_id = F_in(i, (j_min + 1) % 3);
        int v2_id = F_in(i, (j_min + 2) % 3);
        F_in.row(i) << v0_id, v1_id, v2_id;
    }
    F_tmp.resize(0, 0);
    Eigen::VectorXi IF;
    igl::unique_rows(F_in, F_tmp, IF, _);
    F_in = F_tmp;
    std::vector<int> old_input_tags = input_tags;
    input_tags.resize(IF.rows());
    for (int i = 0; i < IF.rows(); i++) {
        input_tags[i] = old_input_tags[IF(i)];
    }
    //
    if (V_in.rows() == 0 || F_in.rows() == 0)
        return false;

    logger().info("remove duplicates: ");
    logger().info("#v: {} -> {}", input_vertices.size(), V_in.rows());
    logger().info("#f: {} -> {}", input_faces.size(), F_in.rows());

    input_vertices.resize(V_in.rows());
    input_faces.clear();
    input_faces.reserve(F_in.rows());
    old_input_tags = input_tags;
    input_tags.clear();
    for (int i = 0; i < V_in.rows(); i++)
        input_vertices[i] = V_in.row(i);
    for (int i = 0; i < F_in.rows(); i++) {
        if (F_in(i, 0) == F_in(i, 1) || F_in(i, 0) == F_in(i, 2) || F_in(i, 2) == F_in(i, 1))
            continue;
        if (i > 0 && (F_in(i, 0) == F_in(i - 1, 0) && F_in(i, 1) == F_in(i - 1, 2) && F_in(i, 2) == F_in(i - 1, 1)))
            continue;
        //check area
        Vector3 u = V_in.row(F_in(i, 1)) - V_in.row(F_in(i, 0));
        Vector3 v = V_in.row(F_in(i, 2)) - V_in.row(F_in(i, 0));
        Vector3 area = u.cross(v);
        if (area.norm() / 2 <= SCALAR_ZERO * params.bbox_diag_length)
            continue;
        input_faces.push_back(F_in.row(i));
        input_tags.push_back(old_input_tags[i]);
    }

    return true;
}

void floatTetWild::collapsing(std::vector<Vector3>& input_vertices, std::vector<Vector3i>& input_faces,
        const AABBWrapper& tree, const Parameters& params,
        std::vector<bool>& v_is_removed, std::vector<bool>& f_is_removed, std::vector<std::unordered_set<int>>& conn_fs){

#ifdef FLOAT_TETWILD_USE_TBB
    std::vector<std::array<int, 2>> edges;
    tbb::concurrent_vector<std::array<int, 2>> edges_tbb;

    const auto build_edges = [&](){
        edges.clear();
        edges.reserve(input_faces.size()*3);

        edges_tbb.clear();
        edges_tbb.reserve(input_faces.size()*3);

        tbb::parallel_for( size_t(0), input_faces.size(), [&](size_t f_id)
        {
            if(f_is_removed[f_id])
                return;

            for (int j = 0; j < 3; j++) {
                std::array<int, 2> e = {{input_faces[f_id][j], input_faces[f_id][mod3(j + 1)]}};
                if (e[0] > e[1])
                    std::swap(e[0], e[1]);
                edges_tbb.push_back(e);
            }
        });


        edges.reserve(edges_tbb.size());
        edges.insert(edges.end(), edges_tbb.begin(), edges_tbb.end());
        assert(edges_tbb.size() == edges.size());
        tbb::parallel_sort(edges.begin(), edges.end());

        edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    };
#else
    std::vector<std::array<int, 2>> edges;
    edges.clear();
    edges.reserve(input_faces.size()*3);
    for(size_t f_id = 0; f_id < input_faces.size(); ++f_id)
    {
        if(f_is_removed[f_id])
            continue;

        const auto& f = input_faces[f_id];
        for (int j = 0; j < 3; j++) {
            std::array<int, 2> e = {{f[j], f[mod3(j + 1)]}};
            if (e[0] > e[1])
                std::swap(e[0], e[1]);
            edges.push_back(e);
        }
    }
    vector_unique(edges);

    std::priority_queue<ElementInQueue, std::vector<ElementInQueue>, cmp_s> sm_queue;
    for (auto& e: edges) {
        Scalar weight = (input_vertices[e[0]] - input_vertices[e[1]]).squaredNorm();
        sm_queue.push(ElementInQueue(e, weight));
        sm_queue.push(ElementInQueue(std::array<int, 2>({{e[1], e[0]}}), weight));
    }
#endif

    auto is_onering_clean = [&](int v_id) {
        std::vector<int> v_ids;
        v_ids.reserve(conn_fs[v_id].size() * 2);
        for (const auto &f_id:conn_fs[v_id]) {
            for (int j = 0; j < 3; j++) {
                if (input_faces[f_id][j] != v_id)
                    v_ids.push_back(input_faces[f_id][j]);
            }
        }
        std::sort(v_ids.begin(), v_ids.end());

        if (v_ids.size() % 2 != 0)
            return false;
        for (int i = 0; i < v_ids.size(); i += 2) {
            if (v_ids[i] != v_ids[i + 1])
                return false;
        }

        return true;
    };

    static const int SUC = 1;
    static const int FAIL_CLEAN = 0;
    static const int FAIL_FLIP = -1;
    static const int FAIL_ENV = -2;

    auto remove_an_edge = [&](int v1_id, int v2_id, const std::vector<int>& n12_f_ids) {
        if (!is_onering_clean(v1_id) || !is_onering_clean(v2_id))
            return FAIL_CLEAN;

//        std::unordered_set<int> new_f_ids;
        std::vector<int> new_f_ids;
        for (int f_id:conn_fs[v1_id]) {
            if (f_id != n12_f_ids[0] && f_id != n12_f_ids[1])
                new_f_ids.push_back(f_id);
        }
        for (int f_id:conn_fs[v2_id]) {
            if (f_id != n12_f_ids[0] && f_id != n12_f_ids[1])
                new_f_ids.push_back(f_id);
        }
        vector_unique(new_f_ids);

        //compute new point
        Vector3 p = (input_vertices[v1_id] + input_vertices[v2_id]) / 2;
        tree.project_to_sf(p);
        // GEO::vec3 geo_p(p[0], p[1], p[2]);
        // GEO::vec3 nearest_p;
        // double _;
        // tree.nearest_facet(geo_p, nearest_p, _);
        // p[0] = nearest_p[0];
        // p[1] = nearest_p[1];
        // p[2] = nearest_p[2];

        //computing normal for checking flipping
        for (int f_id:new_f_ids) {
            Vector3 old_nv = (input_vertices[input_faces[f_id][1]] - input_vertices[input_faces[f_id][2]]).cross(input_vertices[input_faces[f_id][0]] - input_vertices[input_faces[f_id][2]]);

            for (int j = 0; j < 3; j++) {
                if (input_faces[f_id][j] == v1_id || input_faces[f_id][j] == v2_id) {
                    Vector3 new_nv = (input_vertices[input_faces[f_id][mod3(j + 1)]] - input_vertices[input_faces[f_id][mod3(j + 2)]]).cross(p - input_vertices[input_faces[f_id][mod3(j + 2)]]);
                    if (old_nv.dot(new_nv) <= 0)
                        return FAIL_FLIP;
                    //check new tris' area
                    if (new_nv.norm() / 2 <= SCALAR_ZERO_2)
                        return FAIL_FLIP;
                    break;
                }
            }
        }

        //check if go outside of envelop
        for (int f_id:new_f_ids) {
            for (int j = 0; j < 3; j++) {
                if (input_faces[f_id][j] == v1_id || input_faces[f_id][j] == v2_id) {
                    const std::array<Vector3, 3> tri = {{
                        p,
                        input_vertices[input_faces[f_id][mod3(j + 1)]],
                        input_vertices[input_faces[f_id][mod3(j + 2)]]
                    }};
                    if (is_out_envelope(tri, tree, params))
                        return FAIL_ENV;
                    break;
                }
            }
        }

        //real update
//        std::unordered_set<int> n_v_ids;//get this info before real update for later usage
        std::vector<int> n_v_ids;//get this info before real update for later usage
        for (int f_id:new_f_ids) {
            for (int j = 0; j < 3; j++) {
                if (input_faces[f_id][j] != v1_id && input_faces[f_id][j] != v2_id)
                    n_v_ids.push_back(input_faces[f_id][j]);
            }
        }
        vector_unique(n_v_ids);

        v_is_removed[v1_id] = true;
        input_vertices[v2_id] = p;
        for (int f_id:n12_f_ids) {
            f_is_removed[f_id] = true;
#ifndef FLOAT_TETWILD_USE_TBB
            for (int j = 0; j < 3; j++) {//rm conn_fs
                if (input_faces[f_id][j] != v1_id) {
                    conn_fs[input_faces[f_id][j]].erase(f_id);
                }
            }
#endif
        }
        for (int f_id:conn_fs[v1_id]) {//add conn_fs
            if (f_is_removed[f_id])
                continue;
            conn_fs[v2_id].insert(f_id);
            for (int j = 0; j < 3; j++) {
                if (input_faces[f_id][j] == v1_id)
                    input_faces[f_id][j] = v2_id;
            }
        }

#ifndef FLOAT_TETWILD_USE_TBB
        //push new edges into the queue
        for (int v_id:n_v_ids) {
            double weight = (input_vertices[v2_id] - input_vertices[v_id]).squaredNorm();
            sm_queue.push(ElementInQueue(std::array<int, 2>({{v2_id, v_id}}), weight));
            sm_queue.push(ElementInQueue(std::array<int, 2>({{v_id, v2_id}}), weight));
        }
#endif
        return SUC;
    };

#ifdef FLOAT_TETWILD_USE_TBB
    tbb::atomic<int> cnt(0);
    int cnt_suc = 0;
    // tbb::atomic<int> fail_clean(0);
    // tbb::atomic<int> fail_flip(0);
    // tbb::atomic<int> fail_env(0);

    const int stopping = input_vertices.size()/10000.;

    std::vector<int> safe_set;
    do {
        build_edges();
        Mesh::one_ring_edge_set(edges, v_is_removed, f_is_removed, conn_fs, input_vertices, safe_set);
        cnt = 0;

        tbb::parallel_for(size_t(0), safe_set.size(), [&](size_t i) {
//        for (int i = 0; i < safe_set.size(); i++) {
            std::array<int, 2> &v_ids = edges[safe_set[i]];

//            if (v_is_removed[v_ids[0]] || v_is_removed[v_ids[1]])
//                return;

            std::vector<int> n12_f_ids;
            set_intersection(conn_fs[v_ids[0]], conn_fs[v_ids[1]], n12_f_ids);

            if (n12_f_ids.size() != 2)
                return;
//                continue;

            int res = remove_an_edge(v_ids[0], v_ids[1], n12_f_ids);
            if (res == SUC)
                cnt++;
        });
//        }

        //cleanup conn_fs
        tbb::parallel_for(size_t(0), conn_fs.size(), [&](size_t i) {
//        for (int i = 0; i < conn_fs.size(); i++) {
            if (v_is_removed[i])
//                continue;
                return;
            std::vector<int> r_f_ids;
            for (int f_id: conn_fs[i]) {
                if (f_is_removed[f_id])
                    r_f_ids.push_back(f_id);
            }
            for (int f_id:r_f_ids)
                conn_fs[i].erase(f_id);
//        }
        });

        cnt_suc += cnt;
    } while(cnt > stopping);

#else
    int cnt_suc = 0;
    int fail_clean = 0;
    int fail_flip = 0;
    int fail_env = 0;

    while (!sm_queue.empty()) {
        std::array<int, 2> v_ids = sm_queue.top().v_ids;
        Scalar old_weight = sm_queue.top().weight;
        sm_queue.pop();

        if (v_is_removed[v_ids[0]] || v_is_removed[v_ids[1]])
            continue;
        if(old_weight!=(input_vertices[v_ids[0]] - input_vertices[v_ids[1]]).squaredNorm())
            continue;

        std::vector<int> n12_f_ids;
        set_intersection(conn_fs[v_ids[0]], conn_fs[v_ids[1]], n12_f_ids);
        if(n12_f_ids.size()!=2)
            continue;

        int res = remove_an_edge(v_ids[0], v_ids[1], n12_f_ids);
        if (res == SUC)
            cnt_suc++;
        //        else if (res == FAIL_CLEAN)
        //            fail_clean++;
        //        else if (res == FAIL_FLIP)
        //            fail_flip++;
        //        else if (res == FAIL_ENV)
        //            fail_env++;
    }
#endif

    //    cout<<fail_clean<<endl;
    //    cout<<fail_flip<<endl;
    //    cout<<fail_env<<endl;
//    std::cout<<"#v: "<<build_time<<std::endl;
    logger().debug("{}  faces are collapsed!!", cnt_suc);
}

void floatTetWild::swapping(std::vector<Vector3>& input_vertices, std::vector<Vector3i>& input_faces,
        const AABBWrapper& tree, const Parameters& params,
        std::vector<bool>& v_is_removed, std::vector<bool>& f_is_removed, std::vector<std::unordered_set<int>>& conn_fs) {
    std::vector<std::array<int, 2>> edges;
    edges.reserve(input_faces.size() * 6);
    for (int i = 0; i < input_faces.size(); i++) {
        if (f_is_removed[i])
            continue;
        auto &f = input_faces[i];
        for (int j = 0; j < 3; j++) {
            std::array<int, 2> e = {{f[j], f[mod3(j + 1)]}};
            if (e[0] > e[1])
                std::swap(e[0], e[1]);
            edges.push_back(e);
        }
    }
    vector_unique(edges);

    std::priority_queue<ElementInQueue, std::vector<ElementInQueue>, cmp_l> sm_queue;
    for (auto &e: edges) {
        Scalar weight = (input_vertices[e[0]] - input_vertices[e[1]]).squaredNorm();
        sm_queue.push(ElementInQueue(e, weight));
        sm_queue.push(ElementInQueue(std::array<int, 2>({{e[1], e[0]}}), weight));
    }

    int cnt = 0;
    while (!sm_queue.empty()) {
        int v1_id = sm_queue.top().v_ids[0];
        int v2_id = sm_queue.top().v_ids[1];
        sm_queue.pop();

        std::vector<int> n12_f_ids;
        set_intersection(conn_fs[v1_id], conn_fs[v2_id], n12_f_ids);
        if (n12_f_ids.size() != 2)
            continue;

        std::array<int, 2> n_v_ids = {{-1, -1}};
        for (int j = 0; j < 3; j++) {
            if (n_v_ids[0] < 0
                && input_faces[n12_f_ids[0]][j] != v1_id && input_faces[n12_f_ids[0]][j] != v2_id)
                n_v_ids[0] = input_faces[n12_f_ids[0]][j];

            if (n_v_ids[1] < 0
                && input_faces[n12_f_ids[1]][j] != v1_id && input_faces[n12_f_ids[1]][j] != v2_id)
                n_v_ids[1] = input_faces[n12_f_ids[1]][j];
        }

        //check coplanar
        Scalar cos_a0 = get_angle_cos(input_vertices[n_v_ids[0]], input_vertices[v1_id], input_vertices[v2_id]);
        Scalar cos_a1 = get_angle_cos(input_vertices[n_v_ids[1]], input_vertices[v1_id], input_vertices[v2_id]);
        std::array<Vector3, 2> old_nvs;
        for (int f = 0; f < 2; f++) {
            auto &a = input_vertices[input_faces[n12_f_ids[f]][0]];
            auto &b = input_vertices[input_faces[n12_f_ids[f]][1]];
            auto &c = input_vertices[input_faces[n12_f_ids[f]][2]];
            old_nvs[f] = ((b - c).cross(a - c)).normalized();
        }
        if (cos_a0 > -0.999) {//maybe it's for avoiding numerical issue
            if (old_nvs[0].dot(old_nvs[1]) < 1 - 1e-6)//not coplanar
                continue;
        }

        //check inversion
        auto &old_nv = cos_a1 < cos_a0 ? old_nvs[0] : old_nvs[1];
        bool is_filp = false;
        for (int f_id:n12_f_ids) {
            auto &a = input_vertices[input_faces[f_id][0]];
            auto &b = input_vertices[input_faces[f_id][1]];
            auto &c = input_vertices[input_faces[f_id][2]];
            if (old_nv.dot(((b - c).cross(a - c)).normalized()) < 0) {
                is_filp = true;
                break;
            }
        }
        if (is_filp)
            continue;

        //check quality
        Scalar cos_a0_new = get_angle_cos(input_vertices[v1_id], input_vertices[n_v_ids[0]],
                                          input_vertices[n_v_ids[1]]);
        Scalar cos_a1_new = get_angle_cos(input_vertices[v2_id], input_vertices[n_v_ids[0]],
                                          input_vertices[n_v_ids[1]]);
        if (std::min(cos_a0_new, cos_a1_new) <= std::min(cos_a0, cos_a1))
            continue;

        //check envelope
//        bool is_valid = true;
//        for(int v_id: n_v_ids) {
//            if (is_out_envelope({{input_vertices[v_id], input_vertices[v1_id], input_vertices[v2_id]}},
//                                tree, params)) {
//                is_valid = false;
//                break;
//            }
//        }
//        if(!is_valid)
//            continue;
        if (is_out_envelope({{input_vertices[v1_id], input_vertices[n_v_ids[0]], input_vertices[n_v_ids[1]]}}, tree,
                            params)
            || is_out_envelope({{input_vertices[v2_id], input_vertices[n_v_ids[0]], input_vertices[n_v_ids[1]]}}, tree,
                               params)) {
            continue;
        }

        // real update
        for (int j = 0; j < 3; j++) {
            if (input_faces[n12_f_ids[0]][j] == v2_id)
                input_faces[n12_f_ids[0]][j] = n_v_ids[1];
            if (input_faces[n12_f_ids[1]][j] == v1_id)
                input_faces[n12_f_ids[1]][j] = n_v_ids[0];
        }
        conn_fs[v1_id].erase(n12_f_ids[1]);
        conn_fs[v2_id].erase(n12_f_ids[0]);
        conn_fs[n_v_ids[0]].insert(n12_f_ids[1]);
        conn_fs[n_v_ids[1]].insert(n12_f_ids[0]);
        cnt++;

//        cout << v1_id << " " << v2_id << endl;
//        cout << "neighbors " << n_v_ids[0] << " " << n_v_ids[1] << endl;
//        check_surface(input_vertices, input_faces, f_is_removed, tree, params);
    }

    logger().debug("{}  faces are swapped!!", cnt);
    return;

    ///////////////////

    for (int i = 0; i < input_faces.size(); i++) {//todo go over edges!!!
        if(f_is_removed[i])
            continue;

        bool is_swapped = false;
        for (int j = 0; j < 3; j++) {
            int v_id = input_faces[i][j];
            int v1_id = input_faces[i][mod3(j+1)];
            int v2_id = input_faces[i][mod3(j+2)];

            // manifold
            std::vector<int> n12_f_ids;
            set_intersection(conn_fs[v1_id], conn_fs[v2_id], n12_f_ids);
            if (n12_f_ids.size() != 2) {
                continue;
            }
            if (n12_f_ids[1] == i)
                n12_f_ids = {n12_f_ids[1], n12_f_ids[0]};
            int v3_id = -1;
            for (int k = 0; k < 3; k++)
                if (input_faces[n12_f_ids[1]][k] != v1_id && input_faces[n12_f_ids[1]][k] != v2_id) {
                    v3_id = input_faces[n12_f_ids[1]][k];
                    break;
                }

            // check quality
            Scalar cos_a = get_angle_cos(input_vertices[v_id], input_vertices[v1_id], input_vertices[v2_id]);
            Scalar cos_a1 = get_angle_cos(input_vertices[v3_id], input_vertices[v1_id], input_vertices[v2_id]);
            std::array<Vector3, 2> old_nvs;
            for (int f = 0; f < 2; f++) {
                auto& a = input_vertices[input_faces[n12_f_ids[f]][0]];
                auto& b = input_vertices[input_faces[n12_f_ids[f]][1]];
                auto& c = input_vertices[input_faces[n12_f_ids[f]][2]];
                old_nvs[f] = ((b - c).cross(a - c)).normalized();
            }
            if (cos_a > -0.999) {
//                continue;
                if (old_nvs[0].dot(old_nvs[1]) < 1-1e-6)//not coplanar
                    continue;
            }
            Scalar cos_a_new = get_angle_cos(input_vertices[v1_id], input_vertices[v_id], input_vertices[v3_id]);
            Scalar cos_a1_new = get_angle_cos(input_vertices[v2_id], input_vertices[v_id], input_vertices[v3_id]);
            if (std::min(cos_a_new, cos_a1_new) <= std::min(cos_a, cos_a1))
                continue;

            // non flipping
            auto f1_old = input_faces[n12_f_ids[0]];
            auto f2_old = input_faces[n12_f_ids[1]];
            for (int k = 0; k < 3; k++) {
                if (input_faces[n12_f_ids[0]][k] == v2_id)
                    input_faces[n12_f_ids[0]][k] = v3_id;
                if (input_faces[n12_f_ids[1]][k] == v1_id)
                    input_faces[n12_f_ids[1]][k] = v_id;
            }
            auto& old_nv = cos_a1 < cos_a ? old_nvs[0] : old_nvs[1];
            bool is_filp = false;
            for (int f_id:n12_f_ids) {
                auto& a = input_vertices[input_faces[f_id][0]];
                auto& b = input_vertices[input_faces[f_id][1]];
                auto& c = input_vertices[input_faces[f_id][2]];
                if (old_nv.dot(((b - c).cross(a - c)).normalized()) < 0) {
                    is_filp = true;
                    break;
                }
            }
            if (is_filp) {
                input_faces[n12_f_ids[0]] = f1_old;
                input_faces[n12_f_ids[1]] = f2_old;
                continue;
            }

            // non outside envelop
            bool is_valid = true;
            for(int f_id: n12_f_ids) {
                if (is_out_envelope({{input_vertices[input_faces[f_id][0]], input_vertices[input_faces[f_id][1]],
                                     input_vertices[input_faces[f_id][2]]}}, tree, params)) {
                    is_valid = false;
                    break;
                }
            }
            if(!is_valid){
                input_faces[n12_f_ids[0]] = f1_old;
                input_faces[n12_f_ids[1]] = f2_old;
                continue;
            }

            // real update
            conn_fs[v1_id].erase(n12_f_ids[1]);
            conn_fs[v2_id].erase(n12_f_ids[0]);
            conn_fs[v_id].insert(n12_f_ids[1]);
            conn_fs[v3_id].insert(n12_f_ids[0]);
            is_swapped = true;
            break;
        }
        if (is_swapped)
            cnt++;
    }
    logger().debug("{}  faces are swapped!!", cnt);
}

void floatTetWild::flattening(std::vector<Vector3>& input_vertices, std::vector<Vector3i>& input_faces,
        const AABBWrapper& sf_tree, const Parameters& params) {
    std::vector<Vector3> ns(input_faces.size());
    for (int i = 0; i < input_faces.size(); i++) {
        ns[i] = ((input_vertices[input_faces[i][2]] - input_vertices[input_faces[i][0]]).cross(
                input_vertices[input_faces[i][1]] - input_vertices[input_faces[i][0]])).normalized();
    }

    std::vector<std::vector<int>> conn_fs(input_vertices.size());
    std::vector<std::array<int, 2>> edges;
    edges.reserve(input_faces.size() * 6);
    for (int i = 0; i < input_faces.size(); i++) {
        auto &f = input_faces[i];
        for (int j = 0; j < 3; j++) {
            conn_fs[f[j]].push_back(i);
            std::array<int, 2> e = {{f[j], f[mod3(j + 1)]}};
            if (e[0] > e[1])
                std::swap(e[0], e[1]);
            edges.push_back(e);
        }
    }
    vector_unique(edges);

    auto needs_flattening = [](const Vector3 &n1, const Vector3 &n2) {
        if(n1.dot(n2)>0.98) {
            cout << std::setprecision(17) << n1.dot(n2) << endl;
            cout << n1.norm() << " " << n2.norm() << endl;
        }
        return true;

        double d = std::abs(n1.dot(n2) - 1);
        cout<<n1.dot(n2)<<endl;
        if (d > 1e-15 && d < 1e-5)
            return true;
        return false;
    };

    std::vector<int> f_tss(input_faces.size(), 0);
    int ts = 0;
    std::queue<std::array<int, 5>> edge_queue;
    for (auto &e: edges) {
        std::vector<int> n_f_ids;
        set_intersection(conn_fs[e[0]], conn_fs[e[1]], n_f_ids);
        if (n_f_ids.size() != 2)
            continue;
        if (!needs_flattening(ns[n_f_ids[0]], ns[n_f_ids[1]]))
            continue;
        edge_queue.push({{e[0], e[1], n_f_ids[0], n_f_ids[1], ts}});
    }

    while (!edge_queue.empty()) {
        std::array<int, 2> e = {{edge_queue.front()[0], edge_queue.front()[1]}};
        std::array<int, 2> n_f_ids = {{edge_queue.front()[2], edge_queue.front()[3]}};
        int e_ts = edge_queue.front().back();
        edge_queue.pop();

        std::vector<int> all_f_ids;
        for (int j = 0; j < 2; j++)
            all_f_ids.insert(all_f_ids.end(), conn_fs[e[j]].begin(), conn_fs[e[j]].end());
        vector_unique(all_f_ids);

        bool is_valid = true;
        for(int f_id: all_f_ids){
            if(f_tss[f_id]>e_ts) {
                is_valid = false;
                break;
            }
        }
        if(!is_valid)
            continue;

        auto &n1 = ns[n_f_ids[0]];
        auto &n2 = ns[n_f_ids[1]];
//        if(n_f_ids[0] == 61 && n_f_ids[0] == 62 || n_f_ids[0] == 62 && n_f_ids[0] == 61){
//            cout<<n1<<endl;
//            cout<<n2<<endl;
//            cout<<n1.dot(n2)<<endl;
//            pausee();
//        }
        std::vector<int> n_v_ids;
        for (int f_id: n_f_ids) {
            for (int j = 0; j < 3; j++) {
                if (input_faces[f_id][j] != e[0] && input_faces[f_id][j] != e[1]) {
                    n_v_ids.push_back(input_faces[f_id][j]);
                    break;
                }
            }
        }
        assert(n_v_ids.size() == 2 && n_v_ids[0] != n_v_ids[1]);
        Vector3 n = (n1 + n2) / 2;
//        Vector3 p = (input_vertices[e[0]] + input_vertices[e[1]]) / 2;
        Vector3 p = (input_vertices[n_v_ids[0]] + input_vertices[n_v_ids[1]]) / 2;

        std::array<Vector3, 2> old_ps;
        for (int j = 0; j < 2; j++) {
            old_ps[j] = input_vertices[e[j]];
            input_vertices[e[j]] -= n.dot(input_vertices[e[j]] - p) * n;
        }
//        cout<<"ok1"<<endl;

        is_valid = true;
        for (int f_id: all_f_ids) {
            const std::array<Vector3, 3> tri = {{input_vertices[input_faces[f_id][0]],
                                                        input_vertices[input_faces[f_id][1]],
                                                        input_vertices[input_faces[f_id][2]]}};
            if (is_out_envelope(tri, sf_tree, params)) {
                is_valid = false;
                break;
            }
        }
        if (!is_valid) {
            for (int j = 0; j < 2; j++)
                input_vertices[e[j]] = old_ps[j];
        }
//        cout<<"ok2"<<endl;

        ///update
        //
        ns[n_f_ids[0]] = n;
        ns[n_f_ids[1]] = n;
        //
        ts++;
//        for (int f_id: all_f_ids)
//            f_tss[f_id] = ts;
        //
//        std::vector<std::array<int, 2>> new_edges;
//        for (int f_id: all_f_ids) {
//            for (int j = 0; j < 3; j++) {
//                if (input_faces[f_id][j] < input_faces[f_id][(j + 1) % 3])
//                    new_edges.push_back({{input_faces[f_id][j], input_faces[f_id][(j + 1) % 3]}});
//                else
//                    new_edges.push_back({{input_faces[f_id][(j + 1) % 3], input_faces[f_id][j]}});
//            }
//        }
//        vector_unique(new_edges);
//        for(auto& new_e: new_edges){
//            if(new_e == e)
//                continue;
//            std::vector<int> new_n_f_ids;
//            set_intersection(conn_fs[new_e[0]], conn_fs[new_e[1]], new_n_f_ids);
//            if (new_n_f_ids.size() != 2)
//                continue;
//            if (!needs_flattening(ns[new_n_f_ids[0]], ns[new_n_f_ids[1]]))
//                continue;
//            edge_queue.push({{new_e[0], new_e[1], new_n_f_ids[0], new_n_f_ids[1], ts}});
//        }
    }

    cout << "flattening " << ts << " faces" << endl;
}

floatTetWild::Scalar floatTetWild::get_angle_cos(const Vector3& p, const Vector3& p1, const Vector3& p2) {
    Vector3 v1 = p1 - p;
    Vector3 v2 = p2 - p;
    Scalar res = v1.dot(v2) / (v1.norm() * v2.norm());
    if(res > 1)
        return 1;
    if (res < -1)
        return -1;
    return res;
}

bool floatTetWild::is_out_envelope(const std::array<Vector3, 3>& vs, const AABBWrapper& tree, const Parameters& params) {
#ifdef NEW_ENVELOPE
    return tree.is_out_sf_envelope_exact_simplify(vs);
#else
    #ifdef STORE_SAMPLE_POINTS
        std::vector<GEO::vec3> ps;
        sample_triangle(vs, ps, params.dd_simplification);
        return tree.is_out_sf_envelope(ps, params.eps_2_simplification);
    #else
        GEO::index_t prev_facet = GEO::NO_FACET;
        return sample_triangle_and_check_is_out(vs, params.dd_simplification, params.eps_2_simplification, tree, prev_facet);
    #endif
#endif

    // GEO::vec3 init_point(vs[0][0], vs[0][1], vs[0][2]);
    // GEO::vec3 nearest_point;
    // double sq_distg;
    // GEO::index_t prev_facet = tree.nearest_facet(init_point, nearest_point, sq_distg);
    // Scalar sq_dist = sq_distg;
    // if (sq_dist > params.eps_2_simplification)
    //     return true;

    // std::vector<GEO::vec3> ps;
    // sample_triangle(vs, ps, params.dd_simplification);
    // int cnt = 0;
    // const unsigned int ps_size = ps.size();
    // for (unsigned int i = ps_size / 2;; i = (i + 1) % ps_size) {//check from the middle
    //     GEO::vec3 &current_point = ps[i];
    //     sq_distg = current_point.distance2(nearest_point);
    //     tree.nearest_facet_with_hint(current_point, prev_facet, nearest_point, sq_distg);
    //     sq_dist = sq_distg;
    //     if (sq_dist > params.eps_2_simplification)
    //         return true;
    //     cnt++;
    //     if (cnt >= ps_size)
    //         break;
    // }

    // return false;
}

void floatTetWild::check_surface(std::vector<Vector3>& input_vertices, std::vector<Vector3i>& input_faces, const std::vector<bool>& f_is_removed,
                   const AABBWrapper& tree, const Parameters& params) {
    cout<<"checking surface"<<endl;
    bool is_valid = true;
    for (int i = 0; i < input_faces.size(); i++) {
        if(f_is_removed[i])
            continue;
        std::vector<GEO::vec3> ps;
        sample_triangle(
                {{input_vertices[input_faces[i][0]], input_vertices[input_faces[i][1]], input_vertices[input_faces[i][2]]}},
                ps, params.dd_simplification);
        Scalar dist = tree.dist_sf_envelope(ps, params.eps_2);
        if (dist > 0) {
            cout << "is_out_sf_envelope!!" << endl;
            is_valid = false;
//            cout<<input_vertices[input_faces[i][0]].transpose()<<endl;
//            cout<<input_vertices[input_faces[i][1]].transpose()<<endl;
//            cout<<input_vertices[input_faces[i][2]].transpose()<<endl;
            cout<<input_faces[i][0]<<" "<<input_faces[i][1]<<" "<<input_faces[i][2]<<endl;
            cout<<dist<<endl;
//            //pausee();
        }
    }
    // if(!is_valid)
    //     pausee();
}


void floatTetWild::output_component(const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces,
        const std::vector<int>& input_tags){
    return;

    Eigen::MatrixXd V(input_vertices.size(), 3);
    for(int i=0;i<input_vertices.size();i++)
        V.row(i) = input_vertices[i];

    const int C = 2;
    Eigen::MatrixXi F(std::count(input_tags.begin(), input_tags.end(), C), 3);
    int cnt = 0;
    for(int i=0;i<input_tags.size();i++){
        if(input_tags[i] == C)
            F.row(cnt++) = input_faces[i];
    }

    igl::writeOFF("comp.off", V, F);
}
