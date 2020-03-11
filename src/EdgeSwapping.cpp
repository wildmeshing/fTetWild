// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include <floattetwild/EdgeSwapping.h>
#include <floattetwild/LocalOperations.h>

#include <floattetwild/MeshImprovement.h>

namespace floatTetWild {
    bool is_es_check = false;
}

void floatTetWild::edge_swapping(Mesh& mesh) {
    auto &tet_vertices = mesh.tet_vertices;
    auto &tets = mesh.tets;

    int counter = 0;
    int suc_counter3 = 0;
    int suc_counter4 = 0;
    int suc_counter5 = 0;

    mesh.reset_t_empty_start();
    mesh.reset_v_empty_start();

    auto is_swappable = [&](int v1_id, int v2_id, const std::vector<int> &n12_t_ids) {
        if (n12_t_ids.size() < 3 || n12_t_ids.size() > 5)
            return false;
        if (!is_valid_edge(mesh, v1_id, v2_id, n12_t_ids))
            return false;
        if (is_surface_edge(mesh, v1_id, v2_id, n12_t_ids))
            return false;
        if (is_bbox_edge(mesh, v1_id, v2_id, n12_t_ids))
            return false;
        return true;
    };


    std::vector<std::array<int, 2>> edges;
    get_all_edges(mesh, edges);

    std::priority_queue<ElementInQueue, std::vector<ElementInQueue>, cmp_l> es_queue;
    for (auto &e:edges) {
        es_queue.push(ElementInQueue(e, get_edge_length_2(mesh, e[0], e[1])));
    }
    edges.clear();

    while (!es_queue.empty()) {
        std::array<int, 2> v_ids = es_queue.top().v_ids;
        es_queue.pop();

        if(tet_vertices[v_ids[0]].is_freezed && tet_vertices[v_ids[1]].is_freezed)
            continue;

        std::vector<int> n12_t_ids;
        set_intersection(tet_vertices[v_ids[0]].conn_tets, tet_vertices[v_ids[1]].conn_tets, n12_t_ids);
        if (!is_swappable(v_ids[0], v_ids[1], n12_t_ids))
            continue;

        while (!es_queue.empty()) {
            if (es_queue.top().v_ids == v_ids)
                es_queue.pop();
            else
                break;
        }

        bool is_success = false;
        std::vector<std::array<int, 2>> new_edges;
        if (n12_t_ids.size() == 3 && remove_an_edge_32(mesh, v_ids[0], v_ids[1], n12_t_ids, new_edges)) {
            suc_counter3++;
            is_success = true;
//            output_info(mesh);
        }
        if (n12_t_ids.size() == 4 && remove_an_edge_44(mesh, v_ids[0], v_ids[1], n12_t_ids, new_edges)) {
            suc_counter4++;
            is_success = true;
//            output_info(mesh);
        }
        if (n12_t_ids.size() == 5 && remove_an_edge_56(mesh, v_ids[0], v_ids[1], n12_t_ids, new_edges)) {
            suc_counter5++;
            is_success = true;
//            output_info(mesh);
        }

        for (auto &e:new_edges) {
            es_queue.push(ElementInQueue(e, get_edge_length_2(mesh, e[0], e[1])));
        }

        counter++;
    }

    cout << "success3 = " << suc_counter3 << endl;
    cout << "success4 = " << suc_counter4 << endl;
    cout << "success5 = " << suc_counter5 << endl;
    cout << "success = " << (suc_counter3 + suc_counter4 + suc_counter5) << "(" << counter << ")" << endl;
}

bool floatTetWild::remove_an_edge_32(Mesh& mesh, int v1_id, int v2_id, const std::vector<int>& old_t_ids, std::vector<std::array<int, 2>>& new_edges){
    if(old_t_ids.size()!=3)
        return false;

    auto& tet_vertices = mesh.tet_vertices;
    auto& tets = mesh.tets;

    ////construct
    std::array<int, 2> v_ids;
    std::vector<MeshTet> new_tets;
    std::array<int, 2> t_ids;
    int cnt = 0;
    for (int i = 0; i < 4; i++) {
        if (tets[old_t_ids[0]][i] != v1_id && tets[old_t_ids[0]][i] != v2_id) {
            v_ids[cnt++] = tets[old_t_ids[0]][i];
        }
    }
    int i = tets[old_t_ids[1]].find(v_ids[0]);
    if(i>=0){
        new_tets.push_back(tets[old_t_ids[1]]);
        new_tets.push_back(tets[old_t_ids[2]]);
        t_ids = {{old_t_ids[1], old_t_ids[2]}};
    } else {
        new_tets.push_back(tets[old_t_ids[2]]);
        new_tets.push_back(tets[old_t_ids[1]]);
        t_ids = {{old_t_ids[2], old_t_ids[1]}};
    }
    i = new_tets[0].find(v1_id);
    new_tets[0][i] = v_ids[1];
    i = new_tets[1].find(v2_id);
    new_tets[1][i] = v_ids[0];

    ////check
    for(auto& t:new_tets) {
        if (is_inverted(tet_vertices[t[0]], tet_vertices[t[1]], tet_vertices[t[2]], tet_vertices[t[3]]))
            return false;
    }
    std::vector<Scalar> new_qs;
    Scalar old_max_quality = 0;
    for(int t_id: old_t_ids) {
        if (tets[t_id].quality > old_max_quality)
            old_max_quality = tets[t_id].quality;
    }
    for(auto& t:new_tets) {
        Scalar q = get_quality(tet_vertices[t[0]], tet_vertices[t[1]], tet_vertices[t[2]], tet_vertices[t[3]]);
        if (q >= old_max_quality)//or use > ???
            return false;
        new_qs.push_back(q);
    }

    ////real update
    std::vector<std::array<int, 3>> fs;
    std::vector<int> is_sf_fs;
    std::vector<int> sf_tags;
    std::vector<int> is_bx_fs;
    for(int i=0;i<old_t_ids.size();i++) {
        for (int j = 0; j < 4; j++) {
            if (tets[old_t_ids[i]][j] == v1_id || tets[old_t_ids[i]][j] == v2_id) {
                std::array<int, 3> tmp = {{tets[old_t_ids[i]][mod4(j + 1)], tets[old_t_ids[i]][mod4(j + 2)],
                                                  tets[old_t_ids[i]][mod4(j + 3)]}};
                std::sort(tmp.begin(), tmp.end());
                fs.push_back(tmp);
                is_sf_fs.push_back(tets[old_t_ids[i]].is_surface_fs[j]);
                sf_tags.push_back(tets[old_t_ids[i]].surface_tags[j]);
                is_bx_fs.push_back(tets[old_t_ids[i]].is_bbox_fs[j]);
            }
        }
    }

    tets[old_t_ids[0]].is_removed = true;
    tets[t_ids[0]] = new_tets[0];//v2
    tets[t_ids[1]] = new_tets[1];//v1

    for (int i = 0; i < 2; i++)
        tets[t_ids[i]].quality = new_qs[i];

    for(int i=0;i<4;i++) {
        if (tets[t_ids[0]][i] != v2_id) {
            std::array<int, 3> tmp = {{tets[t_ids[0]][mod4(i + 1)], tets[t_ids[0]][mod4(i + 2)],
                                              tets[t_ids[0]][mod4(i + 3)]}};
            std::sort(tmp.begin(), tmp.end());
            auto it = std::find(fs.begin(), fs.end(), tmp);
            tets[t_ids[0]].is_surface_fs[i] = is_sf_fs[it - fs.begin()];
            tets[t_ids[0]].surface_tags[i] = sf_tags[it - fs.begin()];
            tets[t_ids[0]].is_bbox_fs[i] = is_bx_fs[it - fs.begin()];
        } else {
            tets[t_ids[0]].is_surface_fs[i] = NOT_SURFACE;
            tets[t_ids[0]].surface_tags[i] = NO_SURFACE_TAG;
            tets[t_ids[0]].is_bbox_fs[i] = NOT_BBOX;
        }

        if (tets[t_ids[1]][i] != v1_id) {
            std::array<int, 3> tmp = {{tets[t_ids[1]][mod4(i + 1)], tets[t_ids[1]][mod4(i + 2)],
                                              tets[t_ids[1]][mod4(i + 3)]}};
            std::sort(tmp.begin(), tmp.end());
            auto it = std::find(fs.begin(), fs.end(), tmp);
            tets[t_ids[1]].is_surface_fs[i] = is_sf_fs[it - fs.begin()];
            tets[t_ids[1]].surface_tags[i] = sf_tags[it - fs.begin()];
            tets[t_ids[1]].is_bbox_fs[i] = is_bx_fs[it - fs.begin()];
        } else {
            tets[t_ids[1]].is_surface_fs[i] = NOT_SURFACE;
            tets[t_ids[1]].surface_tags[i] = NO_SURFACE_TAG;
            tets[t_ids[1]].is_bbox_fs[i] = NOT_BBOX;
        }
    }

//    tet_vertices[v_ids[0]].conn_tets.erase(old_t_ids[0]);
//    tet_vertices[v_ids[1]].conn_tets.erase(old_t_ids[0]);
//
//    tet_vertices[v_ids[0]].conn_tets.insert(t_ids[1]);
//    tet_vertices[v_ids[1]].conn_tets.insert(t_ids[0]);
//
//    tet_vertices[v1_id].conn_tets.erase(old_t_ids[0]);
//    tet_vertices[v2_id].conn_tets.erase(old_t_ids[0]);
//
//    tet_vertices[v1_id].conn_tets.erase(t_ids[0]);
//    tet_vertices[v2_id].conn_tets.erase(t_ids[1]);

    vector_erase(tet_vertices[v_ids[0]].conn_tets, old_t_ids[0]);
    vector_erase(tet_vertices[v_ids[1]].conn_tets, old_t_ids[0]);

    tet_vertices[v_ids[0]].conn_tets.push_back(t_ids[1]);
    tet_vertices[v_ids[1]].conn_tets.push_back(t_ids[0]);

    vector_erase(tet_vertices[v1_id].conn_tets, old_t_ids[0]);
    vector_erase(tet_vertices[v2_id].conn_tets, old_t_ids[0]);

    vector_erase(tet_vertices[v1_id].conn_tets, t_ids[0]);
    vector_erase(tet_vertices[v2_id].conn_tets, t_ids[1]);

    ////re-push
//    std::unordered_set<int> n12_v_ids;
//    for(int i=0;i<new_tets.size();i++){
//        for(int j=0;j<4;j++){
//            if(new_tets[i][j]!=v1_id && new_tets[i][j]!=v2_id)
//                n12_v_ids.insert(new_tets[i][j]);
//        }
//    }

    new_edges.reserve(new_tets.size()*6);
    for(int i=0;i<new_tets.size();i++) {
        for (int j = 0; j < 3; j++) {
            std::array<int, 2> e = {{new_tets[i][0], new_tets[i][j + 1]}};
            if (e[0] > e[1])
                std::swap(e[0], e[1]);
            new_edges.push_back(e);
            e = {{new_tets[i][j + 1], new_tets[i][mod3(j + 1) + 1]}};
            if (e[0] > e[1])
                std::swap(e[0], e[1]);
            new_edges.push_back(e);
        }
    }
    vector_unique(new_edges);

    return true;
}

bool floatTetWild::remove_an_edge_44(Mesh& mesh, int v1_id, int v2_id, const std::vector<int>& old_t_ids, std::vector<std::array<int, 2>>& new_edges) {
//    bool is_check = false;
//    if(v1_id == 482 && v2_id == 504){
//        is_check = true;
//        pausee("v1_id == 482 && v2_id == 504");
//    }

    const int N = 4;
    if (old_t_ids.size() != N)
        return false;

    auto &tet_vertices = mesh.tet_vertices;
    auto &tets = mesh.tets;

    ////construct
    std::vector<std::array<int, 3>> n12_es;
    n12_es.reserve(old_t_ids.size());
    for (int i = 0; i < old_t_ids.size(); i++) {
        std::array<int, 3> e;
        int cnt = 0;
        for (int j = 0; j < 4; j++)
            if (tets[old_t_ids[i]][j] != v1_id && tets[old_t_ids[i]][j] != v2_id) {
                e[cnt++] = tets[old_t_ids[i]][j];
            }
        e[cnt] = old_t_ids[i];
        n12_es.push_back(e);
    }

    std::vector<int> n12_v_ids;
    std::vector<int> n12_t_ids;
    n12_v_ids.push_back(n12_es[0][0]);
    n12_v_ids.push_back(n12_es[0][1]);
    n12_t_ids.push_back(n12_es[0][2]);
    std::vector<bool> is_visited(N, false);
    is_visited[0] = true;
    for (int i = 0; i < N - 2; i++) {
        for (int j = 0; j < N; j++) {
            if (!is_visited[j]) {
                if (n12_es[j][0] == n12_v_ids.back()) {
                    is_visited[j] = true;
                    n12_v_ids.push_back(n12_es[j][1]);
                } else if (n12_es[j][1] == n12_v_ids.back()) {//else if!!!!!!!!!!
                    is_visited[j] = true;
                    n12_v_ids.push_back(n12_es[j][0]);
                }
                if (is_visited[j]) {
                    n12_t_ids.push_back(n12_es[j][2]);
                    break;
                }
            }
        }
    }
    n12_t_ids.push_back(n12_es[std::find(is_visited.begin(), is_visited.end(), false) - is_visited.begin()][2]);

    ////check
    bool is_valid = false;
//    std::vector<MeshTet> new_tets;
    std::vector<Vector4i> new_tets;
    new_tets.reserve(4);
    std::vector<int> tags;
    std::array<int, 2> v_ids;
    std::vector<Scalar> new_qs;
    Scalar old_max_quality = 0;
    Scalar new_max_quality = 0;
    for (int t_id: old_t_ids) {
        if (tets[t_id].quality > old_max_quality)
            old_max_quality = tets[t_id].quality;
    }
    for (int i = 0; i < 2; i++) {
//        std::vector<MeshTet> tmp_new_tets;
        std::vector<Vector4i> tmp_new_tets;
        std::vector<int> tmp_tags;
        std::array<int, 2> tmp_v_ids;
        tmp_v_ids = {{n12_v_ids[0 + i], n12_v_ids[2 + i]}};
        bool is_break = false;
        for (int j = 0; j < old_t_ids.size(); j++) {
            auto t = tets[old_t_ids[j]];
            int ii = t.find(tmp_v_ids[0]);
            if (ii >= 0) {
                int jt = t.find(v2_id);
                t[jt] = tmp_v_ids[1];
                tmp_tags.push_back(1);
            } else {
                int jt = t.find(v1_id);
                t[jt] = tmp_v_ids[0];
                tmp_tags.push_back(0);
            }
            if (is_inverted(tet_vertices[t[0]], tet_vertices[t[1]], tet_vertices[t[2]], tet_vertices[t[3]])) {
                is_break = true;
                break;
            }
//            tmp_new_tets.push_back(t);
            tmp_new_tets.push_back(t.indices);
        }
        if (is_break)
            continue;

        std::vector<Scalar> tmp_new_qs;
        for (auto &t: tmp_new_tets) {
            Scalar q = get_quality(tet_vertices[t[0]], tet_vertices[t[1]], tet_vertices[t[2]], tet_vertices[t[3]]);
            if (q >= old_max_quality) {
                is_break = true;
                break;
//                return false;
            }
            if (q > new_max_quality)
                new_max_quality = q;
            tmp_new_qs.push_back(q);
        }
        if (is_break)
            continue;

        is_valid = true;
        old_max_quality = new_max_quality;
        new_tets = tmp_new_tets;
        tags = tmp_tags;
        new_qs = tmp_new_qs;
        v_ids = tmp_v_ids;
    }
    if (!is_valid)
        return false;

    ////real update
    std::vector<std::array<int, 3>> fs;
    std::vector<int> is_sf_fs;
    std::vector<int> sf_tags;
    std::vector<int> is_bx_fs;
    for (int i = 0; i < old_t_ids.size(); i++) {
        for (int j = 0; j < 4; j++) {
            if (tets[old_t_ids[i]][j] == v1_id || tets[old_t_ids[i]][j] == v2_id) {
                std::array<int, 3> tmp = {{tets[old_t_ids[i]][mod4(j + 1)], tets[old_t_ids[i]][mod4(j + 2)],
                                                  tets[old_t_ids[i]][mod4(j + 3)]}};
                std::sort(tmp.begin(), tmp.end());
                fs.push_back(tmp);
                is_sf_fs.push_back(tets[old_t_ids[i]].is_surface_fs[j]);
                sf_tags.push_back(tets[old_t_ids[i]].surface_tags[j]);
                is_bx_fs.push_back(tets[old_t_ids[i]].is_bbox_fs[j]);
            }
        }
    }

    for (int j = 0; j < new_tets.size(); j++) {
        if (tags[j] == 0) {
//            tet_vertices[v1_id].conn_tets.erase(old_t_ids[j]);
//            tet_vertices[v_ids[0]].conn_tets.insert(old_t_ids[j]);
            vector_erase(tet_vertices[v1_id].conn_tets, old_t_ids[j]);
            tet_vertices[v_ids[0]].conn_tets.push_back(old_t_ids[j]);
        } else {
//            tet_vertices[v2_id].conn_tets.erase(old_t_ids[j]);
//            tet_vertices[v_ids[1]].conn_tets.insert(old_t_ids[j]);
            vector_erase(tet_vertices[v2_id].conn_tets, old_t_ids[j]);
            tet_vertices[v_ids[1]].conn_tets.push_back(old_t_ids[j]);
        }
//        tets[old_t_ids[j]] = new_tets[j];
        tets[old_t_ids[j]].indices = new_tets[j];
        tets[old_t_ids[j]].quality = new_qs[j];
    }

    for (int i = 0; i < old_t_ids.size(); i++) {//old_t_ids contains new tets
        for (int j = 0; j < 4; j++) {
            tets[old_t_ids[i]].is_surface_fs[j] = NOT_SURFACE;
            tets[old_t_ids[i]].surface_tags[j] = NO_SURFACE_TAG;
            tets[old_t_ids[i]].is_bbox_fs[j] = NOT_BBOX;
            if (tets[old_t_ids[i]][j] == v_ids[0] || tets[old_t_ids[i]][j] == v_ids[1]) {
                std::array<int, 3> tmp = {{tets[old_t_ids[i]][mod4(j + 1)], tets[old_t_ids[i]][mod4(j + 2)],
                                                  tets[old_t_ids[i]][mod4(j + 3)]}};
                std::sort(tmp.begin(), tmp.end());
                auto it = std::find(fs.begin(), fs.end(), tmp);
                tets[old_t_ids[i]].is_surface_fs[j] = is_sf_fs[it - fs.begin()];
                tets[old_t_ids[i]].surface_tags[j] = sf_tags[it - fs.begin()];
                tets[old_t_ids[i]].is_bbox_fs[j] = is_bx_fs[it - fs.begin()];
            }
        }
    }

    ////re-push
    new_edges.reserve(new_tets.size() * 6);
    for (int i = 0; i < new_tets.size(); i++) {
        for (int j = 0; j < 3; j++) {
            std::array<int, 2> e = {{new_tets[i][0], new_tets[i][j + 1]}};
            if (e[0] > e[1])
                std::swap(e[0], e[1]);
            new_edges.push_back(e);
            e = {{new_tets[i][j + 1], new_tets[i][mod3(j + 1) + 1]}};
            if (e[0] > e[1])
                std::swap(e[0], e[1]);
            new_edges.push_back(e);
        }
    }
    vector_unique(new_edges);

    return true;
}

#include <unordered_map>
bool floatTetWild::remove_an_edge_56(Mesh& mesh, int v1_id, int v2_id, const std::vector<int>& old_t_ids, std::vector<std::array<int, 2>>& new_edges) {
    if (old_t_ids.size() != 5)
        return false;

    auto &tet_vertices = mesh.tet_vertices;
    auto &tets = mesh.tets;

    ////construct
    std::vector<std::array<int, 3>> n12_es;
    n12_es.reserve(old_t_ids.size());
    for (int i = 0; i < old_t_ids.size(); i++) {
        std::array<int, 3> e;
        int cnt = 0;
        for (int j = 0; j < 4; j++)
            if (tets[old_t_ids[i]][j] != v1_id && tets[old_t_ids[i]][j] != v2_id) {
                e[cnt++] = tets[old_t_ids[i]][j];
            }
        e[cnt] = old_t_ids[i];
        n12_es.push_back(e);
    }

    std::vector<int> n12_v_ids;
    std::vector<int> n12_t_ids;
    n12_v_ids.push_back(n12_es[0][0]);
    n12_v_ids.push_back(n12_es[0][1]);
    n12_t_ids.push_back(n12_es[0][2]);
    std::vector<bool> is_visited(5, false);
    is_visited[0] = true;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 5; j++) {
            if (!is_visited[j]) {
                if (n12_es[j][0] == n12_v_ids.back()) {
                    is_visited[j] = true;
                    n12_v_ids.push_back(n12_es[j][1]);
                } else if (n12_es[j][1] == n12_v_ids.back()) {//else if!!!!!!!!!!
                    is_visited[j] = true;
                    n12_v_ids.push_back(n12_es[j][0]);
                }
                if (is_visited[j]) {
                    n12_t_ids.push_back(n12_es[j][2]);
                    break;
                }
            }
        }
    }
    n12_t_ids.push_back(n12_es[std::find(is_visited.begin(), is_visited.end(), false) - is_visited.begin()][2]);

    ////check
    Scalar old_max_quality = 0;
    Scalar new_max_quality = 0;
    for (int t_id: old_t_ids) {
        if (tets[t_id].quality > old_max_quality)
            old_max_quality = tets[t_id].quality;
    }

    std::unordered_map<int, std::array<Scalar, 2>> tet_qs;
//    std::unordered_map<int, std::array<MeshTet, 2>> new_tets;
    std::unordered_map<int, std::array<Vector4i, 2>> new_tets;
    std::vector<bool> is_v_valid(5, true);
    for (int i = 0; i < n12_v_ids.size(); i++) {
        if (!is_v_valid[(i + 1) % 5] && !is_v_valid[(i - 1 + 5) % 5])
            continue;

//        std::vector<MeshTet> new_ts;
        std::vector<Vector4i> new_ts;
        new_ts.reserve(6);
        auto t = tets[n12_t_ids[i]];
        int it = t.find(v1_id);
        t[it] = n12_v_ids[(i - 1 + 5) % 5];
        if (is_inverted(tet_vertices[t[0]], tet_vertices[t[1]], tet_vertices[t[2]], tet_vertices[t[3]])) {
            is_v_valid[(i + 1) % 5] = false;
            is_v_valid[(i - 1 + 5) % 5] = false;
            continue;
        }
//        new_ts.push_back(t);
        new_ts.push_back(t.indices);

        t = tets[n12_t_ids[i]];
        it = t.find(v2_id);
        t[it] = n12_v_ids[(i - 1 + 5) % 5];
        if (is_inverted(tet_vertices[t[0]], tet_vertices[t[1]], tet_vertices[t[2]], tet_vertices[t[3]])) {
            is_v_valid[(i + 1) % 5] = false;
            is_v_valid[(i - 1 + 5) % 5] = false;
            continue;
        }
//        new_ts.push_back(t);
        new_ts.push_back(t.indices);
//        new_tets[i] = std::array<MeshTet, 2>({{new_ts[0], new_ts[1]}});
        new_tets[i] = std::array<Vector4i, 2>({{new_ts[0], new_ts[1]}});

        std::vector<Scalar> qs;
        for (auto &t: new_ts) {
            qs.emplace_back(
                    get_quality(tet_vertices[t[0]], tet_vertices[t[1]], tet_vertices[t[2]], tet_vertices[t[3]]));
        }
        tet_qs[i] = std::array<Scalar, 2>({{qs[0], qs[1]}});
    }
    if (std::count(is_v_valid.begin(), is_v_valid.end(), true) == 0)
        return false;

    int selected_id = -1;
    for (int i = 0; i < is_v_valid.size(); i++) {
        if (!is_v_valid[i])
            continue;

//        std::vector<MeshTet> new_ts;
        std::vector<Vector4i> new_ts;
        new_ts.reserve(6);
        auto t = tets[n12_t_ids[(i + 2) % 5]];
        int it = t.find(v1_id);
        t[it] = n12_v_ids[i];
        if (is_inverted(tet_vertices[t[0]], tet_vertices[t[1]], tet_vertices[t[2]], tet_vertices[t[3]]))
            continue;
//        new_ts.push_back(t);
        new_ts.push_back(t.indices);
        t = tets[n12_t_ids[(i + 2) % 5]];
        it = t.find(v2_id);
        t[it] = n12_v_ids[i];
        if (is_inverted(tet_vertices[t[0]], tet_vertices[t[1]], tet_vertices[t[2]], tet_vertices[t[3]]))
            continue;
//        new_ts.push_back(t);
        new_ts.push_back(t.indices);

        std::vector<Scalar> qs;
        for (auto &t:new_ts) {
            qs.push_back(get_quality(tet_vertices[t[0]], tet_vertices[t[1]], tet_vertices[t[2]], tet_vertices[t[3]]));
        }
        for (int j = 0; j < 2; j++) {
            qs.push_back(tet_qs[(i + 1) % 5][j]);
            qs.push_back(tet_qs[(i - 1 + 5) % 5][j]);
        }
        if (qs.size() != 6) {
            assert("qs.size() != 6");
        }
        for (auto &q:qs) {
            if (q > new_max_quality)
                new_max_quality = q;
        }
        if (new_max_quality >= old_max_quality)
            continue;

        old_max_quality = new_max_quality;
        selected_id = i;
        tet_qs[i + 5] = std::array<Scalar, 2>({{qs[0], qs[1]}});
//        new_tets[i + 5] = std::array<MeshTet, 2>({{new_ts[0], new_ts[1]}});
        new_tets[i + 5] = std::array<Vector4i, 2>({{new_ts[0], new_ts[1]}});
    }
    if (selected_id < 0)
        return false;


    ////real update
    //update on surface -- 1
    std::vector<std::array<int, 3>> fs;
    std::vector<int> is_sf_fs;
    std::vector<int> sf_tags;
    std::vector<int> is_bx_fs;
    for (int i = 0; i < old_t_ids.size(); i++) {
        for (int j = 0; j < 4; j++) {
            if (tets[old_t_ids[i]][j] == v1_id || tets[old_t_ids[i]][j] == v2_id) {
                std::array<int, 3> tmp = {{tets[old_t_ids[i]][(j + 1) % 4], tets[old_t_ids[i]][(j + 2) % 4],
                                                  tets[old_t_ids[i]][(j + 3) % 4]}};
                std::sort(tmp.begin(), tmp.end());
                fs.push_back(tmp);
                is_sf_fs.push_back(tets[old_t_ids[i]].is_surface_fs[j]);
                sf_tags.push_back(tets[old_t_ids[i]].surface_tags[j]);
                is_bx_fs.push_back(tets[old_t_ids[i]].is_bbox_fs[j]);
            }
        }
    }

    std::vector<int> new_t_ids = old_t_ids;
    get_new_tet_slots(mesh, 1, new_t_ids);
    tets[new_t_ids.back()].reset();

    for (int i = 0; i < 2; i++) {
        tets[new_t_ids[i]] = new_tets[(selected_id + 1) % 5][i];
        tets[new_t_ids[i + 2]] = new_tets[(selected_id - 1 + 5) % 5][i];
        tets[new_t_ids[i + 4]] = new_tets[selected_id + 5][i];
//        tets[new_t_ids[i]].indices = new_tets[(selected_id + 1) % 5][i];
//        tets[new_t_ids[i + 2]].indices = new_tets[(selected_id - 1 + 5) % 5][i];
//        tets[new_t_ids[i + 4]].indices = new_tets[selected_id + 5][i];

        tets[new_t_ids[i]].quality = tet_qs[(selected_id + 1) % 5][i];
        tets[new_t_ids[i + 2]].quality = tet_qs[(selected_id - 1 + 5) % 5][i];
        tets[new_t_ids[i + 4]].quality = tet_qs[selected_id + 5][i];
    }

    //update on_surface -- 2
    for (int i = 0; i < new_t_ids.size(); i++) {
        for (int j = 0; j < 4; j++) {
            tets[new_t_ids[i]].is_surface_fs[j] = NOT_SURFACE;
            tets[new_t_ids[i]].surface_tags[j] = NO_SURFACE_TAG;
            tets[new_t_ids[i]].is_bbox_fs[j] = NOT_BBOX;
            if (tets[new_t_ids[i]][j] != v1_id && tets[new_t_ids[i]][j] != v2_id
                && tets[new_t_ids[i]][j] != n12_v_ids[(selected_id + 1) % 5]
                && tets[new_t_ids[i]][j] != n12_v_ids[(selected_id - 1 + 5) % 5]) {
                std::array<int, 3> tmp = {{tets[new_t_ids[i]][(j + 1) % 4], tets[new_t_ids[i]][(j + 2) % 4],
                                                  tets[new_t_ids[i]][(j + 3) % 4]}};
                std::sort(tmp.begin(), tmp.end());
                auto it = std::find(fs.begin(), fs.end(), tmp);
                if (it != fs.end()) {
                    tets[new_t_ids[i]].is_surface_fs[j] = is_sf_fs[it - fs.begin()];
                    tets[new_t_ids[i]].surface_tags[j] = sf_tags[it - fs.begin()];
                    tets[new_t_ids[i]].is_bbox_fs[j] = is_bx_fs[it - fs.begin()];
                }
            }
        }
    }

    //update conn_tets
    for (int i = 0; i < n12_v_ids.size(); i++) {
//        tet_vertices[n12_v_ids[i]].conn_tets.erase(n12_t_ids[i]);
//        tet_vertices[n12_v_ids[i]].conn_tets.erase(n12_t_ids[(i - 1 + 5) % 5]);
        vector_erase(tet_vertices[n12_v_ids[i]].conn_tets, n12_t_ids[i]);
        vector_erase(tet_vertices[n12_v_ids[i]].conn_tets, n12_t_ids[(i - 1 + 5) % 5]);
    }
    for (int i = 0; i < n12_t_ids.size(); i++) {
        vector_erase(tet_vertices[v1_id].conn_tets, n12_t_ids[i]);
        vector_erase(tet_vertices[v2_id].conn_tets, n12_t_ids[i]);
    }

    //add
    for (int i = 0; i < new_t_ids.size(); i++) {
        for (int j = 0; j < 4; j++)
//            tet_vertices[tets[new_t_ids[i]][j]].conn_tets.insert(new_t_ids[i]);
            tet_vertices[tets[new_t_ids[i]][j]].conn_tets.push_back(new_t_ids[i]);
    }

    ////re-push
    new_edges.reserve(new_t_ids.size() * 6);
    for (int i = 0; i < new_t_ids.size(); i++) {
        for (int j = 0; j < 3; j++) {
            std::array<int, 2> e = {{tets[new_t_ids[i]][0], tets[new_t_ids[i]][j + 1]}};
            if (e[0] > e[1])
                std::swap(e[0], e[1]);
            new_edges.push_back(e);
            e = {{tets[new_t_ids[i]][j + 1], tets[new_t_ids[i]][mod3(j + 1) + 1]}};
            if (e[0] > e[1])
                std::swap(e[0], e[1]);
            new_edges.push_back(e);
        }
    }
    vector_unique(new_edges);

    return true;
}