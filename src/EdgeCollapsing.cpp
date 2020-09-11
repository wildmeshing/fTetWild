// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include <floattetwild/EdgeCollapsing.h>
#include <floattetwild/LocalOperations.h>

#include <floattetwild/MeshImprovement.h> //todo: tmp


#ifdef FLOAT_TETWILD_USE_TBB
#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include <tbb/atomic.h>

#endif

#define EC_FAIL_INVERSION -1
#define EC_FAIL_QUALITY -2
#define EC_FAIL_ENVELOPE0 -3
#define EC_FAIL_ENVELOPE1 -4
#define EC_FAIL_ENVELOPE2 -5
#define EC_FAIL_ENVELOPE3 -6
#define EC_SUCCESS 1
#define EC_SUCCESS_ENVELOPE 2

#define EC_POSTPROCESS true


namespace floatTetWild {
namespace {
void edge_collapsing_aux(Mesh& mesh, const AABBWrapper& tree, std::vector<std::array<int, 2>> &edges) {
    auto &tets = mesh.tets;
    auto &tet_vertices = mesh.tet_vertices;

    int counter = 0;
    int suc_counter = 0;
    int suc_counter_env = 0;

    ////init
    std::priority_queue<ElementInQueue, std::vector<ElementInQueue>, cmp_s> ec_queue;
    for (auto &e:edges) {
        Scalar l_2 = get_edge_length_2(mesh, e[0], e[1]);
        if (is_collapsable_length(mesh, e[0], e[1], l_2) && is_collapsable_boundary(mesh, e[0], e[1], tree)) {
            ec_queue.push(ElementInQueue(e, l_2));
            ec_queue.push(ElementInQueue({{e[1], e[0]}}, l_2));
        }
    }
    edges.clear();

    ////collapse
    int ts = 0;
    std::vector<std::array<int, 2>> inf_es;
    std::vector<int> inf_e_tss;
    std::vector<int> tet_tss;
    tet_tss.assign(tets.size(), 0);

    do {
        counter = 0;
        suc_counter = 0;
        suc_counter_env = 0;
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

            if(is_edge_freezed(mesh, v_ids[0], v_ids[1]))
                continue;

            if (!is_valid_edge(mesh, v_ids[0], v_ids[1]))
                continue;

            if(! is_collapsable_boundary(mesh, v_ids[0], v_ids[1], tree))
                continue;

            Scalar weight = get_edge_length_2(mesh, v_ids[0], v_ids[1]);
            if (weight != old_weight || !is_collapsable_length(mesh, v_ids[0], v_ids[1], weight))
                continue;

            if (!is_collapsable_bbox(mesh, v_ids[0], v_ids[1]))
                continue;

            std::vector<std::array<int, 2>> new_edges;
            int result = collapse_an_edge(mesh, v_ids[0], v_ids[1], tree, new_edges, ts, tet_tss);
            if (result == EC_SUCCESS || result == EC_SUCCESS_ENVELOPE) {
                if (result == EC_SUCCESS_ENVELOPE)
                    suc_counter_env++;
                suc_counter++;

                for (auto &e:new_edges) {
                    if(is_edge_freezed(mesh, e[0], e[1]))
                        continue;

                    Scalar l_2 = get_edge_length_2(mesh, e[0], e[1]);
                    if (is_collapsable_length(mesh, e[0], e[1], l_2)) {
                        ec_queue.push(ElementInQueue(e, l_2));
                        ec_queue.push(ElementInQueue({{e[1], e[0]}}, l_2));
                    }
                }
            }
#if EC_POSTPROCESS
            else {
//                //fortest
//                int v1_id = v_ids[0];
//                int v2_id = v_ids[1];
////                std::vector<int> n12_t_ids;
////                set_intersection(tet_vertices[v1_id].conn_tets, tet_vertices[v2_id].conn_tets, n12_t_ids);
//
//                Scalar old_max_quality = 0;
//                std::vector<Scalar> new_qs;
//                new_qs.reserve(tet_vertices[v1_id].conn_tets.size());
//                for (int t_id:tet_vertices[v1_id].conn_tets) {
//                    if (tets[t_id].quality > old_max_quality)
//                        old_max_quality = tets[t_id].quality;
//                }
//                if(old_max_quality>5e9){
//                    cout<<"resilt = "<<result<<endl;
//                    cout<<"old_max_quality = "<<old_max_quality<<endl;
//                    cout<<"e = "<<v1_id<<" "<<v2_id<<endl;
//                }
//                Eigen::MatrixXd V;
//                Eigen::MatrixXi T;
//
//                //fortest

//                if(weight<SCALAR_ZERO_2){
//                    cout<<"len = "<<weight<<" but failed "<<result<<" "<<mesh.is_input_all_inserted<<endl;
//                    //pausee();
//                }
//                if(counter < 25){
//                    cout<<"len = "<<weight<<" but failed "<<result<<" "<<mesh.is_input_all_inserted<<endl;
//                    //pausee();
//                }
                inf_es.push_back(v_ids);
                inf_e_tss.push_back(ts);
            }
#endif

            counter++;
        }

        cout << "success(env) = " << suc_counter_env << endl;
        cout << "success = " << suc_counter << "(" << counter << ")" << endl;

#if EC_POSTPROCESS
        if (suc_counter == 0)
#endif
            break;

#if EC_POSTPROCESS
        ////postprocess
        std::vector<std::array<int, 2>> tmp_inf_es;
        const unsigned int inf_es_size = inf_es.size();
        tmp_inf_es.reserve(inf_es_size / 4.0 + 1);
        for (unsigned int i = 0; i < inf_es_size; i++) {
            if(is_edge_freezed(mesh, inf_es[i][0], inf_es[i][1]))
                continue;
            if (!is_valid_edge(mesh, inf_es[i][0], inf_es[i][1]))
                continue;

            Scalar weight = get_edge_length_2(mesh, inf_es[i][0], inf_es[i][1]);
            if (!is_collapsable_length(mesh, inf_es[i][0], inf_es[i][1], weight))
                continue;

            if (!is_collapsable_bbox(mesh, inf_es[i][0], inf_es[i][1]))
                continue;

            bool is_recal = false;
            for (int t_id: tet_vertices[inf_es[i][0]].conn_tets) {
                if (tet_tss[t_id] > inf_e_tss[i]) {
                    is_recal = true;
                    break;
                }
            }
            if (is_recal)
                ec_queue.push(ElementInQueue(inf_es[i], weight));
            else
                tmp_inf_es.push_back(inf_es[i]);
        }
        std::sort(tmp_inf_es.begin(), tmp_inf_es.end());
        tmp_inf_es.erase(std::unique(tmp_inf_es.begin(), tmp_inf_es.end()), tmp_inf_es.end());//it's better
        inf_es = tmp_inf_es;

        ts++;
        inf_e_tss = std::vector<int>(inf_es.size(), ts);
#endif
    } while (suc_counter > 0);
}

}
}

void floatTetWild::edge_collapsing(Mesh& mesh, const AABBWrapper& tree) {
    auto &tet_vertices = mesh.tet_vertices;
    auto &tets = mesh.tets;

    std::vector<std::array<int, 2>> edges;

// #ifdef FLOAT_TETWILD_USE_TBB_bug //TODO: remove bug and fix
//     std::vector<std::vector<int>> partition;
//     int num_partition = tbb::task_scheduler_init::default_num_threads();
//     partition.clear();
//     cout << "num_partition = " << num_partition << endl;
//     mesh.partition(num_partition, partition);

//     for (int p_id = 0; p_id < partition.size(); p_id++) {
//         for (int t_id: partition[p_id])
//             tets[t_id].scalar = p_id;
//     }
//     for (auto &v:tet_vertices) {
//         std::unordered_set<int> p_ids;
//         v.is_freezed = false;
//         for (int t_id: v.conn_tets) {
//             p_ids.insert(tets[t_id].scalar);
//             if (p_ids.size() > 1) {
//                 v.is_freezed = true;
//                 break;
//             }
//         }
//     }

//     const bool skip_freezed = true;
//     tbb::parallel_for(size_t(0), size_t(partition.size()), [&](size_t i) {
//         std::vector<std::array<int, 2>> edges;
//         get_all_edges(mesh, partition[i], edges, skip_freezed);
//         edge_collapsing_aux(mesh, tree, edges);
//     });

// //    for (auto &t: mesh.tets) {
// //        for (int j = 0; j < 3; j++) {
// //            if (mesh.tet_vertices[t[0]].is_freezed || mesh.tet_vertices[t[j + 1]].is_freezed) {
// //                std::array<int, 2> e = {{t[0], t[j + 1]}};
// //                if (e[0] > e[1])
// //                    std::swap(e[0], e[1]);
// //                edges.push_back(e);
// //            }
// //            if (mesh.tet_vertices[t[j + 1]].is_freezed || mesh.tet_vertices[mod3(j + 1) + 1].is_freezed) {
// //                std::array<int, 2> e = {{t[j + 1], t[mod3(j + 1) + 1]}};
// //                if (e[0] > e[1])
// //                    std::swap(e[0], e[1]);
// //                edges.push_back(e);
// //            }
// //        }
// //    }
// //    vector_unique(edges);
//     get_all_edges(mesh, edges);
//     edge_collapsing_aux(mesh, tree, edges);

//     for (auto &v:tet_vertices)
//         v.is_freezed = false;
//     for (auto &t:tets)
//         t.scalar = 0;
// #else
    get_all_edges(mesh, edges);
    edge_collapsing_aux(mesh, tree, edges);
// #endif
}

int floatTetWild::collapse_an_edge(Mesh& mesh, int v1_id, int v2_id, const AABBWrapper& tree,
        std::vector<std::array<int, 2>>& new_edges, int ts, std::vector<int>& tet_tss,
        bool is_check_quality, bool is_update_tss) {
    auto &tet_vertices = mesh.tet_vertices;
    auto &tets = mesh.tets;

    ////check vertices
    //check isolate surface points
    if (tet_vertices[v1_id].is_on_surface && is_isolate_surface_point(mesh, v1_id)) {
        tet_vertices[v1_id].is_on_surface = false;
        tet_vertices[v1_id].is_on_boundary = false;
    }
    //check boundary/surface
    if (tet_vertices[v1_id].is_on_boundary && is_point_out_boundary_envelope(mesh, tet_vertices[v2_id].pos, tree))//todo: you should check/unmark is_on_boundary around here
        return EC_FAIL_ENVELOPE0;
    if (tet_vertices[v1_id].is_on_surface && is_point_out_envelope(mesh, tet_vertices[v2_id].pos, tree))
        return EC_FAIL_ENVELOPE1;


    ////check tets
    std::vector<int> n12_t_ids;
    set_intersection(tet_vertices[v1_id].conn_tets, tet_vertices[v2_id].conn_tets, n12_t_ids);
    if(n12_t_ids.empty())
        return EC_FAIL_INVERSION;
//    std::unordered_set<int> n1_t_ids = tet_vertices[v1_id].conn_tets;//v1.conn_tets - n12_t_ids
//    for (int t_id:n12_t_ids)
//        n1_t_ids.erase(t_id);
    std::vector<int> n1_t_ids;//v1.conn_tets - n12_t_ids
    std::sort(tet_vertices[v1_id].conn_tets.begin(), tet_vertices[v1_id].conn_tets.end());
    std::sort(n12_t_ids.begin(), n12_t_ids.end());
    std::set_difference(tet_vertices[v1_id].conn_tets.begin(), tet_vertices[v1_id].conn_tets.end(),
            n12_t_ids.begin(), n12_t_ids.end(), std::back_inserter(n1_t_ids));

    //inversion
    std::vector<int> js_n1_t_ids;
    for (int t_id:n1_t_ids) {
        int j = tets[t_id].find(v1_id);
        js_n1_t_ids.push_back(j);
        assert(j < 4);
        if (is_inverted(mesh, t_id, j, tet_vertices[v2_id].pos))
            return EC_FAIL_INVERSION;
    }

    //quality
    Scalar old_max_quality = 0;
    if(mesh.is_coarsening){
        old_max_quality = mesh.params.stop_energy;
    } else {
        if (is_check_quality) {
            for (int t_id:tet_vertices[v1_id].conn_tets) {
                if (tets[t_id].quality > old_max_quality)
                    old_max_quality = tets[t_id].quality;
            }
        }
    }
    std::vector<Scalar> new_qs;
    new_qs.reserve(tet_vertices[v1_id].conn_tets.size());
    int ii = 0;
    for (int t_id:n1_t_ids) {
        int j = js_n1_t_ids[ii++];
        Scalar new_q = get_quality(tet_vertices[v2_id], tet_vertices[tets[t_id][mod4(j + 1)]],
                                   tet_vertices[tets[t_id][mod4(j + 2)]], tet_vertices[tets[t_id][mod4(j + 3)]]);
        if (is_check_quality && new_q > old_max_quality)
            return EC_FAIL_QUALITY;
        new_qs.push_back(new_q);
    }

    //envelope
    Scalar l = get_edge_length_2(mesh, v1_id, v2_id);
    if (l > 0) {
        if (tet_vertices[v1_id].is_on_boundary) {
            if (is_out_boundary_envelope(mesh, v1_id, tet_vertices[v2_id].pos, tree))
                return EC_FAIL_ENVELOPE2;
        }
        if (tet_vertices[v1_id].is_on_surface) {
            if (is_out_envelope(mesh, v1_id, tet_vertices[v2_id].pos, tree))
                return EC_FAIL_ENVELOPE3;
        }
    }


    ////real update
    //vertex
    tet_vertices[v1_id].is_removed = true;
    tet_vertices[v2_id].is_on_bbox = tet_vertices[v1_id].is_on_bbox || tet_vertices[v2_id].is_on_bbox;
    tet_vertices[v2_id].is_on_surface = tet_vertices[v1_id].is_on_surface || tet_vertices[v2_id].is_on_surface;
    tet_vertices[v2_id].is_on_boundary = tet_vertices[v1_id].is_on_boundary || tet_vertices[v2_id].is_on_boundary;
    if(tet_vertices[v1_id].on_boundary_e_id >= 0)
        tet_vertices[v2_id].on_boundary_e_id = tet_vertices[v1_id].on_boundary_e_id;

    //tets
    //update quality
    int i=0;
    for (int t_id:n1_t_ids) {
        tets[t_id].quality = new_qs[i++];
    }

    //n_v_id for repush
//    std::vector<int> n12_v_ids;
//    std::vector<int> n1_v_ids;
//    n12_v_ids.reserve(n12_t_ids.size() * 4);
//    for (int t_id:n12_t_ids) {
//        for (int j = 0; j < 4; j++)
//            n12_v_ids.push_back(tets[t_id][j]);
//    }
//    n1_v_ids.reserve(n1_t_ids.size() * 4);
//    for (int t_id:n1_t_ids) {
//        for (int j = 0; j < 4; j++)
//            n1_v_ids.push_back(tets[t_id][j]);
//    }
//    vector_unique(n12_v_ids);
//    vector_unique(n1_v_ids);
//    std::vector<int> n_v_ids;
//    std::set_difference(n1_v_ids.begin(), n1_v_ids.end(), n12_v_ids.begin(), n12_v_ids.end(),
//                        std::inserter(n_v_ids, n_v_ids.begin()));
    std::vector<int> n1_v_ids;
    n1_v_ids.reserve(n1_t_ids.size() * 4);
    for (int t_id:n1_t_ids) {
        for (int j = 0; j < 4; j++)
            n1_v_ids.push_back(tets[t_id][j]);
    }
    vector_unique(n1_v_ids);

    //update tags
//    cout<<"n12_t_ids = ";
//    vector_print(n12_t_ids, " ");

    for (int t_id:n12_t_ids) {
//        cout<<"t_id = "<<t_id<<endl;
        int sf_facing_v1 = NOT_SURFACE;
        int sf_facing_v2 = NOT_SURFACE;
        int tag_facing_v1 = NO_SURFACE_TAG;
        int tag_facing_v2 = NO_SURFACE_TAG;
        int bbox_facing_v1 = NOT_BBOX;
        int bbox_facing_v2 = NOT_BBOX;

        std::array<int, 2> j12;
        for (int j = 0; j < 4; j++) {
            if (tets[t_id][j] == v1_id) {
                sf_facing_v1 = tets[t_id].is_surface_fs[j];
                tag_facing_v1 = tets[t_id].surface_tags[j];
                bbox_facing_v1 = tets[t_id].is_bbox_fs[j];
                j12[0] = j;
            } else if (tets[t_id][j] == v2_id) {
                sf_facing_v2 = tets[t_id].is_surface_fs[j];
                tag_facing_v2 = tets[t_id].surface_tags[j];
                bbox_facing_v2 = tets[t_id].is_bbox_fs[j];
                j12[1] = j;
            }
        }

        std::array<int, 2> sf_connecting_v12 = {{NOT_SURFACE, NOT_SURFACE}};
        std::array<int, 2> tag_connecting_v12 = {{NO_SURFACE_TAG, NO_SURFACE_TAG}};
        if (sf_facing_v2 != NOT_SURFACE && sf_facing_v1 != NOT_SURFACE) {
//            sf_connecting_v12[0] = -sf_facing_v2 + sf_facing_v1;
            int new_tag = -sf_facing_v2 + sf_facing_v1;
            if(new_tag == 0)
//                sf_connecting_v12[0] = NOT_SURFACE;
                sf_connecting_v12[0] = 0;
            else if (new_tag > 0)
                sf_connecting_v12[0] = 1;
            else
                sf_connecting_v12[0] = -1;
            tag_connecting_v12[0] = tag_facing_v1;
        } else if (sf_facing_v2 != NOT_SURFACE) {
            sf_connecting_v12[0] = -sf_facing_v2;
            tag_connecting_v12[0] = tag_facing_v2;
        } else if (sf_facing_v1 != NOT_SURFACE) {
            sf_connecting_v12[0] = sf_facing_v1;
            tag_connecting_v12[0] = tag_facing_v1;
        }
        if(sf_connecting_v12[0] != NOT_SURFACE) {
            sf_connecting_v12[1] = -sf_connecting_v12[0];
            tag_connecting_v12[1] = tag_connecting_v12[0];
        }

        std::array<int, 2> bbox_connecting_v12;
        bbox_connecting_v12[0] = NOT_BBOX;
        if (bbox_facing_v2 != NOT_BBOX)
            bbox_connecting_v12[0] = bbox_facing_v2;
        else if (sf_facing_v1 != NOT_BBOX)
            bbox_connecting_v12[0] = bbox_facing_v1;
        bbox_connecting_v12[1] = bbox_connecting_v12[0];

        for (int i = 0; i < 2; i++) {
            std::vector<int> pair;
            set_intersection(tet_vertices[tets[t_id][mod4(j12[mod2(i + 1)] + 1)]].conn_tets,
                             tet_vertices[tets[t_id][mod4(j12[mod2(i + 1)] + 2)]].conn_tets,
                             tet_vertices[tets[t_id][mod4(j12[mod2(i + 1)] + 3)]].conn_tets, pair);
//            if(!(pair.size() == 1 || pair.size() == 2)) {
//                cout << "!(pair.size() == 1 || pair.size() == 2)" << endl;
//                cout << "******"<<mod4(j12[mod2(i + 1)] + 1) << endl;
//                vector_print(tet_vertices[tets[t_id][mod4(j12[mod2(i + 1)] + 1)]].conn_tets, " ");
//                cout << endl;
//                cout << "******"<<mod4(j12[mod2(i + 1)] + 2) << endl;
//                vector_print(tet_vertices[tets[t_id][mod4(j12[mod2(i + 1)] + 2)]].conn_tets, " ");
//                cout << endl;
//                cout << "******"<<mod4(j12[mod2(i + 1)] + 3) << endl;
//                vector_print(tet_vertices[tets[t_id][mod4(j12[mod2(i + 1)] + 3)]].conn_tets, " ");
//                cout << endl;
//                cout << "******" << endl;
//                vector_print(pair, " ");
//                //pausee();
//            }
            if (pair.size() > 1) {
                int opp_t_id = pair[0] == t_id ? pair[1] : pair[0];
                for (int j = 0; j < 4; j++) {
                    if (tets[opp_t_id][j] != tets[t_id][mod4(j12[mod2(i + 1)] + 1)]
                        && tets[opp_t_id][j] != tets[t_id][mod4(j12[mod2(i + 1)] + 2)]
                        && tets[opp_t_id][j] != tets[t_id][mod4(j12[mod2(i + 1)] + 3)]) {
                        tets[opp_t_id].is_surface_fs[j] = sf_connecting_v12[i];
                        tets[opp_t_id].surface_tags[j] = tag_connecting_v12[i];
                        tets[opp_t_id].is_bbox_fs[j] = bbox_connecting_v12[i];
                        break;
                    }
                }
            }
        }
    }

    //update connectivity
    ts++;
    ii = 0;
    for (int t_id:n1_t_ids) {
        int j = js_n1_t_ids[ii++];
        tets[t_id][j] = v2_id;
//        tet_vertices[v2_id].conn_tets.insert(t_id);
        tet_vertices[v2_id].conn_tets.push_back(t_id);
        if(is_update_tss)
            tet_tss[t_id] = ts;//update timestamp
    }
    for (int t_id: n12_t_ids) {
        tets[t_id].is_removed = true;
        for (int j = 0; j < 4; j++) {
            if (tets[t_id][j] != v1_id)
//                tet_vertices[tets[t_id][j]].conn_tets.erase(t_id);
                vector_erase(tet_vertices[tets[t_id][j]].conn_tets, t_id);
        }
    }

    tet_vertices[v1_id].conn_tets.clear();

    ////re-push
    for (int v_id:n1_v_ids) {
        if(v_id!=v1_id)
            new_edges.push_back({{v2_id, v_id}});
    }

//    for (int v_id:n_v_ids) {
//        new_edges.push_back({{v2_id, v_id}});
//    }

    if(tet_vertices[v1_id].is_on_surface)
        return EC_SUCCESS_ENVELOPE;
    return EC_SUCCESS;
}

bool floatTetWild::is_edge_freezed(Mesh& mesh, int v1_id, int v2_id){
    if(mesh.tet_vertices[v1_id].is_freezed || mesh.tet_vertices[v2_id].is_freezed)
        return true;
    return false;
}

bool floatTetWild::is_collapsable_bbox(Mesh& mesh, int v1_id, int v2_id) {
    if (!mesh.tet_vertices[v1_id].is_on_bbox)
        return true;
    else if (!mesh.tet_vertices[v2_id].is_on_bbox)
        return false;

//    std::unordered_set<int> bbox_fs2;
//    for (int t_id:mesh.tet_vertices[v2_id].conn_tets) {
//        for (int j = 0; j < 4; j++) {
//            if (mesh.tets[t_id][j] != v2_id && mesh.tets[t_id].is_bbox_fs[j] != NOT_BBOX)
//                bbox_fs2.insert(mesh.tets[t_id].is_bbox_fs[j]);
//        }
//    }
//    int old_size = bbox_fs2.size();
//    for (int t_id:mesh.tet_vertices[v1_id].conn_tets) {
//        for (int j = 0; j < 4; j++) {
//            if (mesh.tets[t_id][j] != v1_id && mesh.tets[t_id].is_bbox_fs[j] != NOT_BBOX) {
//                bbox_fs2.insert(mesh.tets[t_id].is_bbox_fs[j]);
//                if (bbox_fs2.size() > old_size)
//                    return false;
//            }
//        }
//    }

    std::vector<int> bbox_fs2;
    for (int t_id:mesh.tet_vertices[v2_id].conn_tets) {
        for (int j = 0; j < 4; j++) {
            if (mesh.tets[t_id][j] != v2_id && mesh.tets[t_id].is_bbox_fs[j] != NOT_BBOX)
                bbox_fs2.push_back(mesh.tets[t_id].is_bbox_fs[j]);
        }
    }
    vector_unique(bbox_fs2);

    for (int t_id:mesh.tet_vertices[v1_id].conn_tets) {
        for (int j = 0; j < 4; j++) {
            if (mesh.tets[t_id][j] != v1_id && mesh.tets[t_id].is_bbox_fs[j] != NOT_BBOX) {
                if(std::find(bbox_fs2.begin(), bbox_fs2.end(), mesh.tets[t_id].is_bbox_fs[j]) == bbox_fs2.end())
                    return false;
            }
        }
    }

    return true;
}

bool floatTetWild::is_collapsable_length(Mesh& mesh, int v1_id, int v2_id, Scalar l_2) {
    Scalar sizing_scalar = (mesh.tet_vertices[v1_id].sizing_scalar + mesh.tet_vertices[v2_id].sizing_scalar) / 2;
    if (l_2 <= mesh.params.collapse_threshold_2 * sizing_scalar * sizing_scalar)
        return true;
    return false;
}

bool floatTetWild::is_collapsable_boundary(Mesh& mesh, int v1_id, int v2_id, const AABBWrapper& tree) {
    if (mesh.tet_vertices[v1_id].is_on_boundary && !is_boundary_edge(mesh, v1_id, v2_id, tree))
        return false;
    return true;


//    if (mesh.tet_vertices[v1_id].on_boundary_e_id >= 0 && mesh.tet_vertices[v2_id].on_boundary_e_id
//        && mesh.tet_vertices[v1_id].on_boundary_e_id != mesh.tet_vertices[v2_id].on_boundary_e_id)
//        return false;
//    return true;
}