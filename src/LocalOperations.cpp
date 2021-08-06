// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include <floattetwild/LocalOperations.h>
#include <floattetwild/Predicates.hpp>

#include <igl/Timer.h>

#ifdef FLOAT_TETWILD_USE_TBB
#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include <tbb/atomic.h>
#include <tbb/concurrent_vector.h>
#include <tbb/parallel_sort.h>
#endif

namespace floatTetWild {
    bool use_old_energy = false;
    std::string envelope_log_csv = "";
    int envelope_log_csv_cnt = 0;
}

using floatTetWild::Scalar;

// void floatTetWild::init_b_tree(const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces,
//         GEO::Mesh& b_mesh) {
// //    std::vector<std::array<int, 2>> edges;
// //    for(int i=0;i<sf_mesh.facets.nb();i++){
// //        for(int j=0;j<3;j++) {
// //            if(sf_mesh.facets.adjacent(i, j)==GEO::NO_FACET){
// //                edges.push_back({{}})
// //            }
// //        }
// //    }

//     std::vector<std::vector<int>> conn_tris(input_vertices.size());
//     std::vector<std::array<int, 2>> all_edges;
//     all_edges.reserve(input_faces.size() * 3);
//     for (int i = 0; i < input_faces.size(); i++) {
//         for (int j = 0; j < 3; j++) {
//             conn_tris[input_faces[i][j]].push_back(i);
//             if (input_faces[i][j] < input_faces[i][(j + 1) % 3])
//                 all_edges.push_back({{input_faces[i][j], input_faces[i][(j + 1) % 3]}});
//             else
//                 all_edges.push_back({{input_faces[i][(j + 1) % 3], input_faces[i][j]}});
//         }
//     }
//     vector_unique(all_edges);

//     std::vector<std::array<int, 2>> b_edges;
//     for (auto &e:all_edges) {
//         std::vector<int> tmp;
//         std::set_intersection(conn_tris[e[0]].begin(), conn_tris[e[0]].end(),
//                               conn_tris[e[1]].begin(), conn_tris[e[1]].end(), std::back_inserter(tmp));
//         if (tmp.size() == 1) {
//             b_edges.push_back(e);
//         }
//     }

//     if (b_edges.empty()) {
//         b_mesh.vertices.clear();
//         b_mesh.vertices.create_vertices(1);
//         b_mesh.vertices.point(0) = GEO::vec3(0, 0, 0);
//         b_mesh.facets.clear();
//         b_mesh.facets.create_triangles(1);
//         b_mesh.facets.set_vertex(0, 0, 0);
//         b_mesh.facets.set_vertex(0, 1, 0);
//         b_mesh.facets.set_vertex(0, 2, 0);
//     } else {
//         b_mesh.vertices.clear();
//         b_mesh.vertices.create_vertices((int) b_edges.size() * 2);
//         int cnt = 0;
//         for (auto &e:b_edges) {
//             for (int j = 0; j < 2; j++) {
//                 GEO::vec3 &p = b_mesh.vertices.point(cnt++);
//                 p[0] = input_vertices[e[j]][0];
//                 p[1] = input_vertices[e[j]][1];
//                 p[2] = input_vertices[e[j]][2];
//             }
//         }
//         b_mesh.facets.clear();
//         b_mesh.facets.create_triangles((int) b_edges.size());
//         for (int i = 0; i < b_edges.size(); i++) {
//             b_mesh.facets.set_vertex(i, 0, i * 2);
//             b_mesh.facets.set_vertex(i, 1, i * 2);
//             b_mesh.facets.set_vertex(i, 2, i * 2 + 1);
//         }
//     }
// }

int floatTetWild::get_opp_t_id(const Mesh& mesh, int t_id, int j) {
    std::vector<int> pair;
    set_intersection(mesh.tet_vertices[mesh.tets[t_id][(j + 1) % 4]].conn_tets,
                     mesh.tet_vertices[mesh.tets[t_id][(j + 2) % 4]].conn_tets,
                     mesh.tet_vertices[mesh.tets[t_id][(j + 3) % 4]].conn_tets, pair);
    if (pair.size() == 2)
        return pair[0] == t_id ? pair[1] : pair[0];
    return OPP_T_ID_BOUNDARY;
}

void floatTetWild::set_opp_t_id(Mesh& mesh, int t_id, int j){
    auto& t = mesh.tets[t_id];
//    static double time = 0;
    const int jp1 = mod4(j+1);
    const int jp2 = mod4(j+2);
    const int jp3 = mod4(j+3);
    assert((j + 1) % 4 == jp1);
    assert((j + 2) % 4 == jp2);
    assert((j + 3) % 4 == jp3);
//    igl::Timer timer;
//    timer.start();
//    std::unordered_set<int> tmp;
//    set_intersection(mesh.tet_vertices[t[(j + 1) % 4]].conn_tets,
//                     mesh.tet_vertices[t[(j + 2) % 4]].conn_tets, tmp);
    static std::vector<int> pair;
    pair.clear();
//    set_intersection(mesh.tet_vertices[t[(j + 3) % 4]].conn_tets, tmp, pair);
    set_intersection(mesh.tet_vertices[t[jp1]].conn_tets, mesh.tet_vertices[t[jp2]].conn_tets, mesh.tet_vertices[t[jp3]].conn_tets, pair);
//    timer.stop();
//    time+=timer.getElapsedTimeInSec();
//    std::cout<<"set_opp_t_id "<<time<<std::endl;
    if (pair.size() == 2) {
        int opp_t_id = pair[0] == t_id ? pair[1] : pair[0];
        t.opp_t_ids[j] = opp_t_id;
        auto &opp_t = mesh.tets[opp_t_id];
        for (int k = 0; k < 4; k++) {
            if (opp_t[k] != t[jp1] && opp_t[k] != t[jp2] && opp_t[k] != t[jp3]) {
                opp_t.opp_t_ids[k] = t_id;
                break;
            }
        }
    }
}

void floatTetWild::get_all_edges(const Mesh& mesh, std::vector<std::array<int, 2>>& edges){
    edges.reserve(mesh.tets.size()*6);

#ifdef FLOAT_TETWILD_USE_TBB
    tbb::concurrent_vector<std::array<int, 2>> edges_tbb;
    tbb::parallel_for( size_t(0), mesh.tets.size(), [&](size_t i)
#else
    for (unsigned int i = 0; i < mesh.tets.size(); i++)
#endif
        {
        if (mesh.tets[i].is_removed){
#ifdef FLOAT_TETWILD_USE_TBB
            return;
#else
            continue;
#endif
        }
        for (int j = 0; j < 3; j++) {
            std::array<int, 2> e = {{mesh.tets[i][0], mesh.tets[i][j + 1]}};
            if (e[0] > e[1])
                std::swap(e[0], e[1]);
#ifdef FLOAT_TETWILD_USE_TBB
            edges_tbb.push_back(e);
#else
            edges.push_back(e);
#endif
            e = {{mesh.tets[i][j + 1], mesh.tets[i][mod3(j + 1) + 1]}};
            if (e[0] > e[1])
                std::swap(e[0], e[1]);
#ifdef FLOAT_TETWILD_USE_TBB
            edges_tbb.push_back(e);
#else
            edges.push_back(e);
#endif
        }
    }
#ifdef FLOAT_TETWILD_USE_TBB
    );
    edges.reserve(edges_tbb.size());
    edges.insert(edges.end(), edges_tbb.begin(), edges_tbb.end());
    assert(edges_tbb.size() == edges.size());
    tbb::parallel_sort(edges.begin(), edges.end());

    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
#else
    vector_unique(edges);
#endif
}

void floatTetWild::get_all_edges(const Mesh& mesh, const std::vector<int>& t_ids, std::vector<std::array<int, 2>>& edges, bool skip_freezed) {
    for (unsigned int i = 0; i < t_ids.size(); i++) {
        auto &t = mesh.tets[t_ids[i]];
        for (int j = 0; j < 3; j++) {
            if (skip_freezed) {
                if (!mesh.tet_vertices[t[0]].is_freezed && !mesh.tet_vertices[t[j + 1]].is_freezed) {
                    std::array<int, 2> e = {{t[0], t[j + 1]}};
                    if (e[0] > e[1])
                        std::swap(e[0], e[1]);
                    edges.push_back(e);
                }
                if (!mesh.tet_vertices[t[j + 1]].is_freezed && !mesh.tet_vertices[mod3(j + 1) + 1].is_freezed) {
                    std::array<int, 2> e = {{t[j + 1], t[mod3(j + 1) + 1]}};
                    if (e[0] > e[1])
                        std::swap(e[0], e[1]);
                    edges.push_back(e);
                }
            } else {
                std::array<int, 2> e = {{t[0], t[j + 1]}};
                if (e[0] > e[1])
                    std::swap(e[0], e[1]);
                edges.push_back(e);
                e = {{t[j + 1], t[mod3(j + 1) + 1]}};
                if (e[0] > e[1])
                    std::swap(e[0], e[1]);
                edges.push_back(e);
            }
        }
    }
    vector_unique(edges);
}

Scalar floatTetWild::get_edge_length(const Mesh& mesh, int v1_id, int v2_id) {
    return (mesh.tet_vertices[v1_id].pos - mesh.tet_vertices[v2_id].pos).norm();
}

Scalar floatTetWild::get_edge_length_2(const Mesh& mesh, int v1_id, int v2_id) {
    return (mesh.tet_vertices[v1_id].pos - mesh.tet_vertices[v2_id].pos).squaredNorm();
}

bool floatTetWild::is_bbox_edge(const Mesh& mesh, int v1_id, int v2_id, const std::vector<int>& n12_t_ids) {
    if (!mesh.tet_vertices[v1_id].is_on_bbox || !mesh.tet_vertices[v2_id].is_on_bbox)
        return false;

    for (int t_id:n12_t_ids) {
        for (int j = 0; j < 4; j++) {
            if (mesh.tets[t_id][j] != v1_id && mesh.tets[t_id][j] != v2_id
                && mesh.tets[t_id].is_bbox_fs[j] != NOT_BBOX)
                return true;
        }
    }
    return false;
}

bool floatTetWild::is_surface_edge(const Mesh& mesh, int v1_id, int v2_id, const std::vector<int>& n12_t_ids){
    if (!mesh.tet_vertices[v1_id].is_on_surface || !mesh.tet_vertices[v2_id].is_on_surface)
        return false;

    for (int t_id:n12_t_ids) {
        for (int j = 0; j < 4; j++) {
            if (mesh.tets[t_id][j] != v1_id && mesh.tets[t_id][j] != v2_id
                && mesh.tets[t_id].is_surface_fs[j] != NOT_SURFACE)
                return true;
        }
    }
    return false;
}

bool floatTetWild::is_boundary_edge(const Mesh& mesh, int v1_id, int v2_id, const AABBWrapper& tree) {
    if (!mesh.tet_vertices[v1_id].is_on_boundary || !mesh.tet_vertices[v2_id].is_on_boundary)
        return false;

#ifdef NEW_ENVELOPE
    if(!mesh.is_input_all_inserted) {
        return !tree.is_out_tmp_b_envelope_exact({{mesh.tet_vertices[v1_id].pos, mesh.tet_vertices[v2_id].pos,
                                                          mesh.tet_vertices[v2_id].pos}});
    } else {
        return !tree.is_out_b_envelope_exact({{mesh.tet_vertices[v1_id].pos, mesh.tet_vertices[v2_id].pos,
                                                mesh.tet_vertices[v2_id].pos}});
    }
#else
    std::vector<GEO::vec3> ps;
    ps.push_back(GEO::vec3(mesh.tet_vertices[v1_id].pos[0], mesh.tet_vertices[v1_id].pos[1],
            mesh.tet_vertices[v1_id].pos[2]));
    int p0_id = 0;
    Scalar l = get_edge_length(mesh, v1_id, v2_id);
    int N = l / mesh.params.dd + 1;
    ps.push_back(GEO::vec3(mesh.tet_vertices[v2_id][0], mesh.tet_vertices[v2_id][1],
                           mesh.tet_vertices[v2_id][2]));
    int p1_id = ps.size() - 1;
    for (Scalar j = 1; j < N - 1; j++) {
        ps.push_back(ps[p0_id] * (j / N) + ps[p1_id] * (1 - j / N));
    }

    if(!mesh.is_input_all_inserted) {
        return !tree.is_out_tmp_b_envelope(ps, mesh.params.eps_2);
    } else {
        return !tree.is_out_b_envelope(ps, mesh.params.eps_2);
    }
#endif

//    if(!mesh.is_input_all_inserted)
//        return true;
//
//    int cnt = 0;
//    for (int t_id: mesh.tet_vertices[v1_id].conn_tets) {
//        std::array<int, 4> opp_js;
//        int ii = 0;
//        for (int j = 0; j < 4; j++) {
//            if (mesh.tets[t_id][j] == v1_id || mesh.tets[t_id][j] == v2_id)
//                continue;
//            opp_js[ii++] = j;
//        }
//        if (ii == 2) {
//            if (mesh.tets[t_id].is_surface_fs[opp_js[0]] != NOT_SURFACE)
//                cnt++;
//            if (mesh.tets[t_id].is_surface_fs[opp_js[1]] != NOT_SURFACE)
//                cnt++;
//            if (cnt > 2)
//                return false;
//        }
//    }
//    if (cnt == 2)
//        return true;
//    return false;
}

bool floatTetWild::is_valid_edge(const Mesh& mesh, int v1_id, int v2_id) {
    if (mesh.tet_vertices[v1_id].is_removed || mesh.tet_vertices[v2_id].is_removed)
        return false;
//    std::vector<int> tmp;
//    set_intersection(mesh.tet_vertices[v1_id].conn_tets, mesh.tet_vertices[v2_id].conn_tets, tmp);
//    if (tmp.empty()) {
//        cout<<"happen"<<endl;
//        //pausee();
//        return false;
//    }

    return true;
}

bool floatTetWild::is_valid_edge(const Mesh& mesh, int v1_id, int v2_id, const std::vector<int>& n12_t_ids) {
    if (mesh.tet_vertices[v1_id].is_removed || mesh.tet_vertices[v2_id].is_removed)
        return false;
    if (n12_t_ids.empty())
        return false;

    return true;
}

bool floatTetWild::is_isolate_surface_point(const Mesh& mesh, int v_id) {
    for (int t_id:mesh.tet_vertices[v_id].conn_tets) {
        for (int j = 0; j < 4; j++) {
            if (mesh.tets[t_id][j] != v_id && mesh.tets[t_id].is_surface_fs[j] != NOT_SURFACE)
                return false;
        }
    }

    return true;
}

bool floatTetWild::is_point_out_envelope(const Mesh& mesh, const Vector3& p, const AABBWrapper& tree){
#ifdef NEW_ENVELOPE
    return tree.is_out_sf_envelope_exact(p);
#else
    GEO::index_t prev_facet;
    return tree.is_out_sf_envelope(p, mesh.params.eps_2, prev_facet);
#endif
//    GEO::vec3 geo_p(p[0], p[1], p[2]);
//    if (sf_tree.squared_distance(geo_p) > mesh.params.eps_2)
//        return true;
//
//    return false;
}

bool floatTetWild::is_point_out_boundary_envelope(const Mesh& mesh, const Vector3& p, const AABBWrapper& tree){
    if(mesh.is_input_all_inserted)
        return false;

    GEO::index_t prev_facet;
    return tree.is_out_tmp_b_envelope(p, mesh.params.eps_2, prev_facet);

//    GEO::vec3 geo_p(p[0], p[1], p[2]);
//    if (b_tree.squared_distance(geo_p) > mesh.params.eps_2)
//        return true;
//
//    return false;
}

Scalar floatTetWild::get_quality(const Mesh& mesh, const MeshTet& t) {
    std::array<Scalar, 12> T;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 3; j++)
            T[i * 3 + j] = mesh.tet_vertices[t[i]].pos[j];
    }

    return AMIPS_energy(T);
//    Scalar q = AMIPS_energy(T);
//    if (q > 1e8)
//        return MAX_ENERGY;
//    else
//        return q;
}

Scalar floatTetWild::get_quality(const Mesh& mesh, int t_id) {
    std::array<Scalar, 12> T;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 3; j++)
            T[i * 3 + j] = mesh.tet_vertices[mesh.tets[t_id][i]].pos[j];
    }
    return AMIPS_energy(T);
//    Scalar q = AMIPS_energy(T);
//    if (q > 1e8)
//        return MAX_ENERGY;
//    else
//        return q;
}

Scalar floatTetWild::get_quality(const MeshVertex& v0, const MeshVertex& v1, const MeshVertex& v2, const MeshVertex& v3) {
    std::array<Scalar, 12> T = {{v0.pos[0], v0.pos[1], v0.pos[2], v1.pos[0], v1.pos[1], v1.pos[2],
                                        v2.pos[0], v2.pos[1], v2.pos[2], v3.pos[0], v3.pos[1], v3.pos[2]}};
    return AMIPS_energy(T);
//    Scalar q = AMIPS_energy(T);
//    if (q > 1e8)
//        return MAX_ENERGY;
//    else
//        return q;
}

Scalar floatTetWild::get_quality(const Vector3& v0, const Vector3& v1, const Vector3& v2, const Vector3& v3){
    std::array<Scalar, 12> T = {{v0[0], v0[1], v0[2], v1[0], v1[1], v1[2],
                                        v2[0], v2[1], v2[2], v3[0], v3[1], v3[2]}};
    return AMIPS_energy(T);
//    Scalar q = AMIPS_energy(T);
//    if (q > 1e8)
//        return MAX_ENERGY;
//    else
//        return q;
}

void floatTetWild::get_max_avg_energy(const Mesh& mesh, Scalar& max_energy, Scalar& avg_energy) {
    max_energy = 0;
    avg_energy = 0;
    int cnt = 0;
    for (auto &t: mesh.tets) {
        if (t.is_removed)
            continue;
        if (t.quality > max_energy)
            max_energy = t.quality;
        avg_energy += t.quality;
        cnt++;
    }
    avg_energy /= cnt;
}

Scalar floatTetWild::get_mid_energy(const Mesh& mesh){
    std::vector<Scalar > tmp;
    for (auto& t:mesh.tets) {
        if (t.is_removed)
            continue;
        tmp.push_back(t.quality);
    }
    std::sort(tmp.begin(), tmp.end());
    return tmp[tmp.size() / 2];
}

bool floatTetWild::is_inverted(const Mesh& mesh, int t_id) {
    if (Predicates::orient_3d(mesh.tet_vertices[mesh.tets[t_id][0]].pos, mesh.tet_vertices[mesh.tets[t_id][1]].pos,
                              mesh.tet_vertices[mesh.tets[t_id][2]].pos, mesh.tet_vertices[mesh.tets[t_id][3]].pos) ==
        Predicates::ORI_POSITIVE)
        return false;
    return true;
}

bool floatTetWild::is_inverted(const Mesh& mesh, int t_id, int j, const Vector3& new_p) {
    int ori;
    if (j == 0) {
        ori = Predicates::orient_3d(new_p, mesh.tet_vertices[mesh.tets[t_id][1]].pos,
                                    mesh.tet_vertices[mesh.tets[t_id][2]].pos,
                                    mesh.tet_vertices[mesh.tets[t_id][3]].pos);
    } else if (j == 1) {
        ori = Predicates::orient_3d(mesh.tet_vertices[mesh.tets[t_id][0]].pos, new_p,
                                    mesh.tet_vertices[mesh.tets[t_id][2]].pos,
                                    mesh.tet_vertices[mesh.tets[t_id][3]].pos);
    } else if (j == 2) {
        ori = Predicates::orient_3d(mesh.tet_vertices[mesh.tets[t_id][0]].pos,
                                    mesh.tet_vertices[mesh.tets[t_id][1]].pos, new_p,
                                    mesh.tet_vertices[mesh.tets[t_id][3]].pos);
    } else {
        ori = Predicates::orient_3d(mesh.tet_vertices[mesh.tets[t_id][0]].pos,
                                    mesh.tet_vertices[mesh.tets[t_id][1]].pos,
                                    mesh.tet_vertices[mesh.tets[t_id][2]].pos, new_p);
    }
    if (ori == Predicates::ORI_POSITIVE)
        return false;
    return true;
}

bool floatTetWild::is_inverted(const MeshVertex& v0, const MeshVertex& v1, const MeshVertex& v2, const MeshVertex& v3){
    if (Predicates::orient_3d(v0.pos, v1.pos, v2.pos, v3.pos) == Predicates::ORI_POSITIVE)
        return false;
    return true;
}

bool floatTetWild::is_inverted(const Vector3& v0, const Vector3& v1, const Vector3& v2, const Vector3& v3){
    if (Predicates::orient_3d(v0, v1, v2, v3) == Predicates::ORI_POSITIVE)
        return false;
    return true;
}

bool floatTetWild::is_degenerate(const Vector3& v0, const Vector3& v1, const Vector3& v2, const Vector3& v3){
    if (Predicates::orient_3d(v0, v1, v2, v3) == Predicates::ORI_ZERO)
        return true;
    return false;
}

bool floatTetWild::is_out_boundary_envelope(const Mesh& mesh, int v_id, const Vector3& new_pos, const AABBWrapper& tree){
    if(mesh.is_input_all_inserted)
        return false;
    if(!mesh.tet_vertices[v_id].is_on_cut)
        return false;

    GEO::index_t prev_facet;
    if(tree.is_out_tmp_b_envelope(new_pos, mesh.params.eps_2/100, prev_facet))
        return true;

    std::vector<int> tmp_b_v_ids;
    for (int t_id:mesh.tet_vertices[v_id].conn_tets) {
        for (int j = 0; j < 4; j++) {
            if (mesh.tets[t_id][j] != v_id && mesh.tets[t_id].is_surface_fs[j] <= 0) {
                for(int k=0;k<3;k++){
                    int b_v_id = mesh.tets[t_id][(j+1+k)%4];
                    if(b_v_id != v_id && mesh.tet_vertices[b_v_id].is_on_boundary)
                        tmp_b_v_ids.push_back(b_v_id);
                }
            }
        }
    }
    vector_unique(tmp_b_v_ids);

    std::vector<int> b_v_ids;
    b_v_ids.reserve(tmp_b_v_ids.size());
    for(int b_v_id:tmp_b_v_ids){
        if(is_boundary_edge(mesh, v_id, b_v_id, tree))
            b_v_ids.push_back(b_v_id);
    }
    if(b_v_ids.empty())
        return false;

    std::vector<GEO::vec3> ps;
    ps.push_back(GEO::vec3(new_pos[0], new_pos[1], new_pos[2]));
    int p0_id = 0;
    for(int b_v_id:b_v_ids) {
        Scalar l = get_edge_length(mesh, v_id, b_v_id);
        int N = l / mesh.params.dd + 1;
        ps.push_back(GEO::vec3(mesh.tet_vertices[b_v_id][0], mesh.tet_vertices[b_v_id][1],
                               mesh.tet_vertices[b_v_id][2]));
        int p1_id = ps.size() - 1;
        for (Scalar j = 1; j < N - 1; j++) {
//            ps.push_back(ps[0] * (j / N) + ps[1] * (1 - j / N));
            ps.push_back(ps[p0_id] * (j / N) + ps[p1_id] * (1 - j / N));
        }
    }

    return tree.is_out_tmp_b_envelope(ps, mesh.params.eps_2/100, prev_facet);

//    GEO::vec3 init_point(new_pos[0], new_pos[1], new_pos[2]);
//    GEO::vec3 nearest_point;
//    double sq_distg;
//    GEO::index_t prev_facet = b_tree.nearest_facet(init_point, nearest_point, sq_distg);
//    Scalar sq_dist = sq_distg;
//    if(sq_dist > mesh.params.eps_2)
//        return true;
//
//    std::vector<int> tmp_b_v_ids;
//    for (int t_id:mesh.tet_vertices[v_id].conn_tets) {
//        for (int j = 0; j < 4; j++) {
//            if (mesh.tets[t_id][j] != v_id && mesh.tets[t_id].is_surface_fs[j] < 0) {
//                for(int k=0;k<3;k++){
//                    int b_v_id = mesh.tets[t_id][(j+1+k)%4];
//                    if(b_v_id != v_id && mesh.tet_vertices[b_v_id].is_on_boundary)
//                        tmp_b_v_ids.push_back(b_v_id);
//                }
//            }
//        }
//    }
//    vector_unique(tmp_b_v_ids);
//
//    std::vector<int> b_v_ids;
//    b_v_ids.reserve(tmp_b_v_ids.size());
//    for(int b_v_id:tmp_b_v_ids){
//        if(is_boundary_edge(mesh, v_id, b_v_id))//todo: can be improved, see meshimprovemnet init()
//            b_v_ids.push_back(b_v_id);
//    }
//    if(b_v_ids.empty())
//        return false;
//
//    std::vector<GEO::vec3> ps;
//    for(int b_v_id:b_v_ids) {
//        Scalar l = get_edge_length(mesh, v_id, b_v_id);
//        int N = l / mesh.params.dd + 1;
////        ps.push_back(GEO::vec3(mesh.tet_vertices[v_id][0], mesh.tet_vertices[v_id][1], mesh.tet_vertices[v_id][2]));
//        ps.push_back(GEO::vec3(new_pos[0], new_pos[1], new_pos[2]));
//        ps.push_back(
//                GEO::vec3(mesh.tet_vertices[b_v_id][0], mesh.tet_vertices[b_v_id][1], mesh.tet_vertices[b_v_id][2]));
//        for (Scalar j = 0; j < N - 1; j++) {
//            ps.push_back(ps[0] * (j / N) + ps[1] * (1 - j / N));
//        }
//    }
//
//    int cnt = 0;
//    const unsigned int ps_size = ps.size();
//    for (unsigned int i = ps_size / 2; ; i = (i + 1) % ps_size) {//check from the middle
//        GEO::vec3 &current_point = ps[i];
//        sq_distg = current_point.distance2(nearest_point);
//        b_tree.nearest_facet_with_hint(current_point, prev_facet, nearest_point, sq_distg);
//        sq_dist = sq_distg;
//        if (sq_dist > mesh.params.eps_2)
//            return true;
//        cnt++;
//        if (cnt >= ps_size)
//            break;
//    }
//
//    return false;
}

#include <sstream>
bool floatTetWild::is_out_envelope(Mesh& mesh, int v_id, const Vector3& new_pos, const AABBWrapper& tree) {
#ifdef NEW_ENVELOPE
    if(tree.is_out_sf_envelope_exact(new_pos))
        return true;
#else
    GEO::index_t prev_facet;
    if(tree.is_out_sf_envelope(new_pos, mesh.params.eps_2, prev_facet))
        return true;
#endif

    std::vector<GEO::vec3> ps;
    for (int t_id:mesh.tet_vertices[v_id].conn_tets) {
        for (int j = 0; j < 4; j++) {
            if (mesh.tets[t_id][j] != v_id && mesh.tets[t_id].is_surface_fs[j] <= 0) {
                std::array<Vector3, 3> vs;
                for (int k = 0; k < 3; k++) {
                    if (mesh.tets[t_id][mod4(j + 1 + k)] == v_id)
                        vs[k] = new_pos;
                    else
                        vs[k] = mesh.tet_vertices[mesh.tets[t_id][mod4(j + 1 + k)]].pos;
                }
#ifdef NEW_ENVELOPE
                bool is_out =tree.is_out_sf_envelope_exact(vs);
                if(!mesh.params.envelope_log.empty()){
                    if(envelope_log_csv_cnt < 1e5) {
                        std::ostringstream ss;
                        ss << std::setprecision(17);
                        for (const auto &v: vs) {
                            ss << v[0] << ',' << v[1] << ',' << v[2] << ',';
                        }
                        ss << is_out << "\n";
                        std::string tmp = ss.str();
                        envelope_log_csv += tmp;
                        envelope_log_csv_cnt += 1;
                    } else {
                        std::ofstream fout(mesh.params.envelope_log);
                        fout << envelope_log_csv;
                        fout.close();
                        mesh.params.envelope_log = "";
                    }
                }
                if (is_out)
                    return true;
#else
    #ifdef STORE_SAMPLE_POINTS
                    ps.clear();
                    sample_triangle(vs, ps, mesh.params.dd);
                    bool is_out = tree.is_out_sf_envelope(ps, mesh.params.eps_2, prev_facet);
    #else
                    bool is_out = sample_triangle_and_check_is_out(vs, mesh.params.dd, mesh.params.eps_2, tree, prev_facet);
    #endif
                    if(!mesh.params.envelope_log.empty()){
                        if(envelope_log_csv_cnt < 1e5) {
                            std::ostringstream ss;
                            ss << std::setprecision(17);
                            for (const auto &v: vs) {
                                ss << v[0] << ',' << v[1] << ',' << v[2] << ',';
                            }
                            ss << is_out << "\n";
                            std::string tmp = ss.str();
                            envelope_log_csv += tmp;
                            envelope_log_csv_cnt += 1;
                        } else {
                            std::ofstream fout(mesh.params.envelope_log);
                            fout << envelope_log_csv;
                            fout.close();
                            mesh.params.envelope_log = "";
                        }
                    }
                    if (is_out)
                        return true;
#endif

//                int cnt = 0;
//                const unsigned int ps_size = ps.size();
//                for (unsigned int i = ps_size / 2; ; i = (i + 1) % ps_size) {//check from the middle
//                    GEO::vec3 &current_point = ps[i];
//                    sq_distg = current_point.distance2(nearest_point);
//                    sf_tree.nearest_facet_with_hint(current_point, prev_facet, nearest_point, sq_distg);
//                    sq_dist = sq_distg;
//                    if (sq_dist > mesh.params.eps_2)
//                        return true;
//                    cnt++;
//                    if (cnt >= ps_size)
//                        break;
//                }
            }
        }
    }

    return false;

//    GEO::vec3 init_point(new_pos[0], new_pos[1], new_pos[2]);
//    GEO::vec3 nearest_point;
//    double sq_distg;
//    GEO::index_t prev_facet = sf_tree.nearest_facet(init_point, nearest_point, sq_distg);
//    Scalar sq_dist = sq_distg;
//    if(sq_dist > mesh.params.eps_2)
//        return true;
//
//    std::vector<GEO::vec3> ps;
//    for (int t_id:mesh.tet_vertices[v_id].conn_tets) {
//        for (int j = 0; j < 4; j++) {
//            if (mesh.tets[t_id][j] != v_id && mesh.tets[t_id].is_surface_fs[j] < 0) {
//                std::array<Vector3, 3> vs;
//                for(int k=0;k<3;k++){
//                    if(mesh.tets[t_id][mod4(j + 1 + k)] == v_id)
//                        vs[k] = new_pos;
//                    else
//                        vs[k] = mesh.tet_vertices[mesh.tets[t_id][(j + 1 + k) % 4]].pos;
//                }
//
//                ps.clear();
////                sample_triangle({{mesh.tet_vertices[mesh.tets[t_id][(j + 1) % 4]].pos,
////                                         mesh.tet_vertices[mesh.tets[t_id][(j + 2) % 4]].pos,
////                                         mesh.tet_vertices[mesh.tets[t_id][(j + 3) % 4]].pos}}, ps, mesh.params.dd);
//                sample_triangle(vs, ps, mesh.params.dd);
//
//                int cnt = 0;
//                const unsigned int ps_size = ps.size();
//                for (unsigned int i = ps_size / 2; ; i = (i + 1) % ps_size) {//check from the middle
//                    GEO::vec3 &current_point = ps[i];
//                    sq_distg = current_point.distance2(nearest_point);
//                    sf_tree.nearest_facet_with_hint(current_point, prev_facet, nearest_point, sq_distg);
//                    sq_dist = sq_distg;
//                    if (sq_dist > mesh.params.eps_2)
//                        return true;
//                    cnt++;
//                    if (cnt >= ps_size)
//                        break;
//                }
//            }
//        }
//    }
//
//    return false;
}

void floatTetWild::sample_triangle(const std::array<Vector3, 3>& vs, std::vector<GEO::vec3>& ps, Scalar sampling_dist) {
    Scalar sqrt3_2 = std::sqrt(3) / 2;

    std::array<Scalar, 3> ls;
    for (int i = 0; i < 3; i++) {
        ls[i] = (vs[i] - vs[mod3(i + 1)]).squaredNorm();
    }
    auto min_max = std::minmax_element(ls.begin(), ls.end());
    int min_i = min_max.first - ls.begin();
    int max_i = min_max.second - ls.begin();
    Scalar N = sqrt(ls[max_i]) / sampling_dist;
    if (N <= 1) {
        for (int i = 0; i < 3; i++)
            ps.push_back(GEO::vec3(vs[i][0], vs[i][1], vs[i][2]));
        return;
    }
    if (N == int(N))
        N -= 1;

    GEO::vec3 v0(vs[max_i][0], vs[max_i][1], vs[max_i][2]);
    GEO::vec3 v1(vs[mod3(max_i + 1)][0], vs[mod3(max_i + 1)][1], vs[mod3(max_i + 1)][2]);
    GEO::vec3 v2(vs[mod3(max_i + 2)][0], vs[mod3(max_i + 2)][1], vs[mod3(max_i + 2)][2]);

    GEO::vec3 n_v0v1 = GEO::normalize(v1 - v0);
    for (int n = 0; n <= N; n++) {
        ps.push_back(v0 + n_v0v1 * sampling_dist * n);
    }
    ps.push_back(v1);

    Scalar h = GEO::distance(GEO::dot((v2 - v0), (v1 - v0)) * (v1 - v0) / ls[max_i] + v0, v2);
    int M = h / (sqrt3_2 * sampling_dist);
    if (M < 1) {
        ps.push_back(v2);
        return;
    }

    GEO::vec3 n_v0v2 = GEO::normalize(v2 - v0);
    GEO::vec3 n_v1v2 = GEO::normalize(v2 - v1);
    Scalar tan_v0, tan_v1, sin_v0, sin_v1;
    sin_v0 = GEO::length(GEO::cross((v2 - v0), (v1 - v0))) / (GEO::distance(v0, v2) * GEO::distance(v0, v1));
    tan_v0 = GEO::length(GEO::cross((v2 - v0), (v1 - v0))) / GEO::dot((v2 - v0), (v1 - v0));
    tan_v1 = GEO::length(GEO::cross((v2 - v1), (v0 - v1))) / GEO::dot((v2 - v1), (v0 - v1));
    sin_v1 = GEO::length(GEO::cross((v2 - v1), (v0 - v1))) / (GEO::distance(v1, v2) * GEO::distance(v0, v1));

    for (int m = 1; m <= M; m++) {
        int n = sqrt3_2 / tan_v0 * m + 0.5;
        int n1 = sqrt3_2 / tan_v0 * m;
        if (m % 2 == 0 && n == n1) {
            n += 1;
        }
        GEO::vec3 v0_m = v0 + m * sqrt3_2 * sampling_dist / sin_v0 * n_v0v2;
        GEO::vec3 v1_m = v1 + m * sqrt3_2 * sampling_dist / sin_v1 * n_v1v2;
        if (GEO::distance(v0_m, v1_m) <= sampling_dist)
            break;

        Scalar delta_d = ((n + (m % 2) / 2.0) - m * sqrt3_2 / tan_v0) * sampling_dist;
        GEO::vec3 v = v0_m + delta_d * n_v0v1;
        int N1 = GEO::distance(v, v1_m) / sampling_dist;
//        ps.push_back(v0_m);
        for (int i = 0; i <= N1; i++) {
            ps.push_back(v + i * n_v0v1 * sampling_dist);
        }
//        ps.push_back(v1_m);
    }
    ps.push_back(v2);

    //sample edges
    N = sqrt(ls[mod3(max_i + 1)]) / sampling_dist;
    if (N > 1) {
        if (N == int(N))
            N -= 1;
        GEO::vec3 n_v1v2 = GEO::normalize(v2 - v1);
        for (int n = 1; n <= N; n++) {
            ps.push_back(v1 + n_v1v2 * sampling_dist * n);
        }
    }

    N = sqrt(ls[mod3(max_i + 2)]) / sampling_dist;
    if (N > 1) {
        if (N == int(N))
            N -= 1;
        GEO::vec3 n_v2v0 = GEO::normalize(v0 - v2);
        for (int n = 1; n <= N; n++) {
            ps.push_back(v2 + n_v2v0 * sampling_dist * n);
        }
    }
}

bool floatTetWild::sample_triangle_and_check_is_out(const std::array<Vector3, 3>& vs, Scalar sampling_dist,
        Scalar eps_2, const AABBWrapper& tree, GEO::index_t& prev_facet){
    GEO::vec3 nearest_point;
    double sq_dist = std::numeric_limits<double>::max();

    Scalar sqrt3_2 = std::sqrt(3) / 2;

    std::array<Scalar, 3> ls;
    for (int i = 0; i < 3; i++) {
        ls[i] = (vs[i] - vs[mod3(i + 1)]).squaredNorm();
    }
    auto min_max = std::minmax_element(ls.begin(), ls.end());
    int min_i = min_max.first - ls.begin();
    int max_i = min_max.second - ls.begin();
    Scalar N = sqrt(ls[max_i]) / sampling_dist;
    if (N <= 1) {
        for (int i = 0; i < 3; i++) {
//            ps.push_back(GEO::vec3(vs[i][0], vs[i][1], vs[i][2]));
            if (tree.is_out_sf_envelope(vs[i], eps_2, prev_facet, sq_dist, nearest_point))
                return true;
        }
//        return;
        return false;
    }
    if (N == int(N))
        N -= 1;

    GEO::vec3 v0(vs[max_i][0], vs[max_i][1], vs[max_i][2]);
    GEO::vec3 v1(vs[mod3(max_i + 1)][0], vs[mod3(max_i + 1)][1], vs[mod3(max_i + 1)][2]);
    GEO::vec3 v2(vs[mod3(max_i + 2)][0], vs[mod3(max_i + 2)][1], vs[mod3(max_i + 2)][2]);

    GEO::vec3 n_v0v1 = GEO::normalize(v1 - v0);
    for (int n = 0; n <= N; n++) {
//        ps.push_back(v0 + n_v0v1 * sampling_dist * n);
        if (tree.is_out_sf_envelope(v0 + n_v0v1 * sampling_dist * n, eps_2, prev_facet, sq_dist, nearest_point))
            return true;
    }
//    ps.push_back(v1);
    if (tree.is_out_sf_envelope(v1, eps_2, prev_facet, sq_dist, nearest_point))
        return true;

    Scalar h = GEO::distance(GEO::dot((v2 - v0), (v1 - v0)) * (v1 - v0) / ls[max_i] + v0, v2);
    int M = h / (sqrt3_2 * sampling_dist);
    if (M < 1) {
//        ps.push_back(v2);
//        return;
        return tree.is_out_sf_envelope(v2, eps_2, prev_facet, sq_dist, nearest_point);
    }

    GEO::vec3 n_v0v2 = GEO::normalize(v2 - v0);
    GEO::vec3 n_v1v2 = GEO::normalize(v2 - v1);
    Scalar tan_v0, tan_v1, sin_v0, sin_v1;
    sin_v0 = GEO::length(GEO::cross((v2 - v0), (v1 - v0))) / (GEO::distance(v0, v2) * GEO::distance(v0, v1));
    tan_v0 = GEO::length(GEO::cross((v2 - v0), (v1 - v0))) / GEO::dot((v2 - v0), (v1 - v0));
    tan_v1 = GEO::length(GEO::cross((v2 - v1), (v0 - v1))) / GEO::dot((v2 - v1), (v0 - v1));
    sin_v1 = GEO::length(GEO::cross((v2 - v1), (v0 - v1))) / (GEO::distance(v1, v2) * GEO::distance(v0, v1));

    for (int m = 1; m <= M; m++) {
        int n = sqrt3_2 / tan_v0 * m + 0.5;
        int n1 = sqrt3_2 / tan_v0 * m;
        if (m % 2 == 0 && n == n1) {
            n += 1;
        }
        GEO::vec3 v0_m = v0 + m * sqrt3_2 * sampling_dist / sin_v0 * n_v0v2;
        GEO::vec3 v1_m = v1 + m * sqrt3_2 * sampling_dist / sin_v1 * n_v1v2;
        if (GEO::distance(v0_m, v1_m) <= sampling_dist)
            break;

        Scalar delta_d = ((n + (m % 2) / 2.0) - m * sqrt3_2 / tan_v0) * sampling_dist;
        GEO::vec3 v = v0_m + delta_d * n_v0v1;
        int N1 = GEO::distance(v, v1_m) / sampling_dist;
        for (int i = 0; i <= N1; i++) {
//            ps.push_back(v + i * n_v0v1 * sampling_dist);
            if (tree.is_out_sf_envelope(v + i * n_v0v1 * sampling_dist, eps_2, prev_facet, sq_dist, nearest_point))
                return true;
        }
    }
//    ps.push_back(v2);
    if (tree.is_out_sf_envelope(v2, eps_2, prev_facet, sq_dist, nearest_point))
        return true;

    //sample edges
    N = sqrt(ls[mod3(max_i + 1)]) / sampling_dist;
    if (N > 1) {
        if (N == int(N))
            N -= 1;
        GEO::vec3 n_v1v2 = GEO::normalize(v2 - v1);
        for (int n = 1; n <= N; n++) {
//            ps.push_back(v1 + n_v1v2 * sampling_dist * n);
            if (tree.is_out_sf_envelope(v1 + n_v1v2 * sampling_dist * n, eps_2, prev_facet, sq_dist, nearest_point))
                return true;
        }
    }

    N = sqrt(ls[mod3(max_i + 2)]) / sampling_dist;
    if (N > 1) {
        if (N == int(N))
            N -= 1;
        GEO::vec3 n_v2v0 = GEO::normalize(v0 - v2);
        for (int n = 1; n <= N; n++) {
//            ps.push_back(v2 + n_v2v0 * sampling_dist * n);
            if (tree.is_out_sf_envelope(v2 + n_v2v0 * sampling_dist * n, eps_2, prev_facet, sq_dist, nearest_point))
                return true;
        }
    }

    return false;
}

void floatTetWild::get_new_tet_slots(Mesh& mesh, int n, std::vector<int>& new_conn_tets) {
    int cnt = 0;
    for (int i = mesh.t_empty_start; i < mesh.tets.size(); i++) {
        if (mesh.tets[i].is_removed) {
            new_conn_tets.push_back(i);
            cnt++;
            if (cnt == n) {
                mesh.t_empty_start = i + 1;
                break;
            }
        }
    }
    if (cnt < n) {
        for (int i = 0; i < n - cnt; i++)
            new_conn_tets.push_back(mesh.tets.size() + i);
        mesh.tets.resize(mesh.tets.size() + n - cnt);
        mesh.t_empty_start = mesh.tets.size();
    }
}

void floatTetWild::set_intersection(const std::unordered_set<int>& s1, const std::unordered_set<int>& s2, std::vector<int>& v) {
    if (s2.size() < s1.size()) {
        set_intersection(s2, s1, v);
        return;
    }
    v.clear();
    v.reserve(std::min(s1.size(), s2.size()));
    for (int x : s1) {
        if (s2.count(x)) {
            v.push_back(x);
        }
    }
//    std::sort(v.begin(), v.end());
}

void floatTetWild::set_intersection(const std::unordered_set<int>& s1, const std::unordered_set<int>& s2,
        std::unordered_set<int>& v) {
    if (s2.size() < s1.size()) {
        set_intersection(s2, s1, v);
        return;
    }
    v.clear();
    v.reserve(std::min(s1.size(), s2.size()));
    for (int x : s1) {
        if (s2.count(x)) {
            v.insert(x);
        }
    }
}

void floatTetWild::set_intersection(const std::unordered_set<int>& s1, const std::unordered_set<int>& s2, const std::unordered_set<int>& s3,
                                    std::vector<int>& v) {
    if (s2.size() < s1.size() && s2.size() < s1.size()) {
        set_intersection(s2, s1, s3, v);
        return;
    }

    if (s3.size() < s1.size() && s3.size() < s2.size()) {
        set_intersection(s3, s1, s2, v);
        return;
    }

    assert(s1.size() <= s2.size());
    assert(s1.size() <= s3.size());

    v.clear();
    v.reserve(s1.size());
    for (int x : s1) {
        if (s2.count(x) && s3.count(x)) {
            v.push_back(x);
            if(v.size() == 2)
                break;
        }
    }
}

void floatTetWild::set_intersection(const std::vector<int>& s11, const std::vector<int>& s22, std::vector<int>& v) {
    std::vector<int> s1 = s11;
    std::vector<int> s2 = s22;
    std::sort(s1.begin(), s1.end());
    std::sort(s2.begin(), s2.end());
    std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), std::back_inserter(v));
}

void floatTetWild::set_intersection(const std::vector<int>& s11, const std::vector<int>& s22, const std::vector<int>& s33, std::vector<int>& v) {
    std::vector<int> s1 = s11;
    std::vector<int> s2 = s22;
    std::vector<int> s3 = s33;
    std::sort(s1.begin(), s1.end());
    std::sort(s2.begin(), s2.end());
    std::sort(s3.begin(), s3.end());
    std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), std::back_inserter(v));
    auto it = std::set_intersection(v.begin(), v.end(), s3.begin(), s3.end(), v.begin());
    v.resize(it - v.begin());
}

void floatTetWild::set_intersection_sorted(const std::vector<int>& s1, const std::vector<int>& s2, const std::vector<int>& s3, std::vector<int>& v) {
    std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), std::back_inserter(v));
    auto it = std::set_intersection(v.begin(), v.end(), s3.begin(), s3.end(), v.begin());
    v.resize(it - v.begin());
}

void floatTetWild::pausee(std::string msg) {
    return;
    // if (!msg.empty())
    //     cout << msg << endl;
    // cout << "Is pausing... (Enter '0' to exit and other characters to continue.)" << endl;
    // char c = ' ';
    // std::cin >> c;
    // if (c == '0')
    //     exit(0);
}

bool floatTetWild::is_energy_unstable(const std::array<Scalar, 12>& T, Scalar res) {
    static const std::vector<std::array<int, 4>> combs = {{{0, 1, 3, 2}},
                                                          {{0, 2, 1, 3}},
                                                          {{0, 2, 3, 1}},
                                                          {{0, 3, 1, 2}},
                                                          {{0, 3, 2, 1}},
                                                          {{1, 0, 2, 3}},
                                                          {{1, 0, 3, 2}},
                                                          {{1, 2, 0, 3}},
                                                          {{1, 2, 3, 0}},
                                                          {{1, 3, 0, 2}},
                                                          {{1, 3, 2, 0}},
                                                          {{2, 0, 1, 3}},
                                                          {{2, 0, 3, 1}},
                                                          {{2, 1, 0, 3}},
                                                          {{2, 1, 3, 0}},
                                                          {{2, 3, 0, 1}},
                                                          {{2, 3, 1, 0}},
                                                          {{3, 0, 1, 2}},
                                                          {{3, 0, 2, 1}},
                                                          {{3, 1, 0, 2}},
                                                          {{3, 1, 2, 0}},
                                                          {{3, 2, 0, 1}},
                                                          {{3, 2, 1, 0}}};
    Scalar res0;
    if (std::isinf(res))
        return true;

    for (int i = 0; i < combs.size(); i++) {
        std::array<Scalar, 12> tmp_T;
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 3; k++)
                tmp_T[j * 3 + k] = T[combs[i][j] * 3 + k];
        }
        Scalar res1 = AMIPS_energy_aux(tmp_T);
        if (std::isinf(res1))
            continue;
        if (res0 == 0)
            res0 = res1;
//        if (res1 - res0 > 10)
        if (abs(res1 - res0)/res0 > 0.01)
            return true;
    }
    return false;
}

int cnt_stable = 0;
int cnt_large = 0;
#include <floattetwild/Rational.h>
Scalar floatTetWild::AMIPS_energy(const std::array<Scalar, 12>& T) {
    Scalar res = AMIPS_energy_aux(T);
    if(use_old_energy) {
        return res;
    }

    if (res > 1e8) {
//        //fortest
//        if (res > 1e10) {
//            cout << std::setprecision(16) << res << endl;
//            for (int i = 0; i < T.size(); i++) {
//                if (i % 3 == 0)
//                    cout << endl;
//                cout << T[i] << ", ";
//            }
//            cout << endl;
//            char c;
//            cin >> c;
//        }
//        //fortest

//        //fortest
//        cnt_large++;
//        if(!is_energy_unstable(T, res)){
//            cout<<(cnt_stable++)<<"/"<<cnt_large<<endl;
//        }
//        //fortest

        if(is_degenerate(Vector3(T[0], T[1], T[2]), Vector3(T[3], T[4], T[5]), Vector3(T[6], T[7], T[8]),
                    Vector3(T[9], T[10], T[11]))) {
            pausee("energy computation degenerate found!!!");
            return std::numeric_limits<double>::infinity();
        }

        std::array<triwild::Rational, 12> r_T;
        for (int j = 0; j < 12; j++)
            r_T[j] = T[j];
        const triwild::Rational twothird = triwild::Rational(2) / triwild::Rational(3);
        triwild::Rational tmp = ((-r_T[1 + 2] + r_T[1 + 5]) * r_T[1 + 1] + r_T[1 + 2] * r_T[1 + 7] +
                                 (r_T[1 + -1] - r_T[1 + 5]) * r_T[1 + 4] - r_T[1 + -1] * r_T[1 + 7]) * r_T[1 + 9] +
                                ((r_T[1 + 2] - r_T[1 + 5]) * r_T[1 + 0] - r_T[1 + 2] * r_T[1 + 6] +
                                 (-r_T[1 + -1] + r_T[1 + 5]) * r_T[1 + 3] + r_T[1 + -1] * r_T[1 + 6]) * r_T[1 + 10] +
                                (-r_T[1 + 2] * r_T[1 + 7] + (-r_T[1 + 8] + r_T[1 + 5]) * r_T[1 + 4] +
                                 r_T[1 + 8] * r_T[1 + 7]) * r_T[1 + 0] +
                                (r_T[1 + 2] * r_T[1 + 6] + (r_T[1 + 8] - r_T[1 + 5]) * r_T[1 + 3] - r_T[1 + 8] * r_T[1 + 6]) *
                                r_T[1 + 1] + (r_T[1 + 3] * r_T[1 + 7] - r_T[1 + 4] * r_T[1 + 6]) * (r_T[1 + -1] - r_T[1 + 8]);
        if(tmp == 0)
            return std::numeric_limits<double>::infinity();

        auto res_r = triwild::Rational(27) / 16 *
                     pow(tmp, -2) * pow(r_T[1 + 9] * r_T[1 + 9] +
                                   (-twothird * r_T[1 + 0] - twothird * r_T[1 + 3] - twothird * r_T[1 + 6]) *
                                   r_T[1 + 9] + r_T[1 + 10] * r_T[1 + 10] +
                                   (-twothird * r_T[1 + 1] - twothird * r_T[1 + 4] - twothird * r_T[1 + 7]) *
                                   r_T[1 + 10] + r_T[1 + 0] * r_T[1 + 0] +
                                   (-twothird * r_T[1 + 3] - twothird * r_T[1 + 6]) * r_T[1 + 0] +
                                   r_T[1 + 1] * r_T[1 + 1] +
                                   (-twothird * r_T[1 + 4] - twothird * r_T[1 + 7]) * r_T[1 + 1] +
                                   r_T[1 + 2] * r_T[1 + 2] +
                                   (-twothird * r_T[1 + -1] - twothird * r_T[1 + 8] - twothird * r_T[1 + 5]) *
                                   r_T[1 + 2] + r_T[1 + 3] * r_T[1 + 3] - twothird * r_T[1 + 3] * r_T[1 + 6] +
                                   r_T[1 + 4] * r_T[1 + 4] - twothird * r_T[1 + 4] * r_T[1 + 7] +
                                   r_T[1 + 5] * r_T[1 + 5] +
                                   (-twothird * r_T[1 + -1] - twothird * r_T[1 + 8]) * r_T[1 + 5] -
                                   twothird * r_T[1 + -1] * r_T[1 + 8] + r_T[1 + -1] * r_T[1 + -1] +
                                   r_T[1 + 8] * r_T[1 + 8] + r_T[1 + 6] * r_T[1 + 6] + r_T[1 + 7] * r_T[1 + 7], 3);
        return std::cbrt(res_r.to_double());
    } else {
        return res;
    }
}

Scalar floatTetWild::AMIPS_energy_aux(const std::array<Scalar, 12>& T) {
    Scalar helper_0[12];
    helper_0[0] = T[0];
    helper_0[1] = T[1];
    helper_0[2] = T[2];
    helper_0[3] = T[3];
    helper_0[4] = T[4];
    helper_0[5] = T[5];
    helper_0[6] = T[6];
    helper_0[7] = T[7];
    helper_0[8] = T[8];
    helper_0[9] = T[9];
    helper_0[10] = T[10];
    helper_0[11] = T[11];
    Scalar helper_1 = helper_0[2];
    Scalar helper_2 = helper_0[11];
    Scalar helper_3 = helper_0[0];
    Scalar helper_4 = helper_0[3];
    Scalar helper_5 = helper_0[9];
    Scalar helper_6 = 0.577350269189626 * helper_3 - 1.15470053837925 * helper_4 + 0.577350269189626 * helper_5;
    Scalar helper_7 = helper_0[1];
    Scalar helper_8 = helper_0[4];
    Scalar helper_9 = helper_0[7];
    Scalar helper_10 = helper_0[10];
    Scalar helper_11 = 0.408248290463863 * helper_10 + 0.408248290463863 * helper_7 + 0.408248290463863 * helper_8 -
                       1.22474487139159 * helper_9;
    Scalar helper_12 = 0.577350269189626 * helper_10 + 0.577350269189626 * helper_7 - 1.15470053837925 * helper_8;
    Scalar helper_13 = helper_0[6];
    Scalar helper_14 = -1.22474487139159 * helper_13 + 0.408248290463863 * helper_3 + 0.408248290463863 * helper_4 +
                       0.408248290463863 * helper_5;
    Scalar helper_15 = helper_0[5];
    Scalar helper_16 = helper_0[8];
    Scalar helper_17 = 0.408248290463863 * helper_1 + 0.408248290463863 * helper_15 - 1.22474487139159 * helper_16 +
                       0.408248290463863 * helper_2;
    Scalar helper_18 = 0.577350269189626 * helper_1 - 1.15470053837925 * helper_15 + 0.577350269189626 * helper_2;
    Scalar helper_19 = 0.5 * helper_13 + 0.5 * helper_4;
    Scalar helper_20 = 0.5 * helper_8 + 0.5 * helper_9;
    Scalar helper_21 = 0.5 * helper_15 + 0.5 * helper_16;
    Scalar helper_22 = (helper_1 - helper_2) * (helper_11 * helper_6 - helper_12 * helper_14) -
                       (-helper_10 + helper_7) * (-helper_14 * helper_18 + helper_17 * helper_6) +
                       (helper_3 - helper_5) * (-helper_11 * helper_18 + helper_12 * helper_17);
    Scalar res = -(helper_1 * (-1.5 * helper_1 + 0.5 * helper_2 + helper_21) +
                   helper_10 * (-1.5 * helper_10 + helper_20 + 0.5 * helper_7) +
                   helper_13 * (-1.5 * helper_13 + 0.5 * helper_3 + 0.5 * helper_4 + 0.5 * helper_5) +
                   helper_15 * (0.5 * helper_1 - 1.5 * helper_15 + 0.5 * helper_16 + 0.5 * helper_2) +
                   helper_16 * (0.5 * helper_1 + 0.5 * helper_15 - 1.5 * helper_16 + 0.5 * helper_2) +
                   helper_2 * (0.5 * helper_1 - 1.5 * helper_2 + helper_21) +
                   helper_3 * (helper_19 - 1.5 * helper_3 + 0.5 * helper_5) +
                   helper_4 * (0.5 * helper_13 + 0.5 * helper_3 - 1.5 * helper_4 + 0.5 * helper_5) +
                   helper_5 * (helper_19 + 0.5 * helper_3 - 1.5 * helper_5) +
                   helper_7 * (0.5 * helper_10 + helper_20 - 1.5 * helper_7) +
                   helper_8 * (0.5 * helper_10 + 0.5 * helper_7 - 1.5 * helper_8 + 0.5 * helper_9) +
                   helper_9 * (0.5 * helper_10 + 0.5 * helper_7 + 0.5 * helper_8 - 1.5 * helper_9))
            / std::cbrt(helper_22*helper_22);
//                 * pow(pow((helper_1 - helper_2) * (helper_11 * helper_6 - helper_12 * helper_14) -
//                         (-helper_10 + helper_7) * (-helper_14 * helper_18 + helper_17 * helper_6) +
//                         (helper_3 - helper_5) * (-helper_11 * helper_18 + helper_12 * helper_17), 2),
//                     -0.333333333333333);
    return res;
}

void floatTetWild::AMIPS_jacobian(const std::array<Scalar, 12>& T, Vector3& result_0){
    Scalar helper_0[12];
    helper_0[0] = T[0];
    helper_0[1] = T[1];
    helper_0[2] = T[2];
    helper_0[3] = T[3];
    helper_0[4] = T[4];
    helper_0[5] = T[5];
    helper_0[6] = T[6];
    helper_0[7] = T[7];
    helper_0[8] = T[8];
    helper_0[9] = T[9];
    helper_0[10] = T[10];
    helper_0[11] = T[11];
    Scalar helper_1 = helper_0[1];
    Scalar helper_2 = helper_0[10];
    Scalar helper_3 = helper_1 - helper_2;
    Scalar helper_4 = helper_0[0];
    Scalar helper_5 = helper_0[3];
    Scalar helper_6 = helper_0[9];
    Scalar helper_7 = 0.577350269189626*helper_4 - 1.15470053837925*helper_5 + 0.577350269189626*helper_6;
    Scalar helper_8 = helper_0[2];
    Scalar helper_9 = 0.408248290463863*helper_8;
    Scalar helper_10 = helper_0[5];
    Scalar helper_11 = 0.408248290463863*helper_10;
    Scalar helper_12 = helper_0[8];
    Scalar helper_13 = 1.22474487139159*helper_12;
    Scalar helper_14 = helper_0[11];
    Scalar helper_15 = 0.408248290463863*helper_14;
    Scalar helper_16 = helper_11 - helper_13 + helper_15 + helper_9;
    Scalar helper_17 = 0.577350269189626*helper_8;
    Scalar helper_18 = 1.15470053837925*helper_10;
    Scalar helper_19 = 0.577350269189626*helper_14;
    Scalar helper_20 = helper_17 - helper_18 + helper_19;
    Scalar helper_21 = helper_0[6];
    Scalar helper_22 = -1.22474487139159*helper_21 + 0.408248290463863*helper_4 + 0.408248290463863*helper_5 + 0.408248290463863*helper_6;
    Scalar helper_23 = helper_16*helper_7 - helper_20*helper_22;
    Scalar helper_24 = -helper_14 + helper_8;
    Scalar helper_25 = 0.408248290463863*helper_1;
    Scalar helper_26 = helper_0[4];
    Scalar helper_27 = 0.408248290463863*helper_26;
    Scalar helper_28 = helper_0[7];
    Scalar helper_29 = 1.22474487139159*helper_28;
    Scalar helper_30 = 0.408248290463863*helper_2;
    Scalar helper_31 = helper_25 + helper_27 - helper_29 + helper_30;
    Scalar helper_32 = helper_31*helper_7;
    Scalar helper_33 = 0.577350269189626*helper_1;
    Scalar helper_34 = 1.15470053837925*helper_26;
    Scalar helper_35 = 0.577350269189626*helper_2;
    Scalar helper_36 = helper_33 - helper_34 + helper_35;
    Scalar helper_37 = helper_22*helper_36;
    Scalar helper_38 = helper_4 - helper_6;
    Scalar helper_39 = helper_23*helper_3 - helper_24*(helper_32 - helper_37) - helper_38*(helper_16*helper_36 - helper_20*helper_31);
    Scalar helper_40 = pow(pow(helper_39, 2), -0.333333333333333);
    Scalar helper_41 = 0.707106781186548*helper_10 - 0.707106781186548*helper_12;
    Scalar helper_42 = 0.707106781186548*helper_26 - 0.707106781186548*helper_28;
    Scalar helper_43 = 0.5*helper_21 + 0.5*helper_5;
    Scalar helper_44 = 0.5*helper_26 + 0.5*helper_28;
    Scalar helper_45 = 0.5*helper_10 + 0.5*helper_12;
    Scalar helper_46 = 0.666666666666667*(helper_1*(-1.5*helper_1 + 0.5*helper_2 + helper_44) + helper_10*(-1.5*helper_10 + 0.5*helper_12 + 0.5*helper_14 + 0.5*helper_8) + helper_12*(0.5*helper_10 - 1.5*helper_12 + 0.5*helper_14 + 0.5*helper_8) + helper_14*(-1.5*helper_14 + helper_45 + 0.5*helper_8) + helper_2*(0.5*helper_1 - 1.5*helper_2 + helper_44) + helper_21*(-1.5*helper_21 + 0.5*helper_4 + 0.5*helper_5 + 0.5*helper_6) + helper_26*(0.5*helper_1 + 0.5*helper_2 - 1.5*helper_26 + 0.5*helper_28) + helper_28*(0.5*helper_1 + 0.5*helper_2 + 0.5*helper_26 - 1.5*helper_28) + helper_4*(-1.5*helper_4 + helper_43 + 0.5*helper_6) + helper_5*(0.5*helper_21 + 0.5*helper_4 - 1.5*helper_5 + 0.5*helper_6) + helper_6*(0.5*helper_4 + helper_43 - 1.5*helper_6) + helper_8*(0.5*helper_14 + helper_45 - 1.5*helper_8))/helper_39;
    Scalar helper_47 = -0.707106781186548*helper_21 + 0.707106781186548*helper_5;
    result_0[0] = -helper_40*(1.0*helper_21 - 3.0*helper_4 + helper_46*(helper_41*(-helper_1 + helper_2) - helper_42*(helper_14 - helper_8) - (-helper_17 + helper_18 - helper_19)*(-helper_25 - helper_27 + helper_29 - helper_30) + (-helper_33 + helper_34 - helper_35)*(-helper_11 + helper_13 - helper_15 - helper_9)) + 1.0*helper_5 + 1.0*helper_6);
    result_0[1] = helper_40*(3.0*helper_1 - 1.0*helper_2 - 1.0*helper_26 - 1.0*helper_28 + helper_46*(helper_23 + helper_24*helper_47 - helper_38*helper_41));
    result_0[2] = helper_40*(-1.0*helper_10 - 1.0*helper_12 - 1.0*helper_14 + helper_46*(-helper_3*helper_47 - helper_32 + helper_37 + helper_38*helper_42) + 3.0*helper_8);
}

void floatTetWild::AMIPS_hessian(const std::array<Scalar, 12>& T, Matrix3& result_0){
    Scalar helper_0[12];
    helper_0[0] = T[0];
    helper_0[1] = T[1];
    helper_0[2] = T[2];
    helper_0[3] = T[3];
    helper_0[4] = T[4];
    helper_0[5] = T[5];
    helper_0[6] = T[6];
    helper_0[7] = T[7];
    helper_0[8] = T[8];
    helper_0[9] = T[9];
    helper_0[10] = T[10];
    helper_0[11] = T[11];
    Scalar helper_1 = helper_0[2];
    Scalar helper_2 = helper_0[11];
    Scalar helper_3 = helper_1 - helper_2;
    Scalar helper_4 = helper_0[0];
    Scalar helper_5 = 0.577350269189626*helper_4;
    Scalar helper_6 = helper_0[3];
    Scalar helper_7 = 1.15470053837925*helper_6;
    Scalar helper_8 = helper_0[9];
    Scalar helper_9 = 0.577350269189626*helper_8;
    Scalar helper_10 = helper_5 - helper_7 + helper_9;
    Scalar helper_11 = helper_0[1];
    Scalar helper_12 = 0.408248290463863*helper_11;
    Scalar helper_13 = helper_0[4];
    Scalar helper_14 = 0.408248290463863*helper_13;
    Scalar helper_15 = helper_0[7];
    Scalar helper_16 = 1.22474487139159*helper_15;
    Scalar helper_17 = helper_0[10];
    Scalar helper_18 = 0.408248290463863*helper_17;
    Scalar helper_19 = helper_12 + helper_14 - helper_16 + helper_18;
    Scalar helper_20 = helper_10*helper_19;
    Scalar helper_21 = 0.577350269189626*helper_11;
    Scalar helper_22 = 1.15470053837925*helper_13;
    Scalar helper_23 = 0.577350269189626*helper_17;
    Scalar helper_24 = helper_21 - helper_22 + helper_23;
    Scalar helper_25 = 0.408248290463863*helper_4;
    Scalar helper_26 = 0.408248290463863*helper_6;
    Scalar helper_27 = helper_0[6];
    Scalar helper_28 = 1.22474487139159*helper_27;
    Scalar helper_29 = 0.408248290463863*helper_8;
    Scalar helper_30 = helper_25 + helper_26 - helper_28 + helper_29;
    Scalar helper_31 = helper_24*helper_30;
    Scalar helper_32 = helper_3*(helper_20 - helper_31);
    Scalar helper_33 = helper_4 - helper_8;
    Scalar helper_34 = 0.408248290463863*helper_1;
    Scalar helper_35 = helper_0[5];
    Scalar helper_36 = 0.408248290463863*helper_35;
    Scalar helper_37 = helper_0[8];
    Scalar helper_38 = 1.22474487139159*helper_37;
    Scalar helper_39 = 0.408248290463863*helper_2;
    Scalar helper_40 = helper_34 + helper_36 - helper_38 + helper_39;
    Scalar helper_41 = helper_24*helper_40;
    Scalar helper_42 = 0.577350269189626*helper_1;
    Scalar helper_43 = 1.15470053837925*helper_35;
    Scalar helper_44 = 0.577350269189626*helper_2;
    Scalar helper_45 = helper_42 - helper_43 + helper_44;
    Scalar helper_46 = helper_19*helper_45;
    Scalar helper_47 = helper_41 - helper_46;
    Scalar helper_48 = helper_33*helper_47;
    Scalar helper_49 = helper_11 - helper_17;
    Scalar helper_50 = helper_10*helper_40;
    Scalar helper_51 = helper_30*helper_45;
    Scalar helper_52 = helper_50 - helper_51;
    Scalar helper_53 = helper_49*helper_52;
    Scalar helper_54 = helper_32 + helper_48 - helper_53;
    Scalar helper_55 = pow(helper_54, 2);
    Scalar helper_56 = pow(helper_55, -0.333333333333333);
    Scalar helper_57 = 1.0*helper_27 - 3.0*helper_4 + 1.0*helper_6 + 1.0*helper_8;
    Scalar helper_58 = 0.707106781186548*helper_13;
    Scalar helper_59 = 0.707106781186548*helper_15;
    Scalar helper_60 = helper_58 - helper_59;
    Scalar helper_61 = helper_3*helper_60;
    Scalar helper_62 = 0.707106781186548*helper_35 - 0.707106781186548*helper_37;
    Scalar helper_63 = helper_49*helper_62;
    Scalar helper_64 = helper_47 + helper_61 - helper_63;
    Scalar helper_65 = 1.33333333333333/helper_54;
    Scalar helper_66 = 1.0/helper_55;
    Scalar helper_67 = 0.5*helper_27 + 0.5*helper_6;
    Scalar helper_68 = -1.5*helper_4 + helper_67 + 0.5*helper_8;
    Scalar helper_69 = 0.5*helper_4 + helper_67 - 1.5*helper_8;
    Scalar helper_70 = -1.5*helper_27 + 0.5*helper_4 + 0.5*helper_6 + 0.5*helper_8;
    Scalar helper_71 = 0.5*helper_27 + 0.5*helper_4 - 1.5*helper_6 + 0.5*helper_8;
    Scalar helper_72 = 0.5*helper_13 + 0.5*helper_15;
    Scalar helper_73 = -1.5*helper_11 + 0.5*helper_17 + helper_72;
    Scalar helper_74 = 0.5*helper_11 - 1.5*helper_17 + helper_72;
    Scalar helper_75 = 0.5*helper_11 + 0.5*helper_13 - 1.5*helper_15 + 0.5*helper_17;
    Scalar helper_76 = 0.5*helper_11 - 1.5*helper_13 + 0.5*helper_15 + 0.5*helper_17;
    Scalar helper_77 = 0.5*helper_35 + 0.5*helper_37;
    Scalar helper_78 = -1.5*helper_1 + 0.5*helper_2 + helper_77;
    Scalar helper_79 = 0.5*helper_1 - 1.5*helper_2 + helper_77;
    Scalar helper_80 = 0.5*helper_1 + 0.5*helper_2 + 0.5*helper_35 - 1.5*helper_37;
    Scalar helper_81 = 0.5*helper_1 + 0.5*helper_2 - 1.5*helper_35 + 0.5*helper_37;
    Scalar helper_82 = helper_1*helper_78 + helper_11*helper_73 + helper_13*helper_76 + helper_15*helper_75 + helper_17*helper_74 + helper_2*helper_79 + helper_27*helper_70 + helper_35*helper_81 + helper_37*helper_80 + helper_4*helper_68 + helper_6*helper_71 + helper_69*helper_8;
    Scalar helper_83 = 0.444444444444444*helper_66*helper_82;
    Scalar helper_84 = helper_66*helper_82;
    Scalar helper_85 = -helper_32 - helper_48 + helper_53;
    Scalar helper_86 = 1.0/helper_85;
    Scalar helper_87 = helper_86*pow(pow(helper_85, 2), -0.333333333333333);
    Scalar helper_88 = 0.707106781186548*helper_6;
    Scalar helper_89 = 0.707106781186548*helper_27;
    Scalar helper_90 = helper_88 - helper_89;
    Scalar helper_91 = 0.666666666666667*helper_10*helper_40 + 0.666666666666667*helper_3*helper_90 - 0.666666666666667*helper_30*helper_45 - 0.666666666666667*helper_33*helper_62;
    Scalar helper_92 = -3.0*helper_11 + 1.0*helper_13 + 1.0*helper_15 + 1.0*helper_17;
    Scalar helper_93 = -helper_11 + helper_17;
    Scalar helper_94 = -helper_1 + helper_2;
    Scalar helper_95 = -helper_21 + helper_22 - helper_23;
    Scalar helper_96 = -helper_34 - helper_36 + helper_38 - helper_39;
    Scalar helper_97 = -helper_42 + helper_43 - helper_44;
    Scalar helper_98 = -helper_12 - helper_14 + helper_16 - helper_18;
    Scalar helper_99 = -0.666666666666667*helper_60*helper_94 + 0.666666666666667*helper_62*helper_93 + 0.666666666666667*helper_95*helper_96 - 0.666666666666667*helper_97*helper_98;
    Scalar helper_100 = helper_3*helper_90;
    Scalar helper_101 = helper_33*helper_62;
    Scalar helper_102 = helper_100 - helper_101 + helper_52;
    Scalar helper_103 = -helper_60*helper_94 + helper_62*helper_93 + helper_95*helper_96 - helper_97*helper_98;
    Scalar helper_104 = 0.444444444444444*helper_102*helper_103*helper_82*helper_86 + helper_57*helper_91 - helper_92*helper_99;
    Scalar helper_105 = 1.85037170770859e-17*helper_1*helper_78 + 1.85037170770859e-17*helper_11*helper_73 + 1.85037170770859e-17*helper_13*helper_76 + 1.85037170770859e-17*helper_15*helper_75 + 1.85037170770859e-17*helper_17*helper_74 + 1.85037170770859e-17*helper_2*helper_79 + 1.85037170770859e-17*helper_27*helper_70 + 1.85037170770859e-17*helper_35*helper_81 + 1.85037170770859e-17*helper_37*helper_80 + 1.85037170770859e-17*helper_4*helper_68 + 1.85037170770859e-17*helper_6*helper_71 + 1.85037170770859e-17*helper_69*helper_8;
    Scalar helper_106 = helper_64*helper_82*helper_86;
    Scalar helper_107 = -0.666666666666667*helper_10*helper_19 + 0.666666666666667*helper_24*helper_30 + 0.666666666666667*helper_33*helper_60 - 0.666666666666667*helper_49*helper_90;
    Scalar helper_108 = -3.0*helper_1 + 1.0*helper_2 + 1.0*helper_35 + 1.0*helper_37;
    Scalar helper_109 = -helper_20 + helper_31 + helper_33*helper_60 - helper_49*helper_90;
    Scalar helper_110 = 0.444444444444444*helper_109*helper_82*helper_86;
    Scalar helper_111 = helper_103*helper_110 + helper_107*helper_57 - helper_108*helper_99;
    Scalar helper_112 = -helper_4 + helper_8;
    Scalar helper_113 = -helper_88 + helper_89;
    Scalar helper_114 = -helper_5 + helper_7 - helper_9;
    Scalar helper_115 = -helper_25 - helper_26 + helper_28 - helper_29;
    Scalar helper_116 = helper_82*helper_86*(helper_112*helper_62 + helper_113*helper_94 + helper_114*helper_96 - helper_115*helper_97);
    Scalar helper_117 = -helper_100 + helper_101 - helper_50 + helper_51;
    Scalar helper_118 = -helper_102*helper_110 + helper_107*helper_92 + helper_108*helper_91;
    Scalar helper_119 = helper_82*helper_86*(helper_112*(-helper_58 + helper_59) - helper_113*helper_93 - helper_114*helper_98 + helper_115*helper_95);
    result_0(0, 0) = helper_56*(helper_57*helper_64*helper_65 - pow(helper_64, 2)*helper_83 + 0.666666666666667*helper_64*helper_84*(-helper_41 + helper_46 - helper_61 + helper_63) + 3.0);
    result_0(0, 1) = helper_87*(helper_104 - helper_105*helper_35 + helper_106*helper_91);
    result_0(0, 2) = helper_87*(helper_106*helper_107 + helper_111);
    result_0(1, 0) = helper_87*(helper_104 + helper_116*helper_99);
    result_0(1, 1) = helper_56*(-pow(helper_117, 2)*helper_83 + helper_117*helper_65*helper_92 + helper_117*helper_84*helper_91 + 3.0);
    result_0(1, 2) = helper_87*(-helper_105*helper_6 - helper_107*helper_116 + helper_118);
    result_0(2, 0) = helper_87*(-helper_105*helper_13 + helper_111 + helper_119*helper_99);
    result_0(2, 1) = helper_87*(helper_118 - helper_119*helper_91);
    result_0(2, 2) = helper_56*(-helper_108*helper_109*helper_65 - 1.11111111111111*pow(helper_109, 2)*helper_84 + 3.0);
}
