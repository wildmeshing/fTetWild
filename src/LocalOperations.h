// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#ifndef FLOATTETWILD_LOCALOPERATIONS_H
#define FLOATTETWILD_LOCALOPERATIONS_H

#include <floattetwild/Mesh.hpp>
#include <floattetwild/AABBWrapper.h>

namespace floatTetWild {
    extern bool use_old_energy;
    extern std::string envelope_log_csv;
    extern int envelope_log_csv_cnt;

    // void init_b_tree(const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces, GEO::Mesh& b_mesh);

    int get_opp_t_id(const Mesh& mesh, int t_id, int j);
    void set_opp_t_id(Mesh& mesh, int t_id, int j);
    inline int get_local_f_id(int t_id, int v1_id, int v2_id, int v3_id, Mesh &mesh) {
        for (int j = 0; j < 4; j++) {
            if (mesh.tets[t_id][j] != v1_id && mesh.tets[t_id][j] != v2_id && mesh.tets[t_id][j] != v3_id)
                return j;
        }
        assert(false);
        return -1;
    }

    inline GEO::vec3 to_geo_p(const Vector3& p){
        return GEO::vec3(p[0], p[1], p[2]);
    }

    void get_all_edges(const Mesh& mesh, std::vector<std::array<int, 2>>& edges);
    void get_all_edges(const Mesh& mesh, const std::vector<int>& t_ids, std::vector<std::array<int, 2>>& edges, bool skip_freezed = false);

    Scalar get_edge_length(const Mesh& mesh, int v1_id, int v2_id);
    Scalar get_edge_length_2(const Mesh& mesh, int v1_id, int v2_id);

    Scalar get_quality(const Mesh& mesh, const MeshTet& t);
    Scalar get_quality(const Mesh& mesh, int t_id);
    Scalar get_quality(const MeshVertex& v0, const MeshVertex& v1, const MeshVertex& v2, const MeshVertex& v3);
    Scalar get_quality(const Vector3& v0, const Vector3& v1, const Vector3& v2, const Vector3& v3);
    void get_max_avg_energy(const Mesh& mesh, Scalar& max_energy, Scalar& avg_energy);
    Scalar get_mid_energy(const Mesh& mesh);

    bool is_inverted(const Mesh& mesh, int t_id);
    bool is_inverted(const Mesh& mesh, int t_id, int j, const Vector3& new_p);
    bool is_inverted(const MeshVertex& v0, const MeshVertex& v1, const MeshVertex& v2, const MeshVertex& v3);
    bool is_inverted(const Vector3& v0, const Vector3& v1, const Vector3& v2, const Vector3& v3);
    bool is_degenerate(const Vector3& v0, const Vector3& v1, const Vector3& v2, const Vector3& v3);

    bool is_out_envelope(Mesh& mesh, int v_id, const Vector3& new_pos, const AABBWrapper& tree);
    bool is_out_boundary_envelope(const Mesh& mesh, int v_id, const Vector3& new_pos, const AABBWrapper& tree);
    void sample_triangle(const std::array<Vector3, 3>& vs, std::vector<GEO::vec3>& ps, Scalar sampling_dist);
    bool sample_triangle_and_check_is_out(const std::array<Vector3, 3>& vs, Scalar sampling_dist,
            Scalar eps_2, const AABBWrapper& tree, GEO::index_t& prev_facet);

    bool is_bbox_edge(const Mesh& mesh, int v1_id, int v2_id, const std::vector<int>& n12_t_ids);
    bool is_surface_edge(const Mesh& mesh, int v1_id, int v2_id, const std::vector<int>& n12_t_ids);
    bool is_boundary_edge(const Mesh& mesh, int v1_id, int v2_id, const AABBWrapper& tree);
    bool is_valid_edge(const Mesh& mesh, int v1_id, int v2_id);
    bool is_valid_edge(const Mesh& mesh, int v1_id, int v2_id, const std::vector<int>& n12_t_ids);

    bool is_isolate_surface_point(const Mesh& mesh, int v_id);
    bool is_point_out_envelope(const Mesh& mesh, const Vector3& p, const AABBWrapper& tree);
    bool is_point_out_boundary_envelope(const Mesh& mesh, const Vector3& p, const AABBWrapper& tree);

    void get_new_tet_slots(Mesh& mesh, int n, std::vector<int>& new_conn_tets);

    inline Scalar get_area(const Vector3& a, const Vector3& b, const Vector3& c) {
        return ((b - c).cross(a - c)).norm();
    }

    template<typename T>
    void vector_unique(std::vector<T>& v){
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
    }
    template<typename T>
    bool vector_erase(std::vector<T>& v, const T& t){
        auto it = std::find(v.begin(), v.end(), t);
        if(it == v.end())
            return false;
        v.erase(it);
        return true;
    }
    template<typename T>
    void vector_print(const std::vector<T>& v, std::string s = " "){
        for(auto i: v)
            cout<<i<<s;
        cout<<endl;
    }
    template<typename T>
    bool vector_contains(const std::vector<T>& v, const T& t){
        if(v.empty())
            return false;
        if(std::find(v.begin(), v.end(), t)!=v.end())
            return true;
        return false;
    }
    void set_intersection(const std::unordered_set<int>& s1, const std::unordered_set<int>& s2, std::vector<int>& v);
    void set_intersection(const std::unordered_set<int>& s1, const std::unordered_set<int>& s2, std::unordered_set<int>& v);
    void set_intersection(const std::unordered_set<int>& s1, const std::unordered_set<int>& s2, const std::unordered_set<int>& s3, std::vector<int>& v);

    void set_intersection(const std::vector<int>& s1, const std::vector<int>& s2, std::vector<int>& v);
    void set_intersection(const std::vector<int>& s1, const std::vector<int>& s2, const std::vector<int>& s3, std::vector<int>& v);
    void set_intersection_sorted(const std::vector<int>& s1, const std::vector<int>& s2, const std::vector<int>& s3, std::vector<int>& v);

    inline int mod4(int j) {
//        assert(j >= 0 && j < 7);
//        static const std::array<int, 8> a = {{0, 1, 2, 3, 0, 1, 2, 3}};
//        assert(a[j] == j % 4);
//        return a[j];

//        return j >= 4 ? j % 4 : j;
        return j%4;
    }

    inline int mod3(int j) {
//        assert(j >= 0 && j < 5);
//        static const std::array<int, 6> a = {{0, 1, 2, 0, 1, 2}};
//        assert(a[j] == j % 3);
//        return a[j];
        return j%3;
    }

    inline int mod2(int j) {
//        assert(j >= 0 && j < 3);
//        static const std::array<int, 4> a = {{0, 1, 0, 1}};
//        assert(a[j] == j % 2);
//        return a[j];
        return j%2;
    }

    void pausee(std::string msg = "");

    ///////////////
    class ElementInQueue{
    public:
        std::array<int, 2> v_ids;
        Scalar weight;

        ElementInQueue(){}
        ElementInQueue(const std::array<int, 2>& ids, Scalar w): v_ids(ids), weight(w){}
    };
    struct cmp_l {
        bool operator()(const ElementInQueue &e1, const ElementInQueue &e2) {
            if (e1.weight == e2.weight)
                return e1.v_ids > e2.v_ids;
            return e1.weight < e2.weight;
        }
    };
    struct cmp_s {
        bool operator()(const ElementInQueue &e1, const ElementInQueue &e2) {
            if (e1.weight == e2.weight)
                return e1.v_ids < e2.v_ids;
            return e1.weight > e2.weight;
        }
    };

    Scalar AMIPS_energy_aux(const std::array<Scalar, 12>& T);
    bool is_energy_unstable(const std::array<Scalar, 12>& T, Scalar res);
    Scalar AMIPS_energy(const std::array<Scalar, 12>& T);
    void AMIPS_jacobian(const std::array<Scalar, 12>& T, Vector3& result_0);
    void AMIPS_hessian(const std::array<Scalar, 12>& T, Matrix3& result_0);
}


#endif //FLOATTETWILD_LOCALOPERATIONS_H
