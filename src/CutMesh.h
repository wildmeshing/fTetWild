//
// Created by Yixin Hu on 9/12/19.
//

#ifndef FLOATTETWILD_CUTMESH_H
#define FLOATTETWILD_CUTMESH_H

#include <floattetwild/Mesh.hpp>
#include <floattetwild/AABBWrapper.h>

namespace floatTetWild {
    class CutMesh {
    public:
        std::vector<int> v_ids;
        std::map<int, int> map_v_ids;
        std::vector<std::array<int, 4>> tets;
        std::vector<std::array<int, 4>> opp_t_ids;

        std::vector<Scalar> to_plane_dists;
        std::vector<bool> is_snapped;

        const Mesh &mesh;
        const Vector3 &p_n;
        const std::array<Vector3, 3> &p_vs;

        CutMesh(const Mesh &_mesh, const Vector3 &_p_n, const std::array<Vector3, 3> &_p_vs) :
                mesh(_mesh), p_n(_p_n), p_vs(_p_vs) {}

        void construct(const std::vector<int> &cut_t_ids);

        bool snap_to_plane();

        void expand(std::vector<int> &cut_t_ids);

        bool get_intersecting_edges_and_points(std::vector<Vector3> &points,
                                               std::map<std::array<int, 2>, int> &map_edge_to_intersecting_point,
                                               std::vector<int> &subdivide_t_ids);

        void get_one_ring_t_ids(std::vector<int> &old_t_ids, std::vector<int> &neighbor_t_ids);

        inline bool is_v_on_plane(int lv_id) {
            if (is_snapped[lv_id] || to_plane_dists[lv_id] == 0)
                return true;
            return false;
        }

        inline Scalar get_to_plane_dist(const Vector3 &p) {
            return p_n.dot(p - p_vs[0]);
        }
    };
}


#endif //FLOATTETWILD_CUTMESH_H
