//
// Created by Yixin Hu on 2019-08-27.
//

#ifndef FLOATTETWILD_TRIANGLEINSERTION_H
#define FLOATTETWILD_TRIANGLEINSERTION_H

#include <floattetwild/Mesh.hpp>
#include <floattetwild/AABBWrapper.h>

#ifdef FLOAT_TETWILD_USE_TBB
#include <tbb/concurrent_vector.h>
#endif

namespace floatTetWild {
    void insert_triangles(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
                          const std::vector<int> &input_tags,
                          Mesh &mesh, std::vector<bool> &is_face_inserted, AABBWrapper &tree,
                          bool is_again);

    void insert_one_triangle(int f_id, const std::vector<Vector3> &input_vertices,
                             const std::vector<Vector3i> &input_faces, const std::vector<int> &input_tags,
                             Mesh &mesh, std::vector<bool> &is_face_inserted, AABBWrapper &tree,
                             bool is_again);

    class CutMesh {
    public:
        std::vector<int> v_ids;
        std::map<int, int> map_v_ids;
        std::vector<Scalar> to_plane_dists;
        std::vector<bool> is_snapped;

        std::vector<std::array<int, 4>> tets;
        std::vector<std::array<int, 4>> opp_t_ids = {{-1, -1, -1, -1}};

        const Mesh &mesh;
        const Vector3 &p_n;
        const std::array<Vector3, 3> &p_vs;

        CutMesh(const Mesh &_mesh, const Vector3 &_p_n, const std::array<Vector3, 3> &_p_vs) :
                mesh(_mesh), p_n(_p_n), p_vs(_p_vs) {}

        void construct(const std::vector<int>& cut_t_ids);
        bool snap_to_plane();
        void expand(std::vector<int>& cut_t_ids);

        inline bool is_v_on_plane(int lv_id){
            if(is_snapped[lv_id] || to_plane_dists[lv_id])
                return true;
            return false;
        }
    };
}


#endif //FLOATTETWILD_TRIANGLEINSERTION_H
