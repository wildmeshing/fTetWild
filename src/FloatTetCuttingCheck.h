#ifndef FLOATTETWILD_FLOATTETCUTTINGCHECK_H
#define FLOATTETWILD_FLOATTETCUTTINGCHECK_H

#include <floattetwild/Mesh.hpp>
#include <floattetwild/AABBWrapper.h>

namespace floatTetWild {
#define COMMON_INPUT_FOR_CHECK input_vertices,input_faces,mesh,cut_f_ids

    void
    find_tets_for_cut_new(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces, int f_id,
                          Mesh &mesh, std::vector<int> &intersection_results, std::vector<int> &oris);

    bool check(const Mesh &mesh);
    void check_is_surface_fs(const Mesh &mesh);
    void check_cut_f_ids(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
                         const Mesh &mesh, std::vector<std::array<std::vector<int>, 4>>& cut_f_ids);

    void plot_cover_for_tetf(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
            const Mesh &mesh, const std::vector<std::array<std::vector<int>, 4>>& cut_f_ids,
            int t_id, int j);
    void plot_cover_for_trif(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
            const Mesh &mesh, const std::vector<std::array<std::vector<int>, 4>>& cut_f_ids,
            int f_id);

    void find_bad_cover_for_tetf(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
            const Mesh &mesh, const std::vector<std::array<std::vector<int>, 4>>& cut_f_ids, const AABBWrapper& tree,
            int t_id, int j, bool is_on_surface);
    void find_bad_cover_for_trif(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
            const Mesh &mesh, const std::vector<std::array<std::vector<int>, 4>>& cut_f_ids, int f_id = -1);

    void get_samples(const std::array<Vector2, 3>& rf, std::vector<Vector2>& cs, Scalar ratio = 0.1);

}

#endif //FLOATTETWILD_FLOATTETCUTTINGCHECK_H
