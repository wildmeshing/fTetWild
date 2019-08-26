#ifndef FLOATTETWILD_FLOATTETCUTTINGPARALLEL_H
#define FLOATTETWILD_FLOATTETCUTTINGPARALLEL_H

#include <floattetwild/Mesh.hpp>
#include <floattetwild/AABBWrapper.h>

#ifdef FLOAT_TETWILD_USE_TBB
#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include <tbb/atomic.h>
#include <floattetwild/FloatTetCuttingParallel.h>
#include <tbb/concurrent_vector.h>
#endif


namespace floatTetWild {
    void generate_coloring_graph(const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces,
                                 const Mesh& mesh, std::vector<std::array<int, 2>>& graph, bool is_again = false);

    void generate_coloring(const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces,
                             const Mesh& mesh, std::vector<int> &colors, bool is_again);

    void box_sets(const int threshold, const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces,
                             const Mesh& mesh, std::vector<std::vector<int>> &concurrent_sets, std::vector<int> &serial_set, bool is_again);

    void partition_mesh(const std::vector<std::vector<int>>& partition_t_ids,
                        const Mesh& mesh, std::vector<Mesh>& sub_meshes, //tbb::concurrent_vector<Mesh>& sub_meshes,
                        const std::vector<std::array<std::vector<int>, 4>>& cut_f_ids,
//                        tbb::concurrent_vector<std::vector<std::array<std::vector<int>, 4>>>& sub_cut_f_ids,
                        std::vector<std::vector<std::array<std::vector<int>, 4>>>& sub_cut_f_ids,
                        std::vector<std::array<int, 2>>& map_input_to_partition,
                        std::vector<std::vector<std::array<int, 2>>>& v_p_ids, bool is_again);

    void merge_meshes(const std::vector<Mesh>& sub_meshes, //const tbb::concurrent_vector<Mesh>& sub_meshes,
            Mesh& mesh,
            const std::vector<std::vector<std::array<std::vector<int>, 4>>>& sub_cut_f_ids,
//            const tbb::concurrent_vector<std::vector<std::array<std::vector<int>, 4>>>& sub_cut_f_ids,
            std::vector<std::array<std::vector<int>, 4>>& cut_f_ids,
            const std::vector<std::array<int, 2>>& map_input_to_partition, std::vector<int>& map_input_to_merged_mesh,
            const std::vector<std::vector<std::array<int, 2>>>& v_p_ids);

    bool is_cutting_cross_partitions(Mesh& sub_mesh, const std::vector<int>& intersection_results_wn);
}

#endif //FLOATTETWILD_FLOATTETCUTTINGPARALLEL_H
