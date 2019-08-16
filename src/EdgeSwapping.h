#ifndef FLOATTETWILD_EDGESWAPPING_H
#define FLOATTETWILD_EDGESWAPPING_H

#include <floattetwild/Mesh.hpp>

namespace floatTetWild {
    extern bool is_es_check;

    void edge_swapping(Mesh& mesh);
    bool remove_an_edge_32(Mesh& mesh, int v1_id, int v2_id, const std::vector<int>& n12_t_ids, std::vector<std::array<int, 2>>& new_edges);
    bool remove_an_edge_44(Mesh& mesh, int v1_id, int v2_id, const std::vector<int>& n12_t_ids, std::vector<std::array<int, 2>>& new_edges);
    bool remove_an_edge_56(Mesh& mesh, int v1_id, int v2_id, const std::vector<int>& n12_t_ids, std::vector<std::array<int, 2>>& new_edges);
}

#endif //FLOATTETWILD_EDGESWAPPING_H
