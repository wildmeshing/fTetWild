#ifndef FLOATTETWILD_EDGECOLLAPSING_H
#define FLOATTETWILD_EDGECOLLAPSING_H

#include <floattetwild/Mesh.hpp>
#include <floattetwild/AABBWrapper.h>

namespace floatTetWild {
    void edge_collapsing(Mesh& mesh, const AABBWrapper& tree);
    int collapse_an_edge(Mesh& mesh, int v1_id, int v2_id, const AABBWrapper& tree,
            std::vector<std::array<int, 2>>& new_edges, int ts, std::vector<int>& tet_tss);

    bool is_edge_freezed(Mesh& mesh, int v1_id, int v2_id);
    bool is_collapsable_bbox(Mesh& mesh, int v1_id, int v2_id);
    bool is_collapsable_length(Mesh& mesh, int v1_id, int v2_id, Scalar l_2);
    bool is_collapsable_boundary(Mesh& mesh, int v1_id, int v2_id);
}

#endif //FLOATTETWILD_EDGECOLLAPSING_H
