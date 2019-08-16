#ifndef FLOATTETWILD_EDGESPLITTING_H
#define FLOATTETWILD_EDGESPLITTING_H

#include <floattetwild/Mesh.hpp>

namespace floatTetWild {
    void edge_splitting(Mesh& mesh);
    bool split_an_edge(Mesh& mesh, int v1_id, int v2_id, bool is_repush, std::vector<std::array<int, 2>>& new_edges);
}

#endif //FLOATTETWILD_EDGESPLITTING_H
