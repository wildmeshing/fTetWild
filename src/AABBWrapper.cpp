#include <floattetwild/LocalOperations.h>
#include <floattetwild/AABBWrapper.h>
#include <geogram/mesh/mesh_reorder.h>
#include <geogram/basic/geometry_nd.h>

void floatTetWild::AABBWrapper::init_b_mesh(const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces) {
    b_mesh.clear(false,false);
//    std::vector<std::array<int, 2>> edges;
//    for(int i=0;i<sf_mesh.facets.nb();i++){
//        for(int j=0;j<3;j++) {
//            if(sf_mesh.facets.adjacent(i, j)==GEO::NO_FACET){
//                edges.push_back({{}})
//            }
//        }
//    }

    std::vector<std::vector<int>> conn_tris(input_vertices.size());
    std::vector<std::array<int, 2>> all_edges;
    all_edges.reserve(input_faces.size() * 3);
    for (int i = 0; i < input_faces.size(); i++) {
        for (int j = 0; j < 3; j++) {
            conn_tris[input_faces[i][j]].push_back(i);
            if (input_faces[i][j] < input_faces[i][(j + 1) % 3])
                all_edges.push_back({{input_faces[i][j], input_faces[i][(j + 1) % 3]}});
            else
                all_edges.push_back({{input_faces[i][(j + 1) % 3], input_faces[i][j]}});
        }
    }
    vector_unique(all_edges);

    std::vector<std::array<int, 2>> b_edges;
    for (auto &e:all_edges) {
        std::vector<int> tmp;
        std::set_intersection(conn_tris[e[0]].begin(), conn_tris[e[0]].end(),
                              conn_tris[e[1]].begin(), conn_tris[e[1]].end(), std::back_inserter(tmp));
        if (tmp.size() == 1) {
            b_edges.push_back(e);
        }
    }

    if (b_edges.empty()) {
        b_mesh.vertices.clear();
        b_mesh.vertices.create_vertices(1);
        b_mesh.vertices.point(0) = GEO::vec3(0, 0, 0);
        b_mesh.facets.clear();
        b_mesh.facets.create_triangles(1);
        b_mesh.facets.set_vertex(0, 0, 0);
        b_mesh.facets.set_vertex(0, 1, 0);
        b_mesh.facets.set_vertex(0, 2, 0);
    } else {
        b_mesh.vertices.clear();
        b_mesh.vertices.create_vertices((int) b_edges.size() * 2);
        int cnt = 0;
        for (auto &e:b_edges) {
            for (int j = 0; j < 2; j++) {
                GEO::vec3 &p = b_mesh.vertices.point(cnt++);
                p[0] = input_vertices[e[j]][0];
                p[1] = input_vertices[e[j]][1];
                p[2] = input_vertices[e[j]][2];
            }
        }
        b_mesh.facets.clear();
        b_mesh.facets.create_triangles((int) b_edges.size());
        for (int i = 0; i < b_edges.size(); i++) {
            b_mesh.facets.set_vertex(i, 0, i * 2);
            b_mesh.facets.set_vertex(i, 1, i * 2);
            b_mesh.facets.set_vertex(i, 2, i * 2 + 1);
        }
    }

    mesh_reorder(b_mesh, GEO::MESH_ORDER_MORTON);
}

void floatTetWild::AABBWrapper::init_tmp_b_mesh_and_tree(const Mesh& mesh, const std::vector<std::array<int, 2>>& b_edges){
    if (b_edges.empty()) {
        tmp_b_mesh.vertices.clear();
        tmp_b_mesh.vertices.create_vertices(1);
        tmp_b_mesh.vertices.point(0) = GEO::vec3(0, 0, 0);
        tmp_b_mesh.facets.clear();
        tmp_b_mesh.facets.create_triangles(1);
        tmp_b_mesh.facets.set_vertex(0, 0, 0);
        tmp_b_mesh.facets.set_vertex(0, 1, 0);
        tmp_b_mesh.facets.set_vertex(0, 2, 0);
    } else {
        tmp_b_mesh.vertices.clear();
        tmp_b_mesh.vertices.create_vertices((int) b_edges.size() * 2);
        int cnt = 0;
        for (auto &e:b_edges) {
            for (int j = 0; j < 2; j++) {
                GEO::vec3 &p = tmp_b_mesh.vertices.point(cnt++);
                p[0] = mesh.tet_vertices[e[j]].pos[0];
                p[1] = mesh.tet_vertices[e[j]].pos[1];
                p[2] = mesh.tet_vertices[e[j]].pos[2];
            }
        }
        tmp_b_mesh.facets.clear();
        tmp_b_mesh.facets.create_triangles((int) b_edges.size());
        for (int i = 0; i < b_edges.size(); i++) {
            tmp_b_mesh.facets.set_vertex(i, 0, i * 2);
            tmp_b_mesh.facets.set_vertex(i, 1, i * 2);
            tmp_b_mesh.facets.set_vertex(i, 2, i * 2 + 1);
        }
    }
    mesh_reorder(tmp_b_mesh, GEO::MESH_ORDER_MORTON);
    tmp_b_tree = std::make_shared<GEO::MeshFacetsAABBWithEps>(tmp_b_mesh);
}

void floatTetWild::AABBWrapper::init_tmp_b_mesh_and_tree(const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces,
                              const std::vector<std::array<int, 2>>& b_edges){
    if (b_edges.empty()) {
        tmp_b_mesh.vertices.clear();
        tmp_b_mesh.vertices.create_vertices(1);
        tmp_b_mesh.vertices.point(0) = GEO::vec3(0, 0, 0);
        tmp_b_mesh.facets.clear();
        tmp_b_mesh.facets.create_triangles(1);
        tmp_b_mesh.facets.set_vertex(0, 0, 0);
        tmp_b_mesh.facets.set_vertex(0, 1, 0);
        tmp_b_mesh.facets.set_vertex(0, 2, 0);
    } else {
        tmp_b_mesh.vertices.clear();
        tmp_b_mesh.vertices.create_vertices((int) b_edges.size() * 2);
        int cnt = 0;
        for (auto &e:b_edges) {
            for (int j = 0; j < 2; j++) {
                GEO::vec3 &p = tmp_b_mesh.vertices.point(cnt++);
                p[0] = input_vertices[e[j]][0];
                p[1] = input_vertices[e[j]][1];
                p[2] = input_vertices[e[j]][2];
            }
        }
        tmp_b_mesh.facets.clear();
        tmp_b_mesh.facets.create_triangles((int) b_edges.size());
        for (int i = 0; i < b_edges.size(); i++) {
            tmp_b_mesh.facets.set_vertex(i, 0, i * 2);
            tmp_b_mesh.facets.set_vertex(i, 1, i * 2);
            tmp_b_mesh.facets.set_vertex(i, 2, i * 2 + 1);
        }
    }
    mesh_reorder(tmp_b_mesh, GEO::MESH_ORDER_MORTON);
    tmp_b_tree = std::make_shared<GEO::MeshFacetsAABBWithEps>(tmp_b_mesh);
}