#include <floattetwild/LocalOperations.h>
#include <floattetwild/AABBWrapper.h>
#include <geogram/mesh/mesh_reorder.h>
#include <geogram/basic/geometry_nd.h>
#include <floattetwild/TriangleInsertion.h>



void floatTetWild::AABBWrapper::init_b_mesh_and_tree(const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces, Mesh& mesh) {
    b_mesh.clear(false, false);
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

    std::vector<std::pair<std::array<int, 2>, std::vector<int>>> _;
    std::vector<std::array<int, 2>> b_edges;
    std::vector<bool> _1;
    find_boundary_edges(input_vertices, input_faces,
                        std::vector<bool>(input_faces.size(), true), std::vector<bool>(input_faces.size(), true),
                        _, _1, b_edges);

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
    b_tree = std::make_shared<MeshFacetsAABBWithEps>(b_mesh);

    if(b_edges.empty())
        mesh.is_closed = true;

#ifdef NEW_ENVELOPE
    std::vector<Vector3> vs;
    std::vector<Vector3i> fs;
    if (b_edges.empty()) {
        vs.push_back(Vector3(0, 0, 0));
        fs.push_back(Vector3i(0, 0, 0));
    } else {
        vs.resize(b_edges.size() * 2);
        fs.resize(b_edges.size());
        for (int i = 0; i < b_edges.size(); i++) {
            vs[i * 2] = input_vertices[b_edges[i][0]];
            vs[i * 2 + 1] = input_vertices[b_edges[i][1]];
            fs[i] = Vector3i(i * 2, i * 2 + 1, i * 2 + 1);
        }
    }
//    b_tree_exact = std::make_shared<fastEnvelope::FastEnvelope>(vs, fs, eps);
    b_tree_exact.init(vs, fs, mesh.params.eps);
#endif
}

void floatTetWild::AABBWrapper::init_tmp_b_mesh_and_tree(const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces,
                              const std::vector<std::array<int, 2>>& b_edges1,
                              const Mesh& mesh, const std::vector<std::array<int, 2>>& b_edges2) {
    if (b_edges1.empty() && b_edges2.empty()) {
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
        tmp_b_mesh.vertices.create_vertices((int) (b_edges1.size() + b_edges2.size()) * 2);
        int cnt = 0;
        for (auto &e:b_edges1) {
            for (int j = 0; j < 2; j++) {
                GEO::vec3 &p = tmp_b_mesh.vertices.point(cnt++);
                p[0] = input_vertices[e[j]][0];
                p[1] = input_vertices[e[j]][1];
                p[2] = input_vertices[e[j]][2];
            }
        }
        for (auto &e:b_edges2) {
            for (int j = 0; j < 2; j++) {
                GEO::vec3 &p = tmp_b_mesh.vertices.point(cnt++);
                p[0] = mesh.tet_vertices[e[j]].pos[0];
                p[1] = mesh.tet_vertices[e[j]].pos[1];
                p[2] = mesh.tet_vertices[e[j]].pos[2];
            }
        }

        tmp_b_mesh.facets.clear();
        tmp_b_mesh.facets.create_triangles((int) b_edges1.size() + b_edges2.size());
        for (int i = 0; i < b_edges1.size(); i++) {
            tmp_b_mesh.facets.set_vertex(i, 0, i * 2);
            tmp_b_mesh.facets.set_vertex(i, 1, i * 2);
            tmp_b_mesh.facets.set_vertex(i, 2, i * 2 + 1);
        }
        for (int i = b_edges1.size(); i < b_edges1.size() + b_edges2.size(); i++) {
            tmp_b_mesh.facets.set_vertex(i, 0, i * 2);
            tmp_b_mesh.facets.set_vertex(i, 1, i * 2);
            tmp_b_mesh.facets.set_vertex(i, 2, i * 2 + 1);
        }
    }
    mesh_reorder(tmp_b_mesh, GEO::MESH_ORDER_MORTON);
    tmp_b_tree = std::make_shared<MeshFacetsAABBWithEps>(tmp_b_mesh);

#ifdef NEW_ENVELOPE
    std::vector<Vector3> vs;
    std::vector<Vector3i> fs;
    if (b_edges1.empty() && b_edges2.empty()) {
        vs.push_back(Vector3(0, 0, 0));
        fs.push_back(Vector3i(0, 0, 0));
    } else {
        vs.resize((b_edges1.size() + b_edges2.size()) * 2);
        fs.resize(b_edges1.size() + b_edges2.size());
        for (int i = 0; i < b_edges1.size(); i++) {
            vs[i * 2] = input_vertices[b_edges1[i][0]];
            vs[i * 2 + 1] = input_vertices[b_edges1[i][1]];
            fs[i] = Vector3i(i * 2, i * 2 + 1, i * 2 + 1);
        }
        for (int i = b_edges1.size(); i < b_edges1.size() + b_edges2.size(); i++) {
            vs[i * 2] = mesh.tet_vertices[b_edges2[i - b_edges1.size()][0]].pos;
            vs[i * 2 + 1] = mesh.tet_vertices[b_edges2[i - b_edges1.size()][1]].pos;
            fs[i] = Vector3i(i * 2, i * 2 + 1, i * 2 + 1);
        }
    }
//    tmp_b_tree_exact = std::make_shared<fastEnvelope::FastEnvelope>(vs, fs, mesh.params.eps_input);
    tmp_b_tree_exact.init(vs, fs, mesh.params.eps);
#endif
}
