#pragma once

#include <floattetwild/mesh_AABB.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_geometry.h>
#include <floattetwild/Mesh.hpp>

#include <memory>

namespace floatTetWild {
    class AABBWrapper {
    public:
        GEO::Mesh b_mesh;
        GEO::Mesh tmp_b_mesh;
        const GEO::Mesh &sf_mesh;

        std::shared_ptr<MeshFacetsAABBWithEps> b_tree;
        std::shared_ptr<MeshFacetsAABBWithEps> tmp_b_tree;
        MeshFacetsAABBWithEps sf_tree;

        //// initialization
        inline Scalar get_sf_diag() const { return GEO::bbox_diagonal(sf_mesh); }

        AABBWrapper(const GEO::Mesh &sf_mesh) : sf_mesh(sf_mesh), sf_tree(sf_mesh) {}

        void init_b_mesh_and_tree(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces, Mesh &mesh);

        void init_tmp_b_mesh_and_tree(const std::vector<Vector3> &input_vertices,
                                      const std::vector<Vector3i> &input_faces,
                                      const std::vector<std::array<int, 2>> &b_edges1,
                                      const Mesh &mesh, const std::vector<std::array<int, 2>> &b_edges2);

        //// projection
        inline Scalar project_to_sf(Vector3 &p) const {
            GEO::vec3 geo_p(p[0], p[1], p[2]);
            GEO::vec3 nearest_p;
            double sq_dist = std::numeric_limits<double>::max(); //??
            sf_tree.nearest_facet(geo_p, nearest_p, sq_dist);
            p[0] = nearest_p[0];
            p[1] = nearest_p[1];
            p[2] = nearest_p[2];

            return sq_dist;
        }

        inline Scalar project_to_b(Vector3 &p) const {
            GEO::vec3 geo_p(p[0], p[1], p[2]);
            GEO::vec3 nearest_p;
            double sq_dist = std::numeric_limits<double>::max(); //?
            b_tree->nearest_facet(geo_p, nearest_p, sq_dist);
            p[0] = nearest_p[0];
            p[1] = nearest_p[1];
            p[2] = nearest_p[2];

            return sq_dist;
        }

        inline Scalar project_to_tmp_b(Vector3 &p) const {
            GEO::vec3 geo_p(p[0], p[1], p[2]);
            GEO::vec3 nearest_p;
            double sq_dist = std::numeric_limits<double>::max(); //?
            tmp_b_tree->nearest_facet(geo_p, nearest_p, sq_dist);
            p[0] = nearest_p[0];
            p[1] = nearest_p[1];
            p[2] = nearest_p[2];

            return sq_dist;
        }

        inline int get_nearest_face_sf(const Vector3 &p) const {
            GEO::vec3 geo_p(p[0], p[1], p[2]);
            GEO::vec3 nearest_p;
            double sq_dist = std::numeric_limits<double>::max(); //??
            return sf_tree.nearest_facet(geo_p, nearest_p, sq_dist);
        }

        inline Scalar get_sq_dist_to_sf(const Vector3 &p) const {
            GEO::vec3 geo_p(p[0], p[1], p[2]);
            GEO::vec3 nearest_p;
            double sq_dist = std::numeric_limits<double>::max(); //??
            sf_tree.nearest_facet(geo_p, nearest_p, sq_dist);
            return sq_dist;
        }

        //// envelope check - triangle
        inline bool is_out_sf_envelope(const std::vector<GEO::vec3> &ps, const Scalar eps_2,
                                       GEO::index_t prev_facet = GEO::NO_FACET) const {
            GEO::vec3 nearest_point;
            double sq_dist = std::numeric_limits<double>::max();

            for (const GEO::vec3 &current_point : ps) {
                if (prev_facet != GEO::NO_FACET) {
                    get_point_facet_nearest_point(sf_mesh, current_point, prev_facet, nearest_point, sq_dist);
                }
                if (Scalar(sq_dist) > eps_2) {
                    sf_tree.facet_in_envelope_with_hint(current_point, eps_2, prev_facet, nearest_point, sq_dist);
                }
                if (Scalar(sq_dist) > eps_2) {
                    return true;
                }
            }

            return false;
        }

        inline bool is_out_b_envelope(const std::vector<GEO::vec3> &ps, const Scalar eps_2,
                                      GEO::index_t prev_facet = GEO::NO_FACET) const {
            GEO::vec3 nearest_point;
            double sq_dist = std::numeric_limits<double>::max();

            for (const GEO::vec3 &current_point : ps) {
                if (prev_facet != GEO::NO_FACET) {
                    get_point_facet_nearest_point(b_mesh, current_point, prev_facet, nearest_point, sq_dist);
                }
                if (Scalar(sq_dist) > eps_2) {
                    b_tree->facet_in_envelope_with_hint(current_point, eps_2, prev_facet, nearest_point, sq_dist);
                }
                if (Scalar(sq_dist) > eps_2) {
                    return true;
                }
            }

            return false;
        }

        inline bool is_out_tmp_b_envelope(const std::vector<GEO::vec3> &ps, const Scalar eps_2,
                                          GEO::index_t prev_facet = GEO::NO_FACET) const {
            GEO::vec3 nearest_point;
            double sq_dist = std::numeric_limits<double>::max();

            for (const GEO::vec3 &current_point : ps) {
                if (prev_facet != GEO::NO_FACET) {
                    get_point_facet_nearest_point(tmp_b_mesh, current_point, prev_facet, nearest_point, sq_dist);
                }
                if (Scalar(sq_dist) > eps_2) {
                    tmp_b_tree->facet_in_envelope_with_hint(current_point, eps_2, prev_facet, nearest_point, sq_dist);
                }
                if (Scalar(sq_dist) > eps_2) {
                    return true;
                }
            }

            return false;
        }

        //// envelope check - point
        inline bool is_out_sf_envelope(const Vector3 &p, const Scalar eps_2, GEO::index_t &prev_facet) const {
            GEO::vec3 nearest_p;
            double sq_dist;
            GEO::vec3 geo_p(p[0], p[1], p[2]);
            prev_facet = sf_tree.facet_in_envelope(geo_p, eps_2, nearest_p, sq_dist);

            if (Scalar(sq_dist) > eps_2)
                return true;
            return false;
        }

        inline bool is_out_b_envelope(const Vector3 &p, const Scalar eps_2, GEO::index_t &prev_facet) const {
            GEO::vec3 nearest_p;
            double sq_dist;
            GEO::vec3 geo_p(p[0], p[1], p[2]);
            prev_facet = b_tree->facet_in_envelope(geo_p, eps_2, nearest_p, sq_dist);

            if (Scalar(sq_dist) > eps_2)
                return true;
            return false;
        }

        inline bool is_out_tmp_b_envelope(const Vector3 &p, const Scalar eps_2, GEO::index_t &prev_facet) const {
            GEO::vec3 nearest_p;
            double sq_dist;
            GEO::vec3 geo_p(p[0], p[1], p[2]);
            prev_facet = tmp_b_tree->facet_in_envelope(geo_p, eps_2, nearest_p, sq_dist);

            if (Scalar(sq_dist) > eps_2)
                return true;
            return false;
        }


        //fortest
        inline Scalar dist_sf_envelope(const std::vector<GEO::vec3> &ps, const Scalar eps_2,
                                       GEO::index_t prev_facet = GEO::NO_FACET) const {///only used for checking correctness
            GEO::vec3 nearest_point;
            double sq_dist = std::numeric_limits<double>::max();

            for (const GEO::vec3 &current_point : ps) {
                if (prev_facet != GEO::NO_FACET) {
                    get_point_facet_nearest_point(sf_mesh, current_point, prev_facet, nearest_point, sq_dist);
                }
                if (Scalar(sq_dist) > eps_2) {
                    sf_tree.facet_in_envelope_with_hint(current_point, eps_2, prev_facet, nearest_point, sq_dist);
                }
                if (Scalar(sq_dist) > eps_2) {
                    return sq_dist;
                }
            }

            return 0;
        }
        //fortest

    };

}
