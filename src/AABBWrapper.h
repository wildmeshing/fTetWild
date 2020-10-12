#pragma once

#include <floattetwild/mesh_AABB.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_geometry.h>
#include <floattetwild/Mesh.hpp>

#include <memory>

//#define NEW_ENVELOPE //fortest

#ifdef NEW_ENVELOPE
#include <fastenvelope/FastEnvelope.h>
#endif

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

#ifdef NEW_ENVELOPE
        fastEnvelope::FastEnvelope b_tree_exact;
        fastEnvelope::FastEnvelope tmp_b_tree_exact;
        fastEnvelope::FastEnvelope sf_tree_exact;
        fastEnvelope::FastEnvelope sf_tree_exact_simplify;
//        std::shared_ptr<fastEnvelope::FastEnvelope> b_tree_exact;
//        std::shared_ptr<fastEnvelope::FastEnvelope> tmp_b_tree_exact;
//        std::shared_ptr<fastEnvelope::FastEnvelope> sf_tree_exact;

        inline void init_sf_tree(const std::vector<Vector3> &vs, const std::vector<Vector3i> &fs, double eps) {
//            sf_tree_exact = std::make_shared<fastEnvelope::FastEnvelope>(vs, fs, eps);
            sf_tree_exact.init(vs, fs, eps);
            sf_tree_exact_simplify.init(vs, fs, 0.8*eps);
        }
        inline void init_sf_tree(const std::vector<Vector3> &vs, const std::vector<Vector3i> &fs,
                std::vector<double>& eps, double bbox_diag_length) {
            for (auto& e: eps)
                e *= bbox_diag_length;
            sf_tree_exact.init(vs, fs, eps);
            std::vector<double> eps_simplify = eps;
            for (auto& e: eps_simplify)
                e *= 0.8;
            sf_tree_exact_simplify.init(vs, fs, eps_simplify);
        }
#endif

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

#ifdef NEW_ENVELOPE
        inline bool is_out_sf_envelope_exact(const std::array<Vector3, 3> &triangle) const {
            return sf_tree_exact.is_outside(triangle);
        }

        inline bool is_out_sf_envelope_exact_simplify(const std::array<Vector3, 3> &triangle) const {
            return sf_tree_exact_simplify.is_outside(triangle);
        }

        inline bool is_out_b_envelope_exact(const std::array<Vector3, 3> &triangle) const {
            return b_tree_exact.is_outside(triangle);
        }

        inline bool is_out_tmp_b_envelope_exact(const std::array<Vector3, 3> &triangle) const {
            return tmp_b_tree_exact.is_outside(triangle);
        }
#endif

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

        inline bool is_out_sf_envelope(const Vector3& p, const Scalar eps_2,
                                       GEO::index_t& prev_facet, double& sq_dist, GEO::vec3& nearest_p) const {
            GEO::vec3 geo_p(p[0], p[1], p[2]);
            return is_out_sf_envelope(geo_p, eps_2, prev_facet, sq_dist, nearest_p);
        }
        inline bool is_out_sf_envelope(const GEO::vec3& geo_p, const Scalar eps_2,
                                       GEO::index_t& prev_facet, double& sq_dist, GEO::vec3& nearest_p) const {
            if (prev_facet != GEO::NO_FACET) {
                get_point_facet_nearest_point(sf_mesh, geo_p, prev_facet, nearest_p, sq_dist);
            }
            if (Scalar(sq_dist) > eps_2) {
                sf_tree.facet_in_envelope_with_hint(geo_p, eps_2, prev_facet, nearest_p, sq_dist);
            }

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

#ifdef NEW_ENVELOPE
        inline bool is_out_sf_envelope_exact(const Vector3& p) const {
            return sf_tree_exact.is_outside(p);
        }

        inline bool is_out_b_envelope_exact(const Vector3& p) const {
            return b_tree_exact.is_outside(p);
        }

        inline bool is_out_tmp_b_envelope_exact(const Vector3& p) const {
            return tmp_b_tree_exact.is_outside(p);
        }
#endif


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
