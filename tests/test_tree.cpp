////////////////////////////////////////////////////////////////////////////////
#include <floattetwild/MeshIO.hpp>
#include <floattetwild/AABB.hpp>


#include <igl/Timer.h>

#include <geogram/basic/logger.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>

#include <geogram/mesh/mesh_AABB.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_geometry.h>

#include <catch2/catch_all.hpp>
#include <iostream>
////////////////////////////////////////////////////////////////////////////////

using namespace floatTetWild;



TEST_CASE("point_triangle_squared_distance", "[tree]") {
	double geo_dist;
	double dist;
	GEO::vec3 geo_nearest_p;
	Vector3 nearest_p;


	{
		const GEO::vec3 p1 = GEO::vec3(-46.4756, 3.75375, 9.27062e-15);
		const GEO::vec3 p2 = GEO::vec3(-46.3798, 2.03115, 1.3222);
		const GEO::vec3 p3 = GEO::vec3(-46.4294, 3.84615, 1.3695);
		const GEO::vec3 p = GEO::vec3(-130.681, -13.2918, 12.0658);
		double lambda1, lambda2, lambda3;  // barycentric coords, not used.
		geo_dist = GEO::Geom::point_triangle_squared_distance(p, p1, p2, p3, geo_nearest_p, lambda1, lambda2, lambda3);
	}


	{
		const Vector3 p1(-46.4756, 3.75375, 9.27062e-15);
		const Vector3 p2(-46.3798, 2.03115, 1.3222);
		const Vector3 p3(-46.4294, 3.84615, 1.3695);
		const Vector3 p(-130.681, -13.2918, 12.0658);
		dist = AABB::point_triangle_squared_distance(p, p1, p2, p3, nearest_p);
	}


	REQUIRE(dist == Approx(geo_dist).margin(SCALAR_ZERO*1e-2));
	REQUIRE(nearest_p[0] == Approx(geo_nearest_p[0]).margin(SCALAR_ZERO*1e-2));
	REQUIRE(nearest_p[1] == Approx(geo_nearest_p[1]).margin(SCALAR_ZERO*1e-2));
	REQUIRE(nearest_p[2] == Approx(geo_nearest_p[2]).margin(SCALAR_ZERO*1e-2));
}

TEST_CASE("aabb", "[tree]") {
#ifndef WIN32
	setenv("GEO_NO_SIGNAL_HANDLER", "1", 1);
#endif

	GEO::initialize();


	// Import standard command line arguments, and custom ones
	GEO::CmdLine::import_arg_group("standard");
	GEO::CmdLine::import_arg_group("pre");
	GEO::CmdLine::import_arg_group("algo");

	Parameters params;

	const std::string input_surface_path = "C:\\Users\\3Dscanner\\Downloads\\Octocat-v2.stl";
	std::vector<Vector3> vertices;
	std::vector<Vector3i> faces;

	GEO::Mesh mesh;
	bool ok = MeshIO::load_mesh(input_surface_path, vertices, faces, mesh);

	if (!ok) {
		INFO("Unable to load mesh")
			return;
	}

	//GEO::MeshFacetsAABB geo_tree(mesh, false);
	AABB tree(vertices, faces);

	Box bbox;
	for (const auto &v : vertices)
		bbox.extend(v);

	igl::Timer timer;

	const int n_pts = 1000000;
	std::vector<Vector3> pts(n_pts);
	//std::vector<GEO::vec3> geo_pts(n_pts);

	//std::vector<Vector3> expected(n_pts);
	//std::vector<int> expected_indices(n_pts);
	//std::vector<Scalar> expected_dst(n_pts);

	pts[0] << -130.68095686117647, -13.29181261183755, 12.065754038319332;
	pts[1] << 15.051636485920454334, -133.67956133032245702, -30.049805265677228761;
	//geo_pts[0] = GEO::vec3(pts[0][0], pts[0][1], pts[0][2]);
	//geo_pts[1] = GEO::vec3(pts[1][0], pts[1][1], pts[1][2]);

	for (int i = 2; i < n_pts; ++i)
	{
		pts[i] = Vector3::Random(3, 1).array()*(bbox.max - bbox.min).array() + bbox.min.array();
	//	geo_pts[i] = GEO::vec3(pts[i][0], pts[i][1], pts[i][2]);
	}

	//// std::cout.precision(20);
	//// std::cout<<pts.back()<<std::endl;
	//// std::cout<<geo_pts.back()<<std::endl;

	//GEO::vec3 geo_nearest_point;
	//double sq_dist_geo;
	Scalar sq_dist;
	//GEO::index_t geo_index;

	//INFO("running geogram aabb")
	//	timer.start();
	//for (int i = 0; i < n_pts; ++i) {
	//	geo_index = geo_tree.nearest_facet(geo_pts[i], geo_nearest_point, sq_dist_geo);
	//}
	//timer.stop();
	//INFO("Done!")
	//	std::cout << "geogram " << timer.getElapsedTimeInSec() << std::endl;


	//for (int i = 0; i < n_pts; ++i) {
	//	geo_index = geo_tree.nearest_facet(geo_pts[i], geo_nearest_point, sq_dist_geo);

	//	expected[i] << geo_nearest_point[0], geo_nearest_point[1], geo_nearest_point[2];
	//	expected_dst[i] = sq_dist_geo;
	//	expected_indices[i] = geo_index;
	//}

	Vector3 nearest_point;
	int index;

	INFO("running new aabb")
		timer.start();
	for (int i = 0; i < n_pts; ++i) {
		index = tree.nearest_facet(pts[i], nearest_point, sq_dist);
	}
	timer.stop();
	INFO("Done!")
	//	std::cout << "us " << timer.getElapsedTimeInSec() << std::endl;



	//for (int i = 0; i < n_pts; ++i) {
	//	index = tree.nearest_facet(pts[i], nearest_point, sq_dist);

	//	// REQUIRE(index == expected_indices[i]);
	//	REQUIRE(sq_dist == Approx(expected_dst[i]).margin(SCALAR_ZERO*1e-2));
	//	REQUIRE(nearest_point[0] == Approx(expected[i][0]).margin(SCALAR_ZERO*1e-2));
	//	REQUIRE(nearest_point[1] == Approx(expected[i][1]).margin(SCALAR_ZERO*1e-2));
	//	REQUIRE(nearest_point[2] == Approx(expected[i][2]).margin(SCALAR_ZERO*1e-2));
	//}

}
