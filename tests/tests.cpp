////////////////////////////////////////////////////////////////////////////////
#include <floattetwild/FloatTetwild.h>
#include <catch2/catch.hpp>

#include <floattetwild/AABBWrapper.h>
#include <floattetwild/FloatTetDelaunay.h>
#include <floattetwild/LocalOperations.h>
#include <floattetwild/MeshImprovement.h>
#include <floattetwild/Simplification.h>
#include <floattetwild/Statistics.h>
#include <floattetwild/TriangleInsertion.h>
#include <floattetwild/CSGTreeParser.hpp>
#include <floattetwild/Mesh.hpp>
#include <floattetwild/MeshIO.hpp>
#include <floattetwild/Predicates.hpp>

#include <Eigen/Dense>
#include <floattetwild/Logger.hpp>

#include <igl/Timer.h>

#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/common.h>
#include <geogram/basic/geometry.h>
#include <geogram/basic/logger.h>
#include <geogram/basic/numeric.h>

#include <geogram/mesh/mesh.h>

#include <bitset>

#include <future>
#include <iostream>
////////////////////////////////////////////////////////////////////////////////

// using namespace floattetwild;

int run_ftetwild(const GEO::Mesh& mesh, int& numVertices, int& numTets)
{
    GEO::Mesh local_mesh;
    local_mesh.copy(mesh);

    // reset random number generator
    //srand(1);

    numTets     = 0;
    numVertices = 0;

    // Set arguments
    floatTetWild::Parameters params;
    floatTetWild::Mesh       tmp_mesh;
    params = tmp_mesh.params;

    params.is_quiet    = true;
    params.log_level   = 3;
    params.stop_energy = 10.0;

    //    params.ideal_edge_length_rel = targetEdgeLength;
    // params.eps_rel               = maxTolerance;

    // Initialize TetWild logger if it doesn't exist yet
    if (!floatTetWild::Logger::logger_) {
        floatTetWild::Logger::init(!params.is_quiet, params.log_path);
        params.log_level = std::max(0, std::min(6, params.log_level));
        spdlog::set_level(static_cast<spdlog::level::level_enum>(params.log_level));
        spdlog::flush_every(std::chrono::seconds(3));
    }

    // Call FloatTetWild
    try {
        Eigen::MatrixXd VO;
        Eigen::MatrixXi TO;
        floatTetWild::tetrahedralization(local_mesh, params, VO, TO);
        numTets     = TO.rows();
        numVertices = VO.rows();
    }
    catch (const std::exception& e) {
        // General panic
        return -1;
    }

    return 0;
}

TEST_CASE("test-case", "[suite-name]")
{
    REQUIRE(1 == Approx(1).margin(1e-10));
    REQUIRE(1 == Approx(1).margin(1e-10));
}

TEST_CASE("Tetrahedralize mesh with ftetwild", "[ftetwild]")
{
    GEO::initialize();

    SECTION("Reproducibility 1")
    {
        std::vector<floatTetWild::Vector3> input_vertices;
        std::vector<Eigen::Vector3i>       input_faces;
        std::vector<int>                   input_tags;
        GEO::Mesh                          mesh;

        floatTetWild::MeshIO::load_mesh(
            "..\\..\\..\\tests\\obj\\bunny.obj", input_vertices, input_faces, mesh, input_tags);

        int  numV = -1, numT = -1;
        auto ret = run_ftetwild(mesh, numV, numT);
        REQUIRE(ret == 0);

        int  numV1, numT1;
        auto ret1 = run_ftetwild(mesh, numV1, numT1);
        REQUIRE(ret1 == 0);

        int  numV2, numT2;
        auto ret2 = run_ftetwild(mesh, numV2, numT2);
        REQUIRE(ret2 == 0);

        REQUIRE(numV1 == numV);
        REQUIRE(numT1 == numT);
        REQUIRE(numV2 == numV);
        REQUIRE(numT2 == numT);

    }

    SECTION("Reproducibility 2")
    {
        std::vector<floatTetWild::Vector3> input_vertices;
        std::vector<Eigen::Vector3i>       input_faces;
        std::vector<int>                   input_tags;
        GEO::Mesh                          mesh;

        floatTetWild::MeshIO::load_mesh(
            "..\\..\\..\\tests\\obj\\cylinder.obj", input_vertices, input_faces, mesh, input_tags);

        int  numV = -1, numT = -1;
        auto ret = run_ftetwild(mesh, numV, numT);
        REQUIRE(ret == 0);

        int  numV1, numT1;
        auto ret1 = run_ftetwild(mesh, numV1, numT1);
        REQUIRE(ret1 == 0);

        int  numV2, numT2;
        auto ret2 = run_ftetwild(mesh, numV2, numT2);
        REQUIRE(ret2 == 0);

        REQUIRE(numV1 == numV);
        REQUIRE(numT1 == numT);
        REQUIRE(numV2 == numV);
        REQUIRE(numT2 == numT);

    }

    SECTION("Reproducibility 3")
    {
        std::vector<floatTetWild::Vector3> input_vertices;
        std::vector<Eigen::Vector3i>       input_faces;
        std::vector<int>                   input_tags;
        GEO::Mesh                          mesh;

        floatTetWild::MeshIO::load_mesh(
            "..\\..\\..\\tests\\obj\\box.stl", input_vertices, input_faces, mesh, input_tags);

        int  numV = -1, numT = -1;
        auto ret = run_ftetwild(mesh, numV, numT);
        REQUIRE(ret == 0);

        int  numV1, numT1;
        auto ret1 = run_ftetwild(mesh, numV1, numT1);
        REQUIRE(ret1 == 0);

        int  numV2, numT2;
        auto ret2 = run_ftetwild(mesh, numV2, numT2);
        REQUIRE(ret2 == 0);

        REQUIRE(numV1 == numV);
        REQUIRE(numT1 == numT);
        REQUIRE(numV2 == numV);
        REQUIRE(numT2 == numT);

    }

    SECTION("Reproducibility 4")
    {
        std::vector<floatTetWild::Vector3> input_vertices;
        std::vector<Eigen::Vector3i>       input_faces;
        std::vector<int>                   input_tags;
        GEO::Mesh                          mesh;

        floatTetWild::MeshIO::load_mesh(
            "..\\..\\..\\tests\\obj\\cylinder.stl", input_vertices, input_faces, mesh, input_tags);

        int  numV = -1, numT = -1;
        auto ret = run_ftetwild(mesh, numV, numT);
        REQUIRE(ret == 0);

        int  numV1, numT1;
        auto ret1 = run_ftetwild(mesh, numV1, numT1);
        REQUIRE(ret1 == 0);

        int  numV2, numT2;
        auto ret2 = run_ftetwild(mesh, numV2, numT2);
        REQUIRE(ret2 == 0);

        REQUIRE(numV1 == numV);
        REQUIRE(numT1 == numT);
        REQUIRE(numV2 == numV);
        REQUIRE(numT2 == numT);

    }

    SECTION("Thread safety")
    {
        // TODO: FIX ftetwild not thread safe
        // auto future1 = std::async(std::launch::async, run_ftetwild);
        // auto future2 = std::async(std::launch::async, run_ftetwild);

        // auto ret1 = future1.get();
        // auto ret2 = future2.get();
        // REQUIRE(ret1 == 0);
        // REQUIRE(ret2 == 0);
    }
}
