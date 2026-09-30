// Standalone reproducer for a nondeterministic bus error / segfault in
// tetrahedralization(). Root cause: collapsing()'s remove_an_edge() lambda
// in src/Simplification.cpp only updates conn_fs bookkeeping inside
// `#ifndef FLOAT_TETWILD_USE_TBB` guards, so under TBB (the default build)
// conn_fs can drift out of sync with input_faces. simplify()'s post-collapse
// vertex cleanup then trusted conn_fs[i].empty() to decide whether vertex i
// survives; when that's stale, a vertex still referenced by a live face gets
// mapped to -1, corrupting input_faces -- which crashes later (observed in
// find_boundary_edges()'s conn_tris[input_faces[i][j]].push_back()).
//
// Usage: repro_dragon <mesh.ply> [iters]
//
// Needs a real-world mesh with enough vertices/triangles for TBB to have
// meaningful parallelism during collapsing() -- we could not reproduce this
// on small/decimated meshes (tried both uniformly decimated and
// geometry-preserving cropped regions, up to ~20K faces, hundreds of
// iterations, zero repros). It reproduces reliably (order of 1-in-5 to
// 1-in-10 runs) on the public "xyzrgb_dragon" model
// (https://graphics.stanford.edu/data/3Dscanrep/, ~360K vertices / 720K
// faces); building with -DFLOAT_TETWILD_ENABLE_TBB=OFF avoids the crash
// entirely (7/7 clean runs in our testing, vs. crashes with TBB on).
#include <floattetwild/FloatTetwild.h>
#include <floattetwild/Logger.hpp>
#include <floattetwild/MeshIO.hpp>
#include <floattetwild/Parameters.h>
#include <floattetwild/Types.hpp>

#include <igl/readPLY.h>

#include <geogram/basic/common.h>
#include <geogram/mesh/mesh.h>

#include <Eigen/Core>
#include <cstdlib>
#include <iostream>

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cerr << "usage: repro_dragon <mesh.ply> [iters]" << std::endl;
        return 1;
    }
    int iters = argc >= 3 ? std::atoi(argv[2]) : 1;

    ::setenv("GEO_NO_SIGNAL_HANDLER", "1", 1);
    GEO::initialize();

    Eigen::MatrixXd VI;
    Eigen::MatrixXi FI;
    Eigen::MatrixXd N, UV;
    if (!igl::readPLY(argv[1], VI, FI, N, UV)) {
        std::cerr << "failed to read " << argv[1] << std::endl;
        return 1;
    }
    std::cout << "VI: " << VI.rows() << " x " << VI.cols()
              << ", FI: " << FI.rows() << " x " << FI.cols() << std::endl;

    for (int it = 0; it < iters; it++) {
        std::cout << "=== iteration " << it << " ===" << std::endl;

        std::vector<floatTetWild::Vector3> input_vertices(VI.rows());
        for (int i = 0; i < VI.rows(); i++) {
            input_vertices[i] = VI.row(i).transpose();
        }
        std::vector<floatTetWild::Vector3i> input_faces(FI.rows());
        for (int i = 0; i < FI.rows(); i++) {
            input_faces[i] = FI.row(i).transpose();
        }
        std::vector<int> input_flags;

        GEO::Mesh sf_mesh;
        floatTetWild::MeshIO::load_mesh(input_vertices, input_faces, sf_mesh, input_flags);

        floatTetWild::Logger::init(true, "");
        spdlog::set_level(spdlog::level::info);

        floatTetWild::Parameters params;
        Eigen::MatrixXd          VO;
        Eigen::MatrixXi          TO;
        int ret = floatTetWild::tetrahedralization(sf_mesh, params, VO, TO);
        std::cout << "ret=" << ret << " nV=" << VO.rows() << " nT=" << TO.rows() << std::endl;
    }

    return 0;
}
