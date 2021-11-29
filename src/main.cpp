// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include <CLI/CLI.hpp>

#ifdef FLOAT_TETWILD_USE_TBB
#include <tbb/task_scheduler_init.h>
#include <thread>
#endif

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

#include <Eigen/Dense>
#include <floattetwild/Logger.hpp>

#include <igl/Timer.h>
#include <igl/write_triangle_mesh.h>

#ifdef LIBIGL_WITH_TETGEN
#include <igl/copyleft/tetgen/tetrahedralize.h>
#endif

#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/logger.h>
#include <geogram/mesh/mesh.h>
#include <bitset>

using namespace floatTetWild;
using namespace Eigen;

class GeoLoggerForward : public GEO::LoggerClient
{
    std::shared_ptr<spdlog::logger> logger_;

  public:
    template<typename T>
    GeoLoggerForward(T logger)
        : logger_(logger)
    {}

  private:
    std::string truncate(const std::string& msg)
    {
        static size_t prefix_len = GEO::CmdLine::ui_feature(" ", false).size();
        return msg.substr(prefix_len, msg.size() - 1 - prefix_len);
    }

  protected:
    void div(const std::string& title) override
    {
        logger_->trace(title.substr(0, title.size() - 1));
    }

    void out(const std::string& str) override { logger_->info(truncate(str)); }

    void warn(const std::string& str) override { logger_->warn(truncate(str)); }

    void err(const std::string& str) override { logger_->error(truncate(str)); }

    void status(const std::string& str) override
    {
        // Errors and warnings are also dispatched as status by geogram, but without
        // the "feature" header. We thus forward them as trace, to avoid duplicated
        // logger info...
        logger_->trace(str.substr(0, str.size() - 1));
    }
};

#include <geogram/basic/common.h>
#include <geogram/basic/geometry.h>
#include <geogram/basic/numeric.h>
#include <floattetwild/Predicates.hpp>

#include <floattetwild/MshLoader.h>
#include <geogram/mesh/mesh_AABB.h>

void connect_2_meshes(std::string m1, std::string m2, std::string m);

// extern "C" void exactinit();
int main(int argc, char** argv)
{
#ifdef STORE_SAMPLE_POINTS
    cout << "STORE_SAMPLE_POINTS defined" << endl;
#endif

    //    connect_2_meshes("/Users/yixinhu/Downloads/test_cutting/100729.stl",
    //                     "/Users/yixinhu/Downloads/test_cutting/37627.stl");

    //    connect_2_meshes("/Users/yixinhu/Downloads/stamp.stl",
    //                     "/Users/yixinhu/Downloads/cookies.stl", "stamp-cookies");

    //    connect_2_meshes("/Users/yixinhu/Downloads/cylinder.stl",
    //                     "/Users/yixinhu/Downloads/metatron.stl", "cylinder-metatron");

    //    return 0;
    //    igl::Timer test_timer;
    //    int n = 1e8;
    //    std::vector<std::array<int, 2>> edges(n);
    //    for(int i=0;i<n;i++){
    ////        edges[i] = {{rand() % 30, rand() % 30}};
    //        edges[i] = {{rand() % 7, rand() % 7}};
    //    }
    //    test_timer.start();
    //    int sum = 0;
    //    for(int i=0;i<n;i++){
    ////        edges[i] = {{edges[i][1], edges[i][0]}};
    //        sum += edges[i][0]%4;
    ////        for(int j=0;j<4;j++)
    ////            sum+=edges[i][mod2(j)];
    //    }
    //    cout<<sum<<endl;
    //    cout<<test_timer.getElapsedTime()<<endl;
    //
    //    for(int i=0;i<n;i++){
    ////        edges[i] = {{rand() % 30, rand() % 30}};
    //        edges[i] = {{rand() % 7, rand() % 7}};
    //    }
    //    test_timer.start();
    //    int sum1 = 0;
    //    for(int i=0;i<n;i++){
    ////        std::swap(edges[i][0], edges[i][1]);
    //        sum1 += mod4(edges[i][0]);
    ////        for(short int j=0;j<4;j++)
    ////            sum1+=edges[i][mod2(j)];
    //    }
    //    cout<<sum1<<endl;
    //    cout<<test_timer.getElapsedTime()<<endl;
    //
    //    return 0;

#ifndef WIN32
    setenv("GEO_NO_SIGNAL_HANDLER", "1", 1);
#endif

    GEO::initialize();
    //    exactinit();

    std::vector<int> indices(20);
    std::iota(std::begin(indices), std::end(indices), 0);
    floatTetWild::Random::shuffle(indices);
    for (int a : indices)
        std::cout << a << " ";
    std::cout << std::endl;

    // Import standard command line arguments, and custom ones
    GEO::CmdLine::import_arg_group("standard");
    GEO::CmdLine::import_arg_group("pre");
    GEO::CmdLine::import_arg_group("algo");

    bool run_tet_gen   = false;
    bool skip_simplify = false;
    bool nobinary      = false;
    bool nocolor       = false;
    bool export_raw    = false;

    Mesh        mesh;
    Parameters& params = mesh.params;

    CLI::App command_line {"float-tetwild"};
    command_line
      .add_option("-i,--input",
                  params.input_path,
                  "Input surface mesh INPUT in .off/.obj/.stl/.ply format. (string, required)")
      ->check(CLI::ExistingFile);
    command_line.add_option("-o,--output",
                            params.output_path,
                            "Output tetmesh OUTPUT in .msh format. (string, optional, default: "
                            "input_file+postfix+'.msh')");

    command_line.add_option("--tag", params.tag_path, "Tag input faces for Boolean operation.")
      ->check(CLI::ExistingFile);
    //    const int UNION = 0;
    //    const int INTERSECTION = 1;
    //    const int DIFFERENCE = 2;
    int         boolean_op = -1;
    std::string csg_file   = "";
    command_line.add_option(
      "--op", boolean_op, "Boolean operation: 0: union, 1: intersection, 2: difference.");

    command_line.add_option(
      "-l,--lr",
      params.ideal_edge_length_rel,
      "ideal_edge_length = diag_of_bbox * L. (double, optional, default: 0.05)");
    command_line.add_option("-e,--epsr",
                            params.eps_rel,
                            "epsilon = diag_of_bbox * EPS. (double, optional, default: 1e-3)");

    command_line.add_option("--max-its", params.max_its, "(for debugging usage only)");
    command_line.add_option(
      "--stop-energy", params.stop_energy, "Stop optimization when max energy is lower than this.");
    command_line.add_option("--stage", params.stage, "(for debugging usage only)");
    command_line.add_option("--stop-p", params.stop_p, "(for debugging usage only)");

    command_line.add_option("--postfix", params.postfix, "(for debugging usage only)");
    command_line.add_option("--log", params.log_path, "Log info to given file.");
    command_line.add_option("--level", params.log_level, "Log level (0 = most verbose, 6 = off).");

    command_line.add_flag("-q,--is-quiet", params.is_quiet, "Mute console output. (optional)");
    command_line.add_flag("--skip-simplify", skip_simplify, "skip preprocessing.");
    command_line.add_flag("--no-binary", nobinary, "export meshes as ascii");
    command_line.add_flag("--no-color", nocolor, "don't export color");
    command_line.add_flag("--not-sort-input", params.not_sort_input, "(for debugging usage only)");
    command_line.add_flag("--correct-surface-orientation",
                          params.correct_surface_orientation,
                          "(for debugging usage only)");

    command_line.add_option("--envelope-log", params.envelope_log, "(for debugging usage only)");
    command_line.add_flag(
      "--smooth-open-boundary", params.smooth_open_boundary, "Smooth the open boundary.");
    command_line.add_flag("--export-raw", export_raw, "Export raw output.");
    command_line.add_flag(
      "--manifold-surface", params.manifold_surface, "Force the output to be manifold.");
    command_line.add_flag("--coarsen", params.coarsen, "Coarsen the output as much as possible.");
    command_line.add_option("--csg", csg_file, "json file containg a csg tree")
      ->check(CLI::ExistingFile);

    command_line.add_flag(
      "--use-old-energy", floatTetWild::use_old_energy, "(for debugging usage only)");  // tmp

    command_line.add_flag(
      "--disable-filtering", params.disable_filtering, "Disable filtering out outside elements.");
    command_line.add_flag(
      "--use-floodfill", params.use_floodfill, "Use flood-fill to extract interior volume.");
    command_line.add_flag("--use-general-wn", params.use_general_wn, "Use general winding number.");
    command_line.add_flag(
      "--use-input-for-wn", params.use_input_for_wn, "Use input surface for winding number.");

    std::string background_mesh = "";
    command_line
      .add_option("--bg-mesh", background_mesh, "Background mesh for sizing field (.msh file).")
      ->check(CLI::ExistingFile);

#ifdef NEW_ENVELOPE
    std::string epsr_tags;
    command_line
      .add_option("--epsr-tags", epsr_tags, "List of envelope size for each input faces.")
      ->check(CLI::ExistingFile);
#endif

#ifdef LIBIGL_WITH_TETGEN
    command_line.add_flag("--tetgen", run_tet_gen, "run tetgen too. (optional)");
#endif
    unsigned int max_threads = std::numeric_limits<unsigned int>::max();
#ifdef FLOAT_TETWILD_USE_TBB
    command_line.add_option("--max-threads", max_threads, "Maximum number of threads used");
#endif

    try {
        command_line.parse(argc, argv);
    }
    catch (const CLI::ParseError& e) {
        return command_line.exit(e);
    }

#ifdef FLOAT_TETWILD_USE_TBB
    const size_t MB          = 1024 * 1024;
    const size_t stack_size  = 64 * MB;
    unsigned int num_threads = std::max(1u, std::thread::hardware_concurrency());
    num_threads              = std::min(max_threads, num_threads);
    params.num_threads       = num_threads;
    std::cout << "TBB threads " << num_threads << std::endl;
    tbb::task_scheduler_init scheduler(num_threads, stack_size);
#endif

    //    if(params.is_quiet){
    //        std::streambuf *orig_buf = cout.rdbuf();
    //        cout.rdbuf(NULL);
    //    }

    Logger::init(!params.is_quiet, params.log_path);
    params.log_level = std::max(0, std::min(6, params.log_level));
    spdlog::set_level(static_cast<spdlog::level::level_enum>(params.log_level));
    spdlog::flush_every(std::chrono::seconds(3));

    GEO::Logger* geo_logger = GEO::Logger::instance();
    geo_logger->unregister_all_clients();
    geo_logger->register_client(new GeoLoggerForward(logger().clone("geogram")));
    geo_logger->set_pretty(false);

    if (params.output_path.empty())
        params.output_path = params.input_path;
    if (params.log_path.empty())
        params.log_path = params.output_path;

    std::string output_mesh_name = params.output_path;
    if (params.output_path.size() > 3 &&
        params.output_path.substr(params.output_path.size() - 3, params.output_path.size()) ==
          "msh")
        output_mesh_name = params.output_path;
    else if (params.output_path.size() > 4 &&
             params.output_path.substr(params.output_path.size() - 4, params.output_path.size()) ==
               "mesh")
        output_mesh_name = params.output_path;
    else
        output_mesh_name = params.output_path + "_" + params.postfix + ".msh";

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// set sizing field
    Eigen::VectorXd V_in;
    Eigen::VectorXi T_in;
    Eigen::VectorXd values;
    if (!background_mesh.empty()) {
        PyMesh::MshLoader mshLoader(background_mesh);
        V_in   = mshLoader.get_nodes();
        T_in   = mshLoader.get_elements();
        values = mshLoader.get_node_field("values");
    }
    if (V_in.rows() != 0 && T_in.rows() != 0 && values.rows() != 0) {
        params.apply_sizing_field = true;

        params.V_sizing_field = V_in;
        params.T_sizing_field = T_in;
        params.values_sizing_field = values;
    }

    /// set input tage
    std::vector<Vector3>  input_vertices;
    std::vector<Vector3i> input_faces;
    std::vector<int>      input_tags;
    //    std::vector<double> input_epsr_tags;

    if (!params.tag_path.empty()) {
        input_tags.reserve(input_faces.size());
        std::string   line;
        std::ifstream fin(params.tag_path);
        if (fin.is_open()) {
            while (getline(fin, line)) {
                input_tags.push_back(std::stoi(line));
            }
            fin.close();
        }
    }

#ifdef NEW_ENVELOPE
    if (!epsr_tags.empty()) {
        std::ifstream fin(epsr_tags);
        std::string   line;
        while (std::getline(fin, line)) {
            params.input_epsr_tags.push_back(std::stod(line));
        }
        fin.close();
    }
#endif

    /// set envelope
    igl::Timer               timer;
    GEO::Mesh                sf_mesh;
    json                     tree_with_ids;
    std::vector<std::string> meshes;
    if (!csg_file.empty()) {
        json          csg_tree = json({});
        std::ifstream file(csg_file);

        if (file.is_open())
            file >> csg_tree;
        else {
            logger().error("unable to open {} file", csg_file);
            return EXIT_FAILURE;
        }
        file.close();

        CSGTreeParser::get_meshes(csg_tree, meshes, tree_with_ids);

        if (!CSGTreeParser::load_and_merge(
              meshes, input_vertices, input_faces, sf_mesh, input_tags))
            return EXIT_FAILURE;

        // To disable the recent modification of using input for wn, use meshes.clear();
    }
    else {
#ifdef NEW_ENVELOPE
        if (!MeshIO::load_mesh(params.input_path,
                               input_vertices,
                               input_faces,
                               sf_mesh,
                               input_tags,
                               params.input_epsr_tags)) {
#else
        if (!MeshIO::load_mesh(
              params.input_path, input_vertices, input_faces, sf_mesh, input_tags)) {
#endif
            logger().error("Unable to load mesh at {}", params.input_path);
            MeshIO::write_mesh(output_mesh_name, mesh, false);
            return EXIT_FAILURE;
        }
        else if (input_vertices.empty() || input_faces.empty()) {
            MeshIO::write_mesh(output_mesh_name, mesh, false);
            return EXIT_FAILURE;
        }

        if (input_tags.size() != input_faces.size()) {
            input_tags.resize(input_faces.size());
            std::fill(input_tags.begin(), input_tags.end(), 0);
        }
    }
    AABBWrapper tree(sf_mesh);
    if (!params.init(tree.get_sf_diag())) {
        return EXIT_FAILURE;
    }

#ifdef NEW_ENVELOPE
    if (!epsr_tags.empty())
        tree.init_sf_tree(
          input_vertices, input_faces, params.input_epsr_tags, params.bbox_diag_length);
    else
        tree.init_sf_tree(input_vertices, input_faces, params.eps);
#endif

#ifdef LIBIGL_WITH_TETGEN
    if (run_tet_gen) {
        Eigen::MatrixXd tetgen_pts(input_vertices.size(), 3);
        Eigen::MatrixXi tetgen_faces(input_faces.size(), 3);

        for (size_t i = 0; i < input_vertices.size(); ++i) {
            tetgen_pts.row(i) = input_vertices[i].cast<double>();
        }

        for (size_t i = 0; i < input_faces.size(); ++i) {
            tetgen_faces.row(i) = input_faces[i];
        }

        std::stringstream buf;
        buf.precision(100);
        buf.setf(std::ios::fixed, std::ios::floatfield);
        buf << "Qpq2.0a"
            << params.ideal_edge_length * params.ideal_edge_length * params.ideal_edge_length *
                 sqrt(2.) / 12.;

        Eigen::MatrixXi tetgen_generated_tets;
        Eigen::MatrixXd tetgen_generated_points;
        Eigen::MatrixXi tetgen_generated_faces;

        timer.start();
        igl::copyleft::tetgen::tetrahedralize(tetgen_pts,
                                              tetgen_faces,
                                              buf.str(),
                                              tetgen_generated_points,
                                              tetgen_generated_tets,
                                              tetgen_generated_faces);
        timer.stop();
        logger().info("Tetgen time {}s", timer.getElapsedTimeInSec());
        stats().record(StateInfo::tetgen_id,
                       timer.getElapsedTimeInSec(),
                       tetgen_generated_points.rows(),
                       tetgen_generated_tets.rows(),
                       0,
                       0);
    }
#endif

    stats().record(StateInfo::init_id, 0, input_vertices.size(), input_faces.size(), -1, -1);

    timer.start();
    simplify(input_vertices, input_faces, input_tags, tree, params, skip_simplify);
    tree.init_b_mesh_and_tree(input_vertices, input_faces, mesh);
    logger().info("preprocessing {}s", timer.getElapsedTimeInSec());
    logger().info("");
    stats().record(StateInfo::preprocessing_id,
                   timer.getElapsedTimeInSec(),
                   input_vertices.size(),
                   input_faces.size(),
                   -1,
                   -1);
    if (params.log_level <= 1)
        output_component(input_vertices, input_faces, input_tags);

    timer.start();
    std::vector<bool> is_face_inserted(input_faces.size(), false);
    FloatTetDelaunay::tetrahedralize(input_vertices, input_faces, tree, mesh, is_face_inserted);
    logger().info("#v = {}", mesh.get_v_num());
    logger().info("#t = {}", mesh.get_t_num());
    logger().info("tetrahedralizing {}s", timer.getElapsedTimeInSec());
    logger().info("");
    stats().record(StateInfo::tetrahedralization_id,
                   timer.getElapsedTimeInSec(),
                   mesh.get_v_num(),
                   mesh.get_t_num(),
                   -1,
                   -1);

    timer.start();
    insert_triangles(input_vertices, input_faces, input_tags, mesh, is_face_inserted, tree, false);
    logger().info("cutting {}s", timer.getElapsedTimeInSec());
    logger().info("");
    stats().record(StateInfo::cutting_id,
                   timer.getElapsedTimeInSec(),
                   mesh.get_v_num(),
                   mesh.get_t_num(),
                   mesh.get_max_energy(),
                   mesh.get_avg_energy(),
                   std::count(is_face_inserted.begin(), is_face_inserted.end(), false));

    //    timer.start();
    ////    cutting(input_vertices, input_faces, mesh, is_face_inserted, tree);
    //    cutting(input_vertices, input_faces, input_tags, mesh, is_face_inserted, tree);
    //    logger().info("cutting {}s", timer.getElapsedTimeInSec());
    //    logger().info("");
    //    stats().record(StateInfo::cutting_id, timer.getElapsedTimeInSec(), mesh.get_v_num(),
    //    mesh.get_t_num(),
    //                                                   mesh.get_max_energy(),
    //                                                   mesh.get_avg_energy(),
    //                                                   std::count(is_face_inserted.begin(),
    //                                                   is_face_inserted.end(), false));

    timer.start();
    optimization(
      input_vertices, input_faces, input_tags, is_face_inserted, mesh, tree, {{1, 1, 1, 1}});
    logger().info("mesh optimization {}s", timer.getElapsedTimeInSec());
    logger().info("");
    stats().record(StateInfo::optimization_id,
                   timer.getElapsedTimeInSec(),
                   mesh.get_v_num(),
                   mesh.get_t_num(),
                   mesh.get_max_energy(),
                   mesh.get_avg_energy());

    timer.start();
    correct_tracked_surface_orientation(mesh, tree);
    logger().info("correct_tracked_surface_orientation done");

    if (export_raw) {
        Eigen::Matrix<Scalar, Eigen::Dynamic, 3> Vt;
        Eigen::Matrix<int, Eigen::Dynamic, 3>    Ft;

        if (!csg_file.empty()) {
            int max_id = CSGTreeParser::get_max_id(tree_with_ids);

            for (int i = 0; i <= max_id; ++i) {
                get_tracked_surface(mesh, Vt, Ft, i);
                igl::write_triangle_mesh(
                  params.output_path + "_" + params.postfix + "_" + std::to_string(i) + "_all.obj",
                  Vt,
                  Ft);
            }
        }
        else {
            get_tracked_surface(mesh, Vt, Ft);
            igl::write_triangle_mesh(
              params.output_path + "_" + params.postfix + "_all.obj", Vt, Ft);
        }
        MeshIO::write_mesh(params.output_path + "_" + params.postfix + "_all.msh", mesh, false);
    }

    if (!csg_file.empty())
        boolean_operation(mesh, tree_with_ids, meshes);
    else if (boolean_op >= 0)
        boolean_operation(mesh, boolean_op);
    else {
        if (params.smooth_open_boundary) {
            smooth_open_boundary(mesh, tree);
            for (auto& t : mesh.tets) {
                if (t.is_outside)
                    t.is_removed = true;
            }
        }
        else {
            if (!params.disable_filtering) {
                if (params.use_floodfill) {
                    filter_outside_floodfill(mesh);
                }
                else if (params.use_input_for_wn) {
                    filter_outside(mesh, input_vertices, input_faces);
                }
                else
                    filter_outside(mesh);
            }
        }
    }
    Eigen::MatrixXd V_sf;
    Eigen::MatrixXi F_sf;
    if (params.manifold_surface) {
        manifold_surface(mesh, V_sf, F_sf);
    }
    else {
        get_surface(mesh, V_sf, F_sf);
    }
    stats().record(StateInfo::wn_id,
                   timer.getElapsedTimeInSec(),
                   mesh.get_v_num(),
                   mesh.get_t_num(),
                   mesh.get_max_energy(),
                   mesh.get_avg_energy());
    logger().info("after winding number");
    logger().info("#v = {}", mesh.get_v_num());
    logger().info("#t = {}", mesh.get_t_num());
    logger().info("winding number {}s", timer.getElapsedTimeInSec());
    logger().info("");

    //    if (params.output_path.size() > 3
    //        && params.output_path.substr(params.output_path.size() - 3, params.output_path.size())
    //        == "msh") MeshIO::write_mesh(params.output_path, mesh, false);
    //    else if (params.output_path.size() > 4
    //             && params.output_path.substr(params.output_path.size() - 4,
    //             params.output_path.size()) == "mesh")
    //        MeshIO::write_mesh(params.output_path, mesh, false);
    //    else
    //        MeshIO::write_mesh(params.output_path + "_" + params.postfix + ".msh", mesh, false);

    // fortest
    std::vector<Scalar> colors;
    if (!nocolor) {
        colors.resize(mesh.tets.size(), -1);
        for (int i = 0; i < mesh.tets.size(); i++) {
            if (mesh.tets[i].is_removed)
                continue;
            colors[i] = mesh.tets[i].quality;
        }
    }
    // fortest
    MeshIO::write_mesh(output_mesh_name, mesh, false, colors, !nobinary, !csg_file.empty());
    igl::write_triangle_mesh(params.output_path + "_" + params.postfix + "_sf.obj", V_sf, F_sf);
    //    MeshIO::write_surface_mesh(params.output_path + "_" + params.postfix + "_sf.obj", mesh,
    //    false);

    std::ofstream fout(params.log_path + "_" + params.postfix + ".csv");
    if (fout.good())
        fout << stats();
    fout.close();

    if (!params.envelope_log.empty()) {
        std::ofstream fout(params.envelope_log);
        fout << envelope_log_csv;
        fout.close();
    }

    return EXIT_SUCCESS;
}

//#include <igl/readSTL.h>
//#include <igl/writeSTL.h>
//#include <igl/writeOFF.h>
// void connect_2_meshes(std::string m1, std::string m2, std::string m) {
//    Eigen::MatrixXd v1, v2, _;
//    Eigen::MatrixXi f1, f2;
//
//    igl::readSTL(m1, v1, f1, _);
//    igl::readSTL(m2, v2, f2, _);
//
//    MatrixXd V(v1.rows() + v2.rows(), v1.cols());
//    V << v1, v2;
//
//    int v1_rows = v1.rows();
//    for (int i = 0; i < f2.rows(); i++) {
//        for (int j = 0; j < 3; j++)
//            f2(i, j) += v1_rows;
//    }
//    MatrixXi F(f1.rows() + f2.rows(), f1.cols());
//    F << f1, f2;
//
////    igl::writeOFF(m+".off", V, F);
//    igl::writeSTL(m+".stl", V, F);
//    std::ofstream fout(m+"_tags.txt");
//    for (int i = 0; i < f1.rows(); i++)
//        fout << 1 << endl;
//    for (int i = 0; i < f2.rows(); i++)
//        fout << 2 << endl;
//    fout.close();
//
//    //pausee();
//}
//
//#include <igl/readMESH.h>
// void test_manifold(std::string& file_name){
//    Eigen::MatrixXd V;
//    Eigen::MatrixXi T, F;
//    igl::readMESH(file_name, V, T, F);
//
//    Mesh mesh;
//
//    Eigen::MatrixXd V_sf;
//    Eigen::MatrixXi F_sf;
//    manifold_surface(mesh, V_sf, F_sf);
//}
