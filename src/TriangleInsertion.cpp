//
// Created by Yixin Hu on 2019-08-27.
//

#include <floattetwild/TriangleInsertion.h>

#include <floattetwild/FloatTetCutting.h>
#include <floattetwild/FloatTetCuttingCheck.h>
#include <floattetwild/FloatTetCuttingParallel.h>
#include <floattetwild/LocalOperations.h>
#include <floattetwild/Predicates.hpp>
#include <floattetwild/MeshIO.hpp>
#include <floattetwild/auto_table.hpp>
#include <floattetwild/Logger.hpp>
#include <floattetwild/intersections.h>

#include <floattetwild/MeshImprovement.h>//fortest

#include <igl/writeSTL.h>
#include <igl/writeOFF.h>
#include <igl/Timer.h>

#ifdef FLOAT_TETWILD_USE_TBB
#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include <tbb/atomic.h>
#include <floattetwild/FloatTetCuttingParallel.h>
#include <tbb/concurrent_queue.h>
#endif

#include <bitset>
#include <numeric>
#include <unordered_map>

#define III -1

//fortest
double time_find_cutting_tets = 0;
double time_find_cutting_tets1 = 0;
double time_find_cutting_tets2 = 0;
double time_find_cutting_tets3 = 0;
double time_find_cutting_tets4 = 0;
double time_cut_mesh = 0;
double time_cut_mesh1 = 0;
double time_cut_mesh2 = 0;
double time_get_intersecting_edges_and_points = 0;
double time_subdivide_tets = 0;
double time_push_new_tets = 0;
double time_push_new_tets1 = 0;
double time_push_new_tets2 = 0;
double time_push_new_tets3 = 0;
int cnt_snapped = 0;

double old_time_find_cutting_tets = 0;
double old_time_cut_mesh = 0;
double old_time_get_intersecting_edges_and_points = 0;
double old_time_subdivide_tets = 0;
double old_time_push_new_tets = 0;

int global_i = 0;
std::vector<std::array<int, 3>> coplanar_fs;//fortest
//fortest

void floatTetWild::sort_input_faces(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
                                    const Mesh &mesh, std::vector<int> &sorted_f_ids) {///use 38416, 232368 as example //todo: think why
    std::vector<Scalar> weights(input_faces.size());
    sorted_f_ids.resize(input_faces.size());
    for (int i = 0; i < input_faces.size(); i++) {
        sorted_f_ids[i] = i;

//        //fortest
//        Vector3 u = input_vertices[input_faces[i][1]] - input_vertices[input_faces[i][0]];
//        Vector3 v = input_vertices[input_faces[i][2]] - input_vertices[input_faces[i][0]];
//        if(u.cross(v).norm()/2 < SCALAR_ZERO_2) {
//            cout << "degenerate input triangle!!" << endl;
//            pausee();
//        }
//        //fortest

//        for (int j = 0; j < 3; j++) {
//            Scalar dis =
//                    (input_vertices[input_faces[i][j]] - input_vertices[input_faces[i][(j + 1) % 3]]).squaredNorm();
//            if (j == 0)
//                weights[i] = dis;
//            else if (dis > weights[i])
//                weights[i] = dis;
//        }
        Vector3 u = input_vertices[input_faces[i][1]] - input_vertices[input_faces[i][0]];
        Vector3 v = input_vertices[input_faces[i][2]] - input_vertices[input_faces[i][0]];
        weights[i] = u.cross(v).squaredNorm();
    }

    if(mesh.params.not_sort_input)
        return;

    std::random_shuffle(sorted_f_ids.begin(), sorted_f_ids.end());
//    std::sort(sorted_f_ids.begin(), sorted_f_ids.end(), [&weights](int a, int b) {
//        return weights[a] < weights[b];
//    });
}

void floatTetWild::insert_triangles(const std::vector<Vector3> &input_vertices,
        const std::vector<Vector3i> &input_faces, const std::vector<int> &input_tags,
        Mesh &mesh, std::vector<bool> &is_face_inserted, AABBWrapper &tree, bool is_again) {

    logger().info("triangle insertion start, #f = {}, #v = {}, #t = {}",
                  input_faces.size(), mesh.tet_vertices.size(), mesh.tets.size());
    /////
    std::vector<std::array<std::vector<int>, 4 >> track_surface_fs(mesh.tets.size());
    if (!is_again) {
        match_surface_fs(mesh, input_vertices, input_faces, is_face_inserted, track_surface_fs);
    }
    int cnt_matched = std::count(is_face_inserted.begin(), is_face_inserted.end(), true);
    logger().info("matched #f = {}, uninserted #f = {}", cnt_matched, is_face_inserted.size() - cnt_matched);

    std::vector<int> sorted_f_ids;
    sort_input_faces(input_vertices, input_faces, mesh, sorted_f_ids);

    /////
    std::vector<Vector3> new_vertices;
    std::vector<std::array<int, 4>> new_tets;
    int cnt_fail = 0;
    int cnt_total = 0;
//    for (int i = 0; i < input_faces.size(); i++) {
//        int f_id = i;
    for (int i = 0; i < sorted_f_ids.size(); i++) {
//        if(i>30000)
//            cout << "fid " << sorted_f_ids[i] << endl;

        //fortest
        global_i = i;
        if (!is_again && i > 0 && i % 1000 == 0) {
            logger().info("inserting f{}... {} failed", i, cnt_fail);
            logger().info("snapped {}/{}", cnt_snapped, cnt_total);
            logger().info("\t- time_find_cutting_tets = {}s (total {}s)",
                          time_push_new_tets - old_time_push_new_tets, time_find_cutting_tets);
//            logger().info("\t\t- time_find_cutting_tets1 = {}s", time_find_cutting_tets1);
//            logger().info("\t\t- time_find_cutting_tets2 = {}s", time_find_cutting_tets2);
//            logger().info("\t\t- time_find_cutting_tets3 = {}s", time_find_cutting_tets3);
//            logger().info("\t\t- time_find_cutting_tets4 = {}s", time_find_cutting_tets4);
            logger().info("\t- time_cut_mesh = {}s (total {}s)",
                          time_cut_mesh - old_time_cut_mesh, time_cut_mesh);
//            logger().info("\t\t- time_cut_mesh1 = {}s", time_cut_mesh1);
//            logger().info("\t\t- time_cut_mesh2 = {}s", time_cut_mesh2);
//            print_times1();
            logger().info("\t- time_get_intersecting_edges_and_points = {}s (total {}s)",
                          time_get_intersecting_edges_and_points - old_time_get_intersecting_edges_and_points,
                          time_get_intersecting_edges_and_points);
            logger().info("\t- time_subdivide_tets = {}s (total {}s)",
                          time_subdivide_tets - old_time_subdivide_tets, time_subdivide_tets);
            logger().info("\t- time_push_new_tets = {}s (total {}s)",
                          time_push_new_tets - old_time_push_new_tets, time_push_new_tets);
//            logger().info("\t\t- time_push_new_tets1 = {}s", time_push_new_tets1);
//            logger().info("\t\t- time_push_new_tets2 = {}s", time_push_new_tets2);
//            logger().info("\t\t- time_push_new_tets3 = {}s", time_push_new_tets3);

            old_time_find_cutting_tets = time_find_cutting_tets;
            old_time_cut_mesh = time_cut_mesh;
            old_time_get_intersecting_edges_and_points = time_get_intersecting_edges_and_points;
            old_time_subdivide_tets = time_subdivide_tets;
            old_time_push_new_tets = time_push_new_tets;
        }
        //fortest

        int f_id = sorted_f_ids[i];
        if (is_face_inserted[f_id])
            continue;

        cnt_total++;
        if (insert_one_triangle(f_id, input_vertices, input_faces, input_tags, mesh, track_surface_fs, tree, is_again))
            is_face_inserted[f_id] = true;
        else
            cnt_fail++;

//        pausee();//fortest
        if (f_id == III)
            break;//fortest
    }
    logger().info("insert_one_triangle * n done, #v = {}, #t = {}", mesh.tet_vertices.size(), mesh.tets.size());
    logger().info("uninserted #f = {}/{}", std::count(is_face_inserted.begin(), is_face_inserted.end(), false),
                  is_face_inserted.size() - cnt_matched);
    logger().info("total timing: {}s", time_find_cutting_tets + time_cut_mesh + time_get_intersecting_edges_and_points +
                                       time_subdivide_tets + time_push_new_tets);

    pair_track_surface_fs(mesh, track_surface_fs);
    logger().info("pair_track_surface_fs done");


    /////
    std::vector<std::pair<std::array<int, 2>, std::vector<int>>> b_edge_infos;
    find_boundary_edges(input_vertices, input_faces, is_face_inserted, b_edge_infos);
    logger().info("find_boundary_edges done");
    std::vector<std::array<int, 2>> b_edges1;
    std::vector<std::array<int, 3>> known_surface_fs;
    std::vector<std::array<int, 3>> known_not_surface_fs;
    mark_boundary_vs(input_vertices, input_faces, track_surface_fs, b_edge_infos, mesh, b_edges1);

//    insert_boundary_edges(input_vertices, input_faces, b_edge_infos, track_surface_fs, mesh, tree, b_edges1,
//                          is_face_inserted, is_again, known_surface_fs, known_not_surface_fs);
//    logger().info("uninserted #f = {}/{}", std::count(is_face_inserted.begin(), is_face_inserted.end(), false),
//                  is_face_inserted.size() - cnt_matched);

    //fortest
    check_track_surface_fs(mesh, track_surface_fs);
    //fortest

    /////
    std::vector<std::array<int, 2>> b_edges2;
    mark_surface_fs(input_vertices, input_faces, track_surface_fs, is_face_inserted,
                    known_surface_fs, known_not_surface_fs, b_edges2, mesh, tree);
    logger().info("mark_surface_fs done");

    /////
    //build b_tree using b_edges
    tree.init_tmp_b_mesh_and_tree(input_vertices, input_faces, b_edges1, mesh, b_edges2);
    if (!is_again) {
        for (auto &t:mesh.tets)
            t.quality = get_quality(mesh, t);
    }
    if(std::count(is_face_inserted.begin(), is_face_inserted.end(), false) == 0) {
        mesh.is_input_all_inserted = true;
    }
    logger().info("#b_edge1 = {}, #b_edges2 = {}", b_edges1.size(), b_edges2.size());

    ///fortest
    Eigen::MatrixXd V(input_vertices.size(), 3);
    Eigen::MatrixXi F(std::count(is_face_inserted.begin(), is_face_inserted.end(), true), 3);
    for (int i = 0; i < input_vertices.size(); i++)
        V.row(i) = input_vertices[i];
    int cnt = 0;
    for (int i = 0; i < input_faces.size(); i++) {
        if (!is_face_inserted[i])
            continue;
        F.row(cnt) << input_faces[i][0], input_faces[i][1], input_faces[i][2];
        cnt++;
    }
    igl::writeSTL("inserted.stl", V, F);
    ///fortest
}

bool floatTetWild::insert_one_triangle(int insert_f_id, const std::vector<Vector3> &input_vertices,
        const std::vector<Vector3i> &input_faces, const std::vector<int> &input_tags,
        Mesh &mesh, std::vector<std::array<std::vector<int>, 4>>& track_surface_fs,
        AABBWrapper &tree, bool is_again) {

    igl::Timer timer;
    std::array<Vector3, 3> vs = {{input_vertices[input_faces[insert_f_id][0]],
                                         input_vertices[input_faces[insert_f_id][1]],
                                         input_vertices[input_faces[insert_f_id][2]]}};
    Vector3 n = (vs[1] - vs[0]).cross(vs[2] - vs[0]);
    n.normalize();
    int t = get_t(vs[0], vs[1], vs[2]);

    /////
    timer.start();
    std::vector<int> cut_t_ids;
    find_cutting_tets(insert_f_id, input_vertices, input_faces, vs, mesh, cut_t_ids, is_again);
    time_find_cutting_tets += timer.getElapsedTime();

    //fortest
    myassert(!cut_t_ids.empty(), "cut_t_ids.empty()!!!");
    if (cut_t_ids.empty()) {
        cout << get_area(vs[0], vs[1], vs[2]) << endl;
        cout << "f" << insert_f_id << ": " << input_faces[insert_f_id][0] << " " << input_faces[insert_f_id][1]
             << " " << input_faces[insert_f_id][2] << endl;
        pausee();
        return false;
    }
    //fortest

    /////
    timer.start();
    igl::Timer timer1;
    timer1.start();
    CutMesh cut_mesh(mesh, n, vs);
    cut_mesh.construct(cut_t_ids);
    time_cut_mesh1 += timer1.getElapsedTime();
    timer1.start();
//    bool is_expanded = false;//fortest

    if (cut_mesh.snap_to_plane()) {
//        //fortest
//        if(!cut_mesh.check()) {
//            cout<<"checking cut_mesh 1"<<endl;
//            pausee();
//        }
//        //fortest

        cnt_snapped++;
//        cout<<cut_t_ids.size()<<"->";
//        cut_mesh.expand(cut_t_ids);
        cut_mesh.expand_new(cut_t_ids);
//        cout<<cut_t_ids.size()<<endl;
//        vector_print(cut_t_ids);
//        pausee();
    }
//    //fortest
//    if(!cut_mesh.check()) {
//        cout<<"checking cut_mesh 2"<<endl;
//        pausee();
//    }
    //fortest
    time_cut_mesh2 += timer1.getElapsedTime();
    time_cut_mesh += timer.getElapsedTime();

    /////
    timer.start();
    std::vector<Vector3> points;
    std::map<std::array<int, 2>, int> map_edge_to_intersecting_point;
    std::vector<int> subdivide_t_ids;
    if (!cut_mesh.get_intersecting_edges_and_points(points, map_edge_to_intersecting_point, subdivide_t_ids)) {
        time_get_intersecting_edges_and_points += timer.getElapsedTime();
//        cout<<"FAIL get_intersecting_edges_and_points"<<endl;
//        if(is_expanded)
//            cout<<"expanded"<<endl;
        return false;
    }
    //have to add all cut_t_ids
    vector_unique(cut_t_ids);
    std::vector<int> tmp;
    std::set_difference(subdivide_t_ids.begin(), subdivide_t_ids.end(), cut_t_ids.begin(), cut_t_ids.end(),
                        std::back_inserter(tmp));
    std::vector<bool> is_mark_surface(cut_t_ids.size(), true);
    cut_t_ids.insert(cut_t_ids.end(), tmp.begin(), tmp.end());
    is_mark_surface.resize(is_mark_surface.size() + tmp.size(), false);
//    cout << "cut_mesh.get_intersecting_edges_and_points OK" << endl;
    time_get_intersecting_edges_and_points += timer.getElapsedTime();

    /////
    timer.start();
    std::vector<MeshTet> new_tets;
    std::vector<std::array<std::vector<int>, 4>> new_track_surface_fs;
    std::vector<int> modified_t_ids;
    std::vector<std::pair<int, int>> covered_tet_fs;
    if (!subdivide_tets(insert_f_id, mesh, cut_mesh, points, map_edge_to_intersecting_point, track_surface_fs,
                        cut_t_ids, is_mark_surface,
                        new_tets, new_track_surface_fs, modified_t_ids,
                        covered_tet_fs)) {
        time_subdivide_tets += timer.getElapsedTime();
//        cout<<"FAIL subdivide_tets"<<endl;
        return false;
    }
//    cout<<covered_tet_fs.size()<<endl;
//    cout<<coplanar_fs.size()<<endl;
//    pausee();
    time_subdivide_tets += timer.getElapsedTime();

    timer.start();
    std::map<int, int> map_t_ids = push_new_tets(mesh, track_surface_fs, points, new_tets, new_track_surface_fs,
                                                 modified_t_ids, is_again);
    for (auto &f: covered_tet_fs)
        f.first = map_t_ids[f.first];
    time_push_new_tets += timer.getElapsedTime();

    check_is_input_face_covered(vs, mesh, covered_tet_fs);//fortest

    //preserve edges
    for (int i = 0; i < 3; i++) {
//        cout << "covered_tet_fs.size = " << covered_tet_fs.size() << endl;
        if (!preserve_edges(mesh, track_surface_fs, insert_f_id, vs, n, t, i, covered_tet_fs))
            return false;
    }
    //track surface faces
    std::array<Vector2, 3> vs_2d = {{to_2d(vs[0], t), to_2d(vs[1], t), to_2d(vs[2], t)}};
    std::vector<std::pair<int, int>> tmp_tet_fs;//fortest
    for (const auto &info: covered_tet_fs) {
        int t_id = info.first;
        int j = info.second;
        Vector3 c = mesh.tet_vertices[mesh.tets[t_id][(j + 1) % 4]].pos
                    + mesh.tet_vertices[mesh.tets[t_id][(j + 2) % 4]].pos
                    + mesh.tet_vertices[mesh.tets[t_id][(j + 3) % 4]].pos;
        c /= 3;
        Vector2 c_2d = to_2d(c, t);
        if (!is_p_inside_tri_2d(c_2d, vs_2d)) {
            if(!vector_erase(track_surface_fs[t_id][j], insert_f_id))
                cout<<"!vector_erase(track_surface_fs[t_id][j], insert_f_id)"<<endl;
        } else
            tmp_tet_fs.push_back(info);
    }


//    //fortest
//    if (global_i > 6300) {
//        cout << global_i << " covered" << endl;
//        //output input/tet triangles in different colors
//        {
//            Eigen::MatrixXd V(tmp_tet_fs.size() * 3, 3), C(tmp_tet_fs.size() * 3, 3);
//            Eigen::MatrixXi F(tmp_tet_fs.size(), 3);
//            for (int i = 0; i < tmp_tet_fs.size(); i++) {
//                int t_id = tmp_tet_fs[i].first;
//                int j = tmp_tet_fs[i].second;
//                std::array<int, 3> f = {{mesh.tets[t_id][(j + 1) % 4], mesh.tets[t_id][(j + 2) % 4],
//                                         mesh.tets[t_id][(j + 3) % 4]}};
//
//                for (int j = 0; j < 3; j++) {
//                    V.row(i * 3 + j) = mesh.tet_vertices[f[j]].pos;
//                    C.row(i * 3 + j) << 0, 0, 255;
//                }
//                F.row(i) << i * 3, i * 3 + 1, i * 3 + 2;
//            }
//            igl::writeOFF("covered_tet_fs_" + std::to_string(global_i) + ".off", V, F, C);
//        }
//        {
//            Eigen::MatrixXd V(1 * 3, 3), C(1 * 3, 3);
//            Eigen::MatrixXi F(1, 3);
//            int i = 0;
//            for (int j = 0; j < 3; j++) {
//                V.row(i * 3 + j) = vs[j];
//                C.row(i * 3 + j) << 255, 0, 0;
//            }
//            F.row(i) << i * 3, i * 3 + 1, i * 3 + 2;
//            igl::writeOFF("covered_input_f_" + std::to_string(global_i) + ".off", V, F, C);
//        }
//    }
//    //fortest

    return true;
}

std::map<int, int> floatTetWild::push_new_tets(Mesh &mesh, std::vector<std::array<std::vector<int>, 4>> &track_surface_fs,
                                 std::vector<Vector3> &points, std::vector<MeshTet> &new_tets,
                                 std::vector<std::array<std::vector<int>, 4>> &new_track_surface_fs,
                                 std::vector<int> &modified_t_ids, bool is_again) {
//    igl::Timer timer;
//    timer.start();
    ///vs
    const int old_v_size = mesh.tet_vertices.size();
    mesh.tet_vertices.resize(mesh.tet_vertices.size() + points.size());
    for (int i = 0; i < points.size(); i++) {
        mesh.tet_vertices[old_v_size + i].pos = points[i];
        //todo: tags???
    }
//    time_push_new_tets1 += timer.getElapsedTime();

    ///tets
//    timer.start();
//    mesh.tets.reserve(mesh.tets.size() + new_tets.size() - modified_t_ids.size());
//    time_push_new_tets2 += timer.getElapsedTime();

//    timer.start();
    std::map<int, int> map_t_ids;
    for (int i = 0; i < new_tets.size(); i++) {
        if (is_again)
            new_tets[i].quality = get_quality(mesh, new_tets[i]);

        if (i < modified_t_ids.size()) {
            for (int j = 0; j < 4; j++) {
                vector_erase(mesh.tet_vertices[mesh.tets[modified_t_ids[i]][j]].conn_tets, modified_t_ids[i]);
            }
            mesh.tets[modified_t_ids[i]] = new_tets[i];
            track_surface_fs[modified_t_ids[i]] = new_track_surface_fs[i];
            for (int j = 0; j < 4; j++) {
                mesh.tet_vertices[mesh.tets[modified_t_ids[i]][j]].conn_tets.push_back(modified_t_ids[i]);
            }
            map_t_ids[i] = modified_t_ids[i];
        } else {
//            mesh.tets.push_back(new_tets[i]);
//            track_surface_fs.push_back(new_track_surface_fs[i]);
//            for (int j = 0; j < 4; j++) {
//                mesh.tet_vertices[mesh.tets.back()[j]].conn_tets.push_back(mesh.tets.size() - 1);
//            }
            for (int j = 0; j < 4; j++) {
                mesh.tet_vertices[new_tets[i][j]].conn_tets.push_back(mesh.tets.size() + i - modified_t_ids.size());
            }
            map_t_ids[i] = mesh.tets.size() + i - modified_t_ids.size();
        }
        //todo: tags???
    }
    mesh.tets.insert(mesh.tets.end(), new_tets.begin() + modified_t_ids.size(), new_tets.end());
    track_surface_fs.insert(track_surface_fs.end(), new_track_surface_fs.begin() + modified_t_ids.size(),
                            new_track_surface_fs.end());
//    time_push_new_tets3 += timer.getElapsedTime();

    return map_t_ids;
}

void floatTetWild::find_cutting_tets(int f_id, const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
                                     const std::array<Vector3, 3>& vs, Mesh &mesh, std::vector<int> &cut_t_ids, bool is_again) {
    //todo: double check why different
//    std::vector<int> all_cut_t_ids;
//    for (int t_id = 0; t_id < mesh.tets.size(); t_id++) {
//        std::array<int, 4> oris;
//        for (int j = 0; j < 4; j++) {
//            oris[j] = Predicates::orient_3d(vs[0], vs[1], vs[2], mesh.tet_vertices[mesh.tets[t_id][j]].pos);
//        }
//
//        for (int j = 0; j < 4; j++) {
//            int cnt_pos = 0;
//            int cnt_neg = 0;
//            int cnt_on = 0;
//            for (int k = 0; k < 3; k++) {
//                if (oris[(j + k + 1) % 4] == Predicates::ORI_ZERO)
//                    cnt_on++;
//                else if (oris[(j + k + 1) % 4] == Predicates::ORI_POSITIVE)
//                    cnt_pos++;
//                else
//                    cnt_neg++;
//            }
//
//            int result = CUT_EMPTY;
//            auto &tp1 = mesh.tet_vertices[mesh.tets[t_id][(j + 1) % 4]].pos;
//            auto &tp2 = mesh.tet_vertices[mesh.tets[t_id][(j + 2) % 4]].pos;
//            auto &tp3 = mesh.tet_vertices[mesh.tets[t_id][(j + 3) % 4]].pos;
//            if (cnt_on == 3) {
//                result = is_tri_tri_cutted_hint(vs[0], vs[1], vs[2], tp1, tp2, tp3, CUT_COPLANAR);
//            } else if (cnt_pos > 0 && cnt_neg > 0) {
//                result = is_tri_tri_cutted_hint(vs[0], vs[1], vs[2], tp1, tp2, tp3, CUT_FACE);
//            }
//            if (result == CUT_EMPTY)
//                continue;
//
//            all_cut_t_ids.push_back(t_id);
//            break;
//        }
//    }
//
//    cut_t_ids = all_cut_t_ids;
//    return;//fortest

    std::vector<bool> is_visited(mesh.tets.size(), false);
    std::queue<int> queue_t_ids;

    if (!is_again) {
        std::vector<int> n_t_ids;
        for (int j = 0; j < 3; j++) {
            n_t_ids.insert(n_t_ids.end(), mesh.tet_vertices[input_faces[f_id][j]].conn_tets.begin(),
                           mesh.tet_vertices[input_faces[f_id][j]].conn_tets.end());
        }
        vector_unique(n_t_ids);

        for (int t_id: n_t_ids) {
            is_visited[t_id] = true;
            queue_t_ids.push(t_id);
        }
    } else {
        Vector3 min_f, max_f;
        get_bbox_face(input_vertices[input_faces[f_id][0]], input_vertices[input_faces[f_id][1]],
                      input_vertices[input_faces[f_id][2]], min_f, max_f);

        for (int t_id = 0; t_id < mesh.tets.size(); t_id++) {
            if (mesh.tets[t_id].is_removed)
                continue;

            Vector3 min_t, max_t;
            get_bbox_tet(mesh.tet_vertices[mesh.tets[t_id][0]].pos, mesh.tet_vertices[mesh.tets[t_id][1]].pos,
                         mesh.tet_vertices[mesh.tets[t_id][2]].pos, mesh.tet_vertices[mesh.tets[t_id][3]].pos,
                         min_t, max_t);
            if (!is_bbox_intersected(min_f, max_f, min_t, max_t))
                continue;

            queue_t_ids.push(t_id);
            is_visited[t_id] = true;
        }
    }

    while (!queue_t_ids.empty()) {
        int t_id = queue_t_ids.front();
        queue_t_ids.pop();
        std::array<int, 4> oris;
        int cnt_pos = 0;
        int cnt_neg = 0;
        int cnt_on = 0;
        for (int j = 0; j < 4; j++) {
            oris[j] = Predicates::orient_3d(vs[0], vs[1], vs[2], mesh.tet_vertices[mesh.tets[t_id][j]].pos);
            if (oris[j] == Predicates::ORI_ZERO)
                cnt_on++;
            else if (oris[j] == Predicates::ORI_POSITIVE)
                cnt_pos++;
            else
                cnt_neg++;
        }
        if (cnt_pos == 0 && cnt_neg == 0 && cnt_on < 3)
            continue;

        bool is_cutted = false;
        std::vector<bool> is_cut_vs = {{false, false, false, false}}; /// is v on cut face
        for (int j = 0; j < 4; j++) {
            int cnt_pos = 0;
            int cnt_neg = 0;
            int cnt_on = 0;
            for (int k = 0; k < 3; k++) {
                if (oris[(j + k + 1) % 4] == Predicates::ORI_ZERO)
                    cnt_on++;
                else if (oris[(j + k + 1) % 4] == Predicates::ORI_POSITIVE)
                    cnt_pos++;
                else
                    cnt_neg++;
            }

            int result = CUT_EMPTY;
            auto &tp1 = mesh.tet_vertices[mesh.tets[t_id][(j + 1) % 4]].pos;
            auto &tp2 = mesh.tet_vertices[mesh.tets[t_id][(j + 2) % 4]].pos;
            auto &tp3 = mesh.tet_vertices[mesh.tets[t_id][(j + 3) % 4]].pos;
            if (cnt_on == 3) {
                if (is_tri_tri_cutted_hint(vs[0], vs[1], vs[2], tp1, tp2, tp3, CUT_COPLANAR) == CUT_COPLANAR) {
                    is_cutted = true;
                    is_cut_vs[(j + 1) % 4] = true;
                    is_cut_vs[(j + 2) % 4] = true;
                    is_cut_vs[(j + 3) % 4] = true;
                }
            } else if (cnt_pos > 0 && cnt_neg > 0) {
                if (is_tri_tri_cutted_hint(vs[0], vs[1], vs[2], tp1, tp2, tp3, CUT_FACE) == CUT_FACE) {
                    is_cutted = true;
                    is_cut_vs[(j + 1) % 4] = true;
                    is_cut_vs[(j + 2) % 4] = true;
                    is_cut_vs[(j + 3) % 4] = true;
                }
            } else if (cnt_on == 2 && oris[(j + 1) % 4] == Predicates::ORI_ZERO
                       && oris[(j + 2) % 4] == Predicates::ORI_ZERO) {
                if (is_tri_tri_cutted_hint(vs[0], vs[1], vs[2], tp1, tp2, tp3, CUT_EDGE_0) == CUT_EDGE_0) {
                    is_cut_vs[(j + 1) % 4] = true;
                    is_cut_vs[(j + 2) % 4] = true;
                }
            } else if (cnt_on == 2 && oris[(j + 2) % 4] == Predicates::ORI_ZERO
                       && oris[(j + 3) % 4] == Predicates::ORI_ZERO) {
                if (is_tri_tri_cutted_hint(vs[0], vs[1], vs[2], tp1, tp2, tp3, CUT_EDGE_1) == CUT_EDGE_1) {
                    is_cut_vs[(j + 2) % 4] = true;
                    is_cut_vs[(j + 3) % 4] = true;
                }
            } else if (cnt_on == 2 && oris[(j + 3) % 4] == Predicates::ORI_ZERO
                       && oris[(j + 1) % 4] == Predicates::ORI_ZERO) {
                if (is_tri_tri_cutted_hint(vs[0], vs[1], vs[2], tp1, tp2, tp3, CUT_EDGE_2) == CUT_EDGE_2) {
                    is_cut_vs[(j + 3) % 4] = true;
                    is_cut_vs[(j + 1) % 4] = true;
                }
            }
            if (is_cut_vs[0] && is_cut_vs[1] && is_cut_vs[2] && is_cut_vs[3])
                break;
        }
        if (is_cutted)
            cut_t_ids.push_back(t_id);

        for (int j = 0; j < 4; j++) {
            if (!is_cut_vs[j])
                continue;
            for (int n_t_id: mesh.tet_vertices[mesh.tets[t_id][j]].conn_tets) {
                if (!is_visited[n_t_id]) {
                    is_visited[n_t_id] = true;
                    queue_t_ids.push(n_t_id);
                }
            }
        }
    }

//    //fortest
//    std::sort(cut_t_ids.begin(), cut_t_ids.end());
//    if(cut_t_ids == all_cut_t_ids)
//        return;
//    std::vector<int> tmp;
//    std::set_difference(all_cut_t_ids.begin(), all_cut_t_ids.end(), cut_t_ids.begin(), cut_t_ids.end(),
//                        std::back_inserter(tmp));
//    cout<<"diff: ";
//    vector_print(tmp);
//    pausee();
//    //fortest
}

bool floatTetWild::subdivide_tets_edge(int insert_f_id, Mesh &mesh, std::vector<Vector3> &points,
                    std::map<std::array<int, 2>, int> &map_edge_to_intersecting_point,
                    std::vector<std::array<std::vector<int>, 4>> &track_surface_fs,
                    std::vector<int> &subdivide_t_ids,
                    std::vector<MeshTet> &new_tets,
                    std::vector<std::array<std::vector<int>, 4>> &new_track_surface_fs,
                    std::vector<int> &modified_t_ids,
                    std::vector<std::pair<int, int>>& new_covered_tet_fs) {
    std::vector<bool> is_mark_surface(subdivide_t_ids.size(), false);
    Vector3 n;
    std::array<Vector3, 3> vs;
    CutMesh cut_mesh(mesh, n, vs);
    return subdivide_tets(insert_f_id, mesh, cut_mesh, points, map_edge_to_intersecting_point, track_surface_fs,
                          subdivide_t_ids, is_mark_surface, new_tets, new_track_surface_fs, modified_t_ids,
                          new_covered_tet_fs, true);
}

bool floatTetWild::subdivide_tets(int insert_f_id, Mesh& mesh, CutMesh& cut_mesh,
        std::vector<Vector3>& points, std::map<std::array<int, 2>, int>& map_edge_to_intersecting_point,
        std::vector<std::array<std::vector<int>, 4>>& track_surface_fs,
        std::vector<int>& subdivide_t_ids, std::vector<bool>& is_mark_surface,
        std::vector<MeshTet>& new_tets, std::vector<std::array<std::vector<int>, 4>>& new_track_surface_fs,
        std::vector<int>& modified_t_ids,
        std::vector<std::pair<int, int>>& new_covered_tet_fs,
        bool is_track_covered_tet_fs) {

    static const std::array<std::array<int, 2>, 6> t_es = {{{{0, 1}}, {{1, 2}}, {{2, 0}}, {{0, 3}}, {{1, 3}}, {{2, 3}}}};
    static const std::array<std::array<int, 3>, 4> t_f_es = {{{{1, 5, 4}}, {{5, 3, 2}}, {{3, 0, 4}}, {{0, 1, 2}}}};
    static const std::array<std::array<int, 3>, 4> t_f_vs = {{{{3, 1, 2}}, {{0, 2, 3}}, {{1, 3, 0}}, {{2, 0, 1}}}};

    //fortest
//    for (auto m: map_edge_to_intersecting_point)
//        cout << (m.first[0]) << " " << (m.first[1]) << ": " << mesh.tet_vertices.size() + m.second << endl;
    //fortest

    coplanar_fs.clear();
    for (int I = 0; I < subdivide_t_ids.size(); I++) {
        int t_id = subdivide_t_ids[I];
        bool is_mark_sf = is_mark_surface[I];

        //fortest
//        cout << endl << "t_id = " << t_id << endl;
//        cout << "is_mark_sf = " << is_mark_sf << endl;
//        cout << mesh.tets[t_id][0] << " " << mesh.tets[t_id][1] << " " << mesh.tets[t_id][2] << " "
//             << mesh.tets[t_id][3] << endl;
        //fortest

        /////
        std::bitset<6> config_bit;
        std::array<std::pair<int, int>, 6> on_edge_p_ids;
        int cnt = 4;
        for (int i = 0; i < t_es.size(); i++) {
            const auto &le = t_es[i];
            std::array<int, 2> e = {{mesh.tets[t_id][le[0]], mesh.tets[t_id][le[1]]}};
            if (e[0] > e[1])
                std::swap(e[0], e[1]);
            if (map_edge_to_intersecting_point.find(e) == map_edge_to_intersecting_point.end()) {
                on_edge_p_ids[i].first = -1;
                on_edge_p_ids[i].second = -1;
            } else {
                on_edge_p_ids[i].first = cnt++;
                on_edge_p_ids[i].second = map_edge_to_intersecting_point[e];
                config_bit.set(i);
            }
        }
        int config_id = config_bit.to_ulong();
//        cout << "config_id = " << config_id << endl;//fortest
        if (config_id == 0) { //no intersection
            if (is_mark_sf) {
//                cout << "is_mark_surface" << endl;//fortest
                for (int j = 0; j < 4; j++) {
                    int cnt_on = 0;
                    for (int k = 0; k < 3; k++) {
                        myassert(cut_mesh.map_v_ids.find(mesh.tets[t_id][(j + k + 1) % 4]) != cut_mesh.map_v_ids.end(),
                                "cut_mesh.map_v_ids.find(mesh.tets[t_id][(j + k + 1) % 4]) != cut_mesh.map_v_ids.end()!!");//fortest
                        if (cut_mesh.is_v_on_plane(cut_mesh.map_v_ids[mesh.tets[t_id][(j + k + 1) % 4]])) {
                            cnt_on++;
//                            cout << (j + k + 1) % 4 << endl;//fortest
                        }
                    }
//                    cout << "cnt_on = " << cnt_on << endl;//fortest
                    //fortest
                    if(cnt_on == 4){
                        cout<<"cnt_on==4!!"<<endl;
                    }
                    //fortest

                    if (cnt_on == 3) {
                        new_tets.push_back(mesh.tets[t_id]);
                        new_track_surface_fs.push_back(track_surface_fs[t_id]);
                        (new_track_surface_fs.back())[j].push_back(insert_f_id);
                        modified_t_ids.push_back(t_id);

                        coplanar_fs.push_back({{mesh.tets[t_id][(j+1)%4],mesh.tets[t_id][(j+2)%4],
                                                mesh.tets[t_id][(j+3)%4]}});//fortest

                        new_covered_tet_fs.push_back(std::make_pair(new_tets.size()-1, j));

//                        //fortest
//                        int opp_t_id = get_opp_t_id(t_id, j, mesh);
//                        if(std::find(subdivide_t_ids.begin(), subdivide_t_ids.end(), opp_t_id) == subdivide_t_ids.end()) {
//                            cout << "cannot find opp_t_id in subdivide_t_ids!" << endl;
//                            cout << "t_id j: " << t_id << " " << j << endl;
//                            mesh.tets[t_id].print();
//                            cout << "opp_t_id: " << opp_t_id << endl;
//                            mesh.tets[opp_t_id].print();
//                            for (int k = 0; k < 3; k++) {
//                                int lv_id = cut_mesh.map_v_ids[mesh.tets[t_id][(j + k + 1) % 4]];
//                                cout << lv_id << ": " << cut_mesh.to_plane_dists[lv_id] << " "
//                                     << cut_mesh.is_snapped[lv_id] << endl;
//                            }
//                            pausee();
//                        }
//                        //fortest

                        break;
                    }
                }
            }
            continue;
        }

        /////
        std::vector<Vector2i> my_diags;
        for (int j = 0; j < 4; j++) {
            std::vector<int> le_ids;
            for (int k = 0; k < 3; k++) {
                if (on_edge_p_ids[t_f_es[j][k]].first < 0)
                    continue;
                le_ids.push_back(k);
            }
            if (le_ids.size() != 2)//no ambiguity
                continue;

            my_diags.emplace_back();
            auto &diag = my_diags.back();
            if (on_edge_p_ids[t_f_es[j][le_ids[0]]].second > on_edge_p_ids[t_f_es[j][le_ids[1]]].second) {
                diag << on_edge_p_ids[t_f_es[j][le_ids[0]]].first, t_f_vs[j][le_ids[0]];
            } else {
                diag << on_edge_p_ids[t_f_es[j][le_ids[1]]].first, t_f_vs[j][le_ids[1]];
            }
            if (diag[0] > diag[1])
                std::swap(diag[0], diag[1]);
        }
        std::sort(my_diags.begin(), my_diags.end(), [](const Vector2i &a, const Vector2i &b) {
            return std::make_tuple(a[0], a[1]) < std::make_tuple(b[0], b[1]);
        });

        /////
        std::map<int, int> map_lv_to_v_id;
        const int v_size = mesh.tet_vertices.size();
        const int vp_size = mesh.tet_vertices.size() + points.size();
        for (int i = 0; i < 4; i++)
            map_lv_to_v_id[i] = mesh.tets[t_id][i];
        cnt = 0;
        for (int i = 0; i < t_es.size(); i++) {
            if (config_bit[i] == 0)
                continue;
            map_lv_to_v_id[4 + cnt] = v_size + on_edge_p_ids[i].second;
            cnt++;
        }

        /////
        auto get_centroid = [&](const std::vector<Vector4i> &config, int lv_id, Vector3 &c) {
            std::vector<int> n_ids;
            for (const auto &tet: config) {
                std::vector<int> tmp;
                for (int j = 0; j < 4; j++) {
                    if (tet[j] != lv_id)
                        tmp.push_back(tet[j]);
                }
                if (tmp.size() == 4)
                    continue;
                n_ids.insert(n_ids.end(), tmp.begin(), tmp.end());
            }
            vector_unique(n_ids);
            c << 0, 0, 0;
            for (int n_id:n_ids) {
                int v_id = map_lv_to_v_id[n_id];
                if (v_id < v_size)
                    c += mesh.tet_vertices[v_id].pos;
                else
                    c += points[v_id - v_size];
            }
            c /= n_ids.size();
        };

        auto check_config = [&](int diag_config_id, std::vector<std::pair<int, Vector3>> &centroids) {
            const std::vector<Vector4i> &config = CutTable::get_tet_conf(config_id, diag_config_id);
            Scalar min_q = -666;
            int cnt = 0;
            std::map<int, int> map_lv_to_c;
            for (const auto &tet: config) {
                std::array<Vector3, 4> vs;
                for (int j = 0; j < 4; j++) {
                    if (map_lv_to_v_id.find(tet[j]) == map_lv_to_v_id.end()) {
                        if (map_lv_to_c.find(tet[j]) == map_lv_to_c.end()) {
                            get_centroid(config, tet[j], vs[j]);
                            map_lv_to_c[tet[j]] = centroids.size();
                            centroids.push_back(std::make_pair(tet[j], vs[j]));
                        }
                        vs[j] = centroids[map_lv_to_c[tet[j]]].second;
                    } else {
                        int v_id = map_lv_to_v_id[tet[j]];
                        if (v_id < v_size)
                            vs[j] = mesh.tet_vertices[v_id].pos;
                        else
                            vs[j] = points[v_id - v_size];
                    }
                }

                Scalar volume = Predicates::orient_3d_volume(vs[0], vs[1], vs[2], vs[3]);

//                //fortest
//                if(volume==0) {
//                    cout<<std::setprecision(16)<<endl;
//                    cout << "volume = " << volume << ",ori = " << Predicates::orient_3d(vs[0], vs[1], vs[2], vs[3]) << endl;
//                    cout << "centroids.size = " << centroids.size() << endl;
//                    cout << "config_id = " << config_id << endl;
//
//                    cout<<"vertices"<<endl;
//                    for(int k=0;k<4;k++){
//                        cout<<"v"<<mesh.tets[t_id][k]<<": "<<mesh.tet_vertices[mesh.tets[t_id][k]].pos.transpose()<<endl;
//                    }
//                    cout<<"intersecting points"<<endl;
//                    for (int i = 0; i < t_es.size(); i++) {
//                        const auto &le = t_es[i];
//                        std::array<int, 2> e = {{mesh.tets[t_id][le[0]], mesh.tets[t_id][le[1]]}};
//                        if (e[0] > e[1])
//                            std::swap(e[0], e[1]);
//                        if (map_edge_to_intersecting_point.find(e) != map_edge_to_intersecting_point.end()) {
//                            cout << e[0] << " " << e[1] << ": p" << map_edge_to_intersecting_point[e] << " "
//                                 << points[map_edge_to_intersecting_point[e]].transpose() << endl;
//                            cout<<"is_mark_sf = "<<is_mark_sf<<endl;
//                            if(cut_mesh.map_v_ids.find(e[0])!=cut_mesh.map_v_ids.end()) {
//                                cout << "e[0] to_plane_dists = " << cut_mesh.to_plane_dists[cut_mesh.map_v_ids[e[0]]]
//                                     << endl;
//                                cout << "e[0] is_snapped = " << cut_mesh.is_snapped[cut_mesh.map_v_ids[e[0]]]
//                                     << endl;
//                            }
//                            if(cut_mesh.map_v_ids.find(e[1])!=cut_mesh.map_v_ids.end()) {
//                                cout << "e[1] to_plane_dists = " << cut_mesh.to_plane_dists[cut_mesh.map_v_ids[e[1]]]
//                                     << endl;
//                                cout << "e[1] is_snapped = " << cut_mesh.is_snapped[cut_mesh.map_v_ids[e[1]]]
//                                     << endl;
//                            }
//
////                            Vector3 p;
////                            Scalar _;
////                            seg_plane_intersection(mesh.tet_vertices[e[0]].pos, mesh.tet_vertices[e[1]].pos,
////                                                   cut_mesh.p_vs[0], cut_mesh.p_n, p, _);
////                            cout << "new p:" << p.transpose() << endl;
////                            cout << (p - mesh.tet_vertices[e[1]].pos).squaredNorm() << " " << SCALAR_ZERO_2 << endl;
//                        }
//                    }
//
//                    pausee();
//                }
//                //fortest

                if (cnt == 0)
                    min_q = volume;
                else if (volume < min_q)
                    min_q = volume;
                cnt++;
            }

            return min_q;
        };

        int diag_config_id = 0;
        std::vector<std::pair<int, Vector3>> centroids;
        if (!my_diags.empty()) {
            const auto &all_diags = CutTable::get_diag_confs(config_id);
            std::vector<std::pair<int, Scalar>> min_qualities(all_diags.size());
            std::vector<std::vector<std::pair<int, Vector3>>> all_centroids(all_diags.size());
            for (int i = 0; i < all_diags.size(); i++) {
                if (my_diags != all_diags[i]) {
                    min_qualities[i] = std::make_pair(i, -1);
                    continue;
                }

                std::vector<std::pair<int, Vector3>> tmp_centroids;
                Scalar min_q = check_config(i, tmp_centroids);
//                if (min_q < SCALAR_ZERO_3)
//                    continue;
                min_qualities[i] = std::make_pair(i, min_q);
                all_centroids[i] = tmp_centroids;
            }
            std::sort(min_qualities.begin(), min_qualities.end(),
                      [](const std::pair<int, Scalar> &a, const std::pair<int, Scalar> &b) {
                          return a.second < b.second;
                      });

            if (min_qualities.back().second < SCALAR_ZERO_3) { // if tet quality is too bad
//                cout<<std::setprecision(16)<<"return 1 "<<min_qualities.back().second<<endl;
                return false;
            }

            diag_config_id = min_qualities.back().first;
            centroids = all_centroids[diag_config_id];
        } else {
            Scalar min_q = check_config(diag_config_id, centroids);
            if (min_q < SCALAR_ZERO_3) {
//                cout<<std::setprecision(16)<<"return 2 "<<min_q<<endl;
                return false;
            }
        }

        for (int i = 0; i < centroids.size(); i++) {
            map_lv_to_v_id[centroids[i].first] = vp_size + i;
            points.push_back(centroids[i].second);
        }

        //add new tets
        const auto &config = CutTable::get_tet_conf(config_id, diag_config_id);
        const auto &new_is_surface_fs = CutTable::get_surface_conf(config_id, diag_config_id);
        const auto &new_local_f_ids = CutTable::get_face_id_conf(config_id, diag_config_id);
        for (int i = 0; i < config.size(); i++) {
            const auto &t = config[i];
            new_tets.push_back(MeshTet(map_lv_to_v_id[t[0]], map_lv_to_v_id[t[1]],
                                       map_lv_to_v_id[t[2]], map_lv_to_v_id[t[3]]));

            //fortest
//            cout << map_lv_to_v_id[t[0]] << " " << map_lv_to_v_id[t[1]] << " " << map_lv_to_v_id[t[2]] << " "
//                 << map_lv_to_v_id[t[3]] << endl;

            new_track_surface_fs.emplace_back();
            for (int j = 0; j < 4; j++) {
                if (new_is_surface_fs[i][j] && is_mark_sf) {
                    (new_track_surface_fs.back())[j].push_back(insert_f_id);

                    //fortest
                    coplanar_fs.push_back(
                            {{new_tets.back()[(j + 1) % 4], new_tets.back()[(j + 3) % 4], new_tets.back()[(j + 2) %
                                                                                                          4]}});
                    //fortest

                    new_covered_tet_fs.push_back(std::make_pair(new_tets.size() - 1, j));

//                    //fortest
//                    cout << "sf " << i << " " << j << endl;
//                    for (int k = 0; k < 3; k++) {
//                        int v_id = new_tets.back()[(j + k + 1) % 4];
//                        Vector3 p;
//                        if (v_id < v_size)
//                            p = mesh.tet_vertices[v_id].pos;
//                        else
//                            p = points[v_id - v_size];
//                        double dist = cut_mesh.get_to_plane_dist(p);
//                        cout << "v_id = " << v_id << endl;
//                        if (dist > mesh.params.eps_2_coplanar) {
//                            cout << "dist > mesh.params.eps_2_coplanar " << dist << endl;
//                            pausee();
//                        }
//                    }
//                    //fortest
                }

                int old_local_f_id = new_local_f_ids[i][j];
                if (old_local_f_id < 0)
                    continue;
                (new_track_surface_fs.back())[j].insert((new_track_surface_fs.back())[j].end(),
                                                        track_surface_fs[t_id][old_local_f_id].begin(),
                                                        track_surface_fs[t_id][old_local_f_id].end());
                (new_tets.back()).is_bbox_fs[j] = mesh.tets[t_id].is_bbox_fs[old_local_f_id];
                (new_tets.back()).is_surface_fs[j] = mesh.tets[t_id].is_surface_fs[old_local_f_id];
                (new_tets.back()).surface_tags[j] = mesh.tets[t_id].surface_tags[old_local_f_id];

                if (is_track_covered_tet_fs) {
                    if (vector_contains((new_track_surface_fs.back())[j], insert_f_id))
                        new_covered_tet_fs.push_back(std::make_pair(new_tets.size() - 1, j));
                }
            }
        }
        modified_t_ids.push_back(t_id);
    }

    return true;
}

void floatTetWild::pair_track_surface_fs(Mesh &mesh, std::vector<std::array<std::vector<int>, 4>> &track_surface_fs) {
    std::vector<std::array<bool, 4>> is_visited(track_surface_fs.size(), {{false, false, false, false}});
    for (int t_id = 0; t_id < track_surface_fs.size(); t_id++) {
        for (int j = 0; j < 4; j++) {
            //
            if (is_visited[t_id][j])
                continue;
            is_visited[t_id][j] = true;
            if (track_surface_fs[t_id][j].empty())
                continue;
            //
            int opp_t_id = get_opp_t_id(t_id, j, mesh);
            if (opp_t_id < 0)
                continue;
            int k = get_local_f_id(opp_t_id, mesh.tets[t_id][(j + 1) % 4], mesh.tets[t_id][(j + 2) % 4],
                                   mesh.tets[t_id][(j + 3) % 4], mesh);
            is_visited[opp_t_id][k] = true;
            //
            std::sort(track_surface_fs[t_id][j].begin(), track_surface_fs[t_id][j].end());
            std::sort(track_surface_fs[opp_t_id][k].begin(), track_surface_fs[opp_t_id][k].end());
            std::vector<int> f_ids;
            if (track_surface_fs[t_id][j] != track_surface_fs[opp_t_id][k]) {
                std::set_union(track_surface_fs[t_id][j].begin(), track_surface_fs[t_id][j].end(),
                               track_surface_fs[opp_t_id][k].begin(), track_surface_fs[opp_t_id][k].end(),
                               std::back_inserter(f_ids));
                track_surface_fs[t_id][j] = f_ids;
                track_surface_fs[opp_t_id][k] = f_ids;
            }
        }
    }
}

void floatTetWild::find_boundary_edges(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
                                       const std::vector<bool> &is_face_inserted,
                                       std::vector<std::pair<std::array<int, 2>, std::vector<int>>>& b_edge_infos) {
    std::vector<std::array<int, 2>> edges;
    std::vector<std::vector<int>> conn_tris(input_vertices.size());
    for (int i = 0; i < input_faces.size(); i++) {
        if(!is_face_inserted[i])///use currently inserted faces as mesh
            continue;
        const auto &f = input_faces[i];
        for (int j = 0; j < 3; j++) {
            //edges
            std::array<int, 2> e = {{f[j], f[(j + 1) % 3]}};
            if (e[0] > e[1])
                std::swap(e[0], e[1]);
            edges.push_back(e);
            //conn_tris
            conn_tris[input_faces[i][j]].push_back(i);
        }
    }
    vector_unique(edges);

    int cnt1 = 0;
    int cnt2 = 0;
    for (const auto &e: edges) {
        std::vector<int> n12_f_ids;
        std::set_intersection(conn_tris[e[0]].begin(), conn_tris[e[0]].end(),
                              conn_tris[e[1]].begin(), conn_tris[e[1]].end(), std::back_inserter(n12_f_ids));

        if (n12_f_ids.size() == 1) {//open boundary
            b_edge_infos.push_back(std::make_pair(e, n12_f_ids));
            cnt1++;
        } else {
            int f_id = n12_f_ids[0];
            int j = 0;
            for (; j < 3; j++) {
                if ((input_faces[f_id][j] == e[0] && input_faces[f_id][mod3(j + 1)] == e[1])
                    || (input_faces[f_id][j] == e[1] && input_faces[f_id][mod3(j + 1)] == e[0]))
                    break;
            }

            Vector3 n = get_normal(input_vertices[input_faces[f_id][0]], input_vertices[input_faces[f_id][1]],
                                   input_vertices[input_faces[f_id][2]]);
            int t = get_t(input_vertices[input_faces[f_id][0]], input_vertices[input_faces[f_id][1]],
                          input_vertices[input_faces[f_id][2]]);

            bool is_fine = false;
            for (int k = 0; k < n12_f_ids.size(); k++) {
                if (n12_f_ids[k] == f_id)
                    continue;
                Vector3 n1 = get_normal(input_vertices[input_faces[n12_f_ids[k]][0]],
                                        input_vertices[input_faces[n12_f_ids[k]][1]],
                                        input_vertices[input_faces[n12_f_ids[k]][2]]);
                if (abs(n1.dot(n)) < 1 - SCALAR_ZERO) {
                    is_fine = true;
                    break;
                }
            }
            if (is_fine)
                continue;

            is_fine = false;
            int ori;
            for (int k = 0; k < n12_f_ids.size(); k++) {
                for (int r = 0; r < 3; r++) {
                    if (input_faces[n12_f_ids[k]][r] != input_faces[f_id][j] &&
                        input_faces[n12_f_ids[k]][r] != input_faces[f_id][mod3(j + 1)]) {
                        if (k == 0) {
                            ori = Predicates::orient_2d(to_2d(input_vertices[input_faces[f_id][j]], t),
                                                        to_2d(input_vertices[input_faces[f_id][mod3(j + 1)]], t),
                                                        to_2d(input_vertices[input_faces[n12_f_ids[k]][r]], t));
                            break;
                        }
                        int new_ori = Predicates::orient_2d(to_2d(input_vertices[input_faces[f_id][j]], t),
                                                            to_2d(input_vertices[input_faces[f_id][mod3(j + 1)]],
                                                                  t),
                                                            to_2d(input_vertices[input_faces[n12_f_ids[k]][r]], t));
                        if (new_ori != ori)
                            is_fine = true;
                        break;
                    }
                }
                if (is_fine)
                    break;
            }
            if (is_fine)
                continue;

            cnt2++;
            b_edge_infos.push_back(std::make_pair(e, n12_f_ids));
        }
    }

    cout << "#boundary_e1 = " << cnt1 << endl;
    cout << "#boundary_e2 = " << cnt2 << endl;
}

bool floatTetWild::insert_boundary_edges(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
                                         std::vector<std::pair<std::array<int, 2>, std::vector<int>>>& b_edge_infos,
                                         std::vector<std::array<std::vector<int>, 4>> &track_surface_fs, Mesh& mesh,
                                         AABBWrapper &tree, std::vector<std::array<int, 2>>& b_edges,
                                         std::vector<bool> &is_face_inserted, bool is_again,
                                         std::vector<std::array<int, 3>>& known_surface_fs,
                                         std::vector<std::array<int, 3>>& known_not_surface_fs) {

    auto mark_known_surface_fs = [&](const std::array<int, 3>& f, int tag){
        std::vector<int> n_t_ids;
        set_intersection(mesh.tet_vertices[f[0]].conn_tets, mesh.tet_vertices[f[1]].conn_tets,
                         mesh.tet_vertices[f[2]].conn_tets, n_t_ids);
        if(n_t_ids.size()!=2)//todo:?????
            return;

        for(int t_id:n_t_ids){
            int j = get_local_f_id(t_id, f[0], f[1], f[2], mesh);
            mesh.tets[t_id].is_surface_fs[j] = tag;
        }
    };

    auto record_boundary_info = [&](const std::vector<Vector3>& points, const std::vector<int>& snapped_v_ids,
            const std::array<int, 2>& e){
        const int tet_vertices_size = mesh.tet_vertices.size();
        for (int i = points.size(); i > 0; i--)
            mesh.tet_vertices[tet_vertices_size - i].is_on_boundary = true;

        for (int v_id: snapped_v_ids) {
            Scalar dis_2 = p_seg_squared_dist_3d(mesh.tet_vertices[v_id].pos, input_vertices[e[0]],
                                                 input_vertices[e[1]]);
            if (dis_2 <= mesh.params.eps_2)
                mesh.tet_vertices[v_id].is_on_boundary = true;
        }

        b_edges.push_back(e);
    };

    bool is_all_inserted = true;
    int cnt = 0;
    for (int I = 0; I < b_edge_infos.size(); I++) {
        const auto &e = b_edge_infos[I].first;
        auto &n_f_ids = b_edge_infos[I].second;///it is sorted

        ///double check neighbors
        for (int i = 0; i < n_f_ids.size(); i++) {
            if (!is_face_inserted[n_f_ids[i]]) {
                n_f_ids.erase(n_f_ids.begin() + i);
                i--;
                break;
            }
        }
        if (n_f_ids.empty()) {
            cout << "FAIL n_f_ids.empty()" << endl;
            continue;
        }

        ///compute intersection
        std::vector<Vector3> points;
        std::map<std::array<int, 2>, int> map_edge_to_intersecting_point;
        std::vector<int> snapped_v_ids;
        std::vector<std::array<int, 3>> cut_fs;
        if (!insert_boundary_edges_get_intersecting_edges_and_points(input_vertices, input_faces, e, n_f_ids,
                                                                     track_surface_fs,
                                                                     mesh, points, map_edge_to_intersecting_point,
                                                                     snapped_v_ids, cut_fs,
                                                                     is_again)) {
            for (int f_id: n_f_ids)
                is_face_inserted[f_id] = false;
            is_all_inserted = false;

            cout << "FAIL insert_boundary_edges_get_intersecting_edges_and_points" << endl;
            continue;
        }
        if (points.empty()) { ///if all snapped
            record_boundary_info(points, snapped_v_ids, e);
            cnt++;
            continue;
        }

        ///subdivision
        std::vector<int> cut_t_ids;
        for (const auto &m: map_edge_to_intersecting_point) {
            const auto &tet_e = m.first;
            std::vector<int> tmp;
            set_intersection(mesh.tet_vertices[tet_e[0]].conn_tets, mesh.tet_vertices[tet_e[1]].conn_tets, tmp);
            cut_t_ids.insert(cut_t_ids.end(), tmp.begin(), tmp.end());
        }
        vector_unique(cut_t_ids);
        std::vector<bool> is_mark_surface(cut_t_ids.size(), false);
        CutMesh empty_cut_mesh(mesh, Vector3(0, 0, 0), std::array<Vector3, 3>());
        //
        std::vector<MeshTet> new_tets;
        std::vector<std::array<std::vector<int>, 4>> new_track_surface_fs;
        std::vector<int> modified_t_ids;
        std::vector<std::pair<int, int>> _;
        if (!subdivide_tets(-1, mesh, empty_cut_mesh, points, map_edge_to_intersecting_point, track_surface_fs,
                            cut_t_ids, is_mark_surface,
                            new_tets, new_track_surface_fs, modified_t_ids, _)) {
            bool is_inside_envelope = true;
            for (auto &f: cut_fs) {
                std::vector<GEO::vec3> ps;
                sample_triangle({{mesh.tet_vertices[f[0]].pos, mesh.tet_vertices[f[1]].pos,
                                         mesh.tet_vertices[f[2]].pos}}, ps, mesh.params.dd);
                if (tree.is_out_sf_envelope(ps, mesh.params.eps_2)) {
                    is_inside_envelope = false;
                    break;
                }
            }
            if (!is_inside_envelope) {
                for (int f_id: n_f_ids)
                    is_face_inserted[f_id] = false;
                for (auto &f: cut_fs) {
                    mark_known_surface_fs(f, KNOWN_NOT_SURFACE);
                }
                known_not_surface_fs.insert(known_not_surface_fs.end(), cut_fs.begin(), cut_fs.end());
                cout << "FAIL subdivide_tets" << endl;
            } else {
                for (auto &f: cut_fs)
                    mark_known_surface_fs(f, KNOWN_SURFACE);
                known_surface_fs.insert(known_surface_fs.end(), cut_fs.begin(), cut_fs.end());
            }

            is_all_inserted = false;//unless now
            continue;
        }

        //
        push_new_tets(mesh, track_surface_fs, points, new_tets, new_track_surface_fs, modified_t_ids, is_again);

        //
        ///mark boundary vertices
        ///notice, here we assume points list is inserted in the end of mesh.tet_vertices
        record_boundary_info(points, snapped_v_ids, e);
        cnt++;
    }

    logger().info("uninsert boundary #e = {}/{}", b_edge_infos.size() - cnt, b_edge_infos.size());

    return is_all_inserted;
}

bool floatTetWild::insert_boundary_edges_get_intersecting_edges_and_points(
        const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
        const std::array<int, 2> &e, const std::vector<int> &n_f_ids,
        std::vector<std::array<std::vector<int>, 4>> &track_surface_fs, Mesh &mesh,
        std::vector<Vector3>& points, std::map<std::array<int, 2>, int>& map_edge_to_intersecting_point,
        std::vector<int>& snapped_v_ids, std::vector<std::array<int, 3>>& cut_fs,
        bool is_again){

    auto is_cross = [](int a, int b) {
        if ((a == Predicates::ORI_POSITIVE && b == Predicates::ORI_NEGATIVE)
            || (a == Predicates::ORI_NEGATIVE && b == Predicates::ORI_POSITIVE))
            return true;
        return false;
    };

    int t = get_t(input_vertices[input_faces[n_f_ids.front()][0]],
                  input_vertices[input_faces[n_f_ids.front()][1]],
                  input_vertices[input_faces[n_f_ids.front()][2]]);
    std::array<Vector2, 2> evs_2d = {{to_2d(input_vertices[e[0]], t), to_2d(input_vertices[e[1]], t)}};

    std::vector<bool> is_visited(mesh.tets.size(), false);
    std::queue<int> t_ids_queue;
    ///find seed t_ids
    if (!is_again) {
        std::vector<int> t_ids;
        for (int v_id: e)
            t_ids.insert(t_ids.end(), mesh.tet_vertices[v_id].conn_tets.begin(),
                         mesh.tet_vertices[v_id].conn_tets.end());
        vector_unique(t_ids);
        for (int t_id: t_ids) {
            t_ids_queue.push(t_id);
            is_visited[t_id] = true;
        }
    } else {
        Vector3 min_e, max_e;
        get_bbox_face(input_vertices[e[0]], input_vertices[e[0]], input_vertices[e[1]], min_e, max_e);

        for (int t_id = 0; t_id < mesh.tets.size(); t_id++) {
            if (mesh.tets[t_id].is_removed)
                continue;

            Vector3 min_t, max_t;
            get_bbox_tet(mesh.tet_vertices[mesh.tets[t_id][0]].pos, mesh.tet_vertices[mesh.tets[t_id][1]].pos,
                         mesh.tet_vertices[mesh.tets[t_id][2]].pos, mesh.tet_vertices[mesh.tets[t_id][3]].pos,
                         min_t, max_t);
            if (!is_bbox_intersected(min_e, max_e, min_t, max_t))
                continue;

            t_ids_queue.push(t_id);
            is_visited[t_id] = true;
        }
    }

    std::vector<std::array<int, 3>> f_oris;
    while (!t_ids_queue.empty()) {
        int t_id = t_ids_queue.front();
        t_ids_queue.pop();

        std::array<bool, 4> is_cut_vs = {{false, false, false, false}};
        for (int j = 0; j < 4; j++) {
            ///check if contains
            if(track_surface_fs[t_id][j].empty())
                continue;
            std::sort(track_surface_fs[t_id][j].begin(), track_surface_fs[t_id][j].end());
            std::vector<int> tmp;
            std::set_intersection(track_surface_fs[t_id][j].begin(), track_surface_fs[t_id][j].end(),
                                  n_f_ids.begin(), n_f_ids.end(), std::back_inserter(tmp));
            if (tmp.empty())
                continue;

            ///check if cut through
            //check tri side of seg
            std::array<int, 3> f_v_ids = {{mesh.tets[t_id][(j + 1) % 4], mesh.tets[t_id][(j + 2) % 4],
                                                  mesh.tets[t_id][(j + 3) % 4]}};
            std::array<Vector2, 3> fvs_2d = {{to_2d(mesh.tet_vertices[f_v_ids[0]].pos, t),
                                                     to_2d(mesh.tet_vertices[f_v_ids[1]].pos, t),
                                                     to_2d(mesh.tet_vertices[f_v_ids[2]].pos, t)}};
            int cnt_pos = 0;
            int cnt_neg = 0;
            int cnt_on = 0;
            std::array<int, 3> oris;
            for (int k = 0; k < 3; k++) {
                oris[k] = Predicates::orient_2d(evs_2d[0], evs_2d[1], fvs_2d[k]);
                if (oris[k] == Predicates::ORI_ZERO) {
                    cnt_on++;
                } else {
                    Scalar dis_2 = p_line_squared_dist_3d(mesh.tet_vertices[f_v_ids[k]].pos, input_vertices[e[0]],
                                                         input_vertices[e[1]]);
                    if (dis_2 < mesh.params.eps_2_coplanar) {
                        oris[k] = Predicates::ORI_ZERO;
                        cnt_on++;
                        continue;
                    }
                    if (oris[k] == Predicates::ORI_POSITIVE)
                        cnt_pos++;
                    else
                        cnt_neg++;
                }
            }
            if (cnt_on >= 2) {
                cut_fs.push_back(f_v_ids);
                std::sort(cut_fs.back().begin(), cut_fs.back().end());
                f_oris.push_back(oris);
                continue;
            }
            if (cnt_neg == 0 || cnt_pos == 0)
                continue;

            //check tri edge - seg intersection
            bool is_intersected = false;
            for (auto &p: evs_2d) { ///first check if endpoints are contained inside the triangle
                if (is_p_inside_tri_2d(p, fvs_2d)) {
                    is_intersected = true;
                    break;
                }
            }
            if(!is_intersected) { ///then check if there's intersection
                for (int k = 0; k < 3; k++) {
                    //if cross
                    if (!is_cross(oris[k], oris[(k + 1) % 3]))
                        continue;
                    //if already know intersect
                    std::array<int, 2> tri_e = {{f_v_ids[k], f_v_ids[(k + 1) % 3]}};
                    if (tri_e[0] > tri_e[1])
                        std::swap(tri_e[0], tri_e[1]);
                    if (map_edge_to_intersecting_point.find(tri_e) != map_edge_to_intersecting_point.end()) {
                        is_intersected = true;
                        break;
                    }
                    //if intersect
                    double t2 = -1;
                    if (seg_seg_intersection_2d(evs_2d, {{fvs_2d[k], fvs_2d[(k + 1) % 3]}}, t2)) {
                        Vector3 p = (1 - t2) * mesh.tet_vertices[f_v_ids[k]].pos
                                    + t2 * mesh.tet_vertices[f_v_ids[(k + 1) % 3]].pos;
                        double dis1 = (p - mesh.tet_vertices[f_v_ids[k]].pos).squaredNorm();
                        double dis2 = (p - mesh.tet_vertices[f_v_ids[(k + 1) % 3]].pos).squaredNorm();
                        if (dis1 < SCALAR_ZERO_2) {
                            oris[k] = Predicates::ORI_ZERO;
                            is_intersected = true;
                            break;
                        }
                        if (dis2 < SCALAR_ZERO_2) {
                            oris[k] = Predicates::ORI_ZERO;
                            is_intersected = true;
                            break;
                        }
                        points.push_back(p);
                        map_edge_to_intersecting_point[tri_e] = points.size() - 1;
                        is_intersected = true;
                        break;
                    }
                }
                if (!is_intersected) /// no need to return false here
                    continue;
            }

            cut_fs.push_back(f_v_ids);
            std::sort(cut_fs.back().begin(), cut_fs.back().end());
            f_oris.push_back(oris);
            is_cut_vs[(j + 1) % 4] = true;
            is_cut_vs[(j + 2) % 4] = true;
            is_cut_vs[(j + 3) % 4] = true;
        }
        for (int j = 0; j < 4; j++) {
            if (!is_cut_vs[j])
                continue;
            for (int n_t_id: mesh.tet_vertices[mesh.tets[t_id][j]].conn_tets) {
                if (!is_visited[n_t_id]) {
                    t_ids_queue.push(n_t_id);
                    is_visited[n_t_id] = true;
                }
            }
        }
    }
//    cout<<"cut_fs.size = "<<cut_fs.size()<<endl;
    {//remove duplicated elements
        std::vector<int> indices(cut_fs.size());
        for (int i = 0; i < cut_fs.size(); i++)
            indices[i] = i;
        std::sort(indices.begin(), indices.end(), [&cut_fs](int i1, int i2) {
            return cut_fs[i1] < cut_fs[i2];
        });
        indices.erase(std::unique(indices.begin(), indices.end(), [&cut_fs](int i1, int i2) {
            return cut_fs[i1] < cut_fs[i2];
        }), indices.end());
        std::vector<std::array<int, 3>> new_cut_fs(indices.size());
        std::vector<std::array<int, 3>> new_f_oris(indices.size());
        for (int i = 0; i < indices.size(); i++) {
            new_cut_fs[i] = cut_fs[indices[i]];
            new_f_oris[i] = f_oris[indices[i]];
        }
        cut_fs = new_cut_fs;
        f_oris = new_f_oris;
    }

    for (int i = 0; i < cut_fs.size(); i++) {
        std::array<bool, 3> is_e_intersected = {{false, false, false}};
        int cnt = 0;
        for (int j = 0; j < 3; j++) {
            if (f_oris[i][j] == Predicates::ORI_ZERO) {
                cnt++;
                is_e_intersected[j] = true;
                snapped_v_ids.push_back(cut_fs[i][j]);
                continue;
            }
            std::array<int, 2> tri_e = {{cut_fs[i][j], cut_fs[i][(j + 1) % 3]}};
            if (tri_e[0] > tri_e[1])
                std::swap(tri_e[0], tri_e[1]);
            if (map_edge_to_intersecting_point.find(tri_e) != map_edge_to_intersecting_point.end()) {
                cnt++;
                is_e_intersected[j] = true;
                continue;
            }

            if (f_oris[i][(j + 1) % 3] == Predicates::ORI_ZERO)
                is_e_intersected[j] = true;
        }
//        if (cnt == 2)
        if (cnt >= 2)
            continue;

        //line - tri edges intersection
        for (int j = 0; j < 3; j++) {
            if (is_e_intersected[j])
                continue;
            if (!is_cross(f_oris[i][j], f_oris[i][(j + 1) % 3]))
                continue;

            std::array<Vector2, 2> tri_evs_2d = {{to_2d(mesh.tet_vertices[cut_fs[i][j]].pos, t),
                                                         to_2d(mesh.tet_vertices[cut_fs[i][(j + 1) % 3]].pos, t)}};
            Scalar t_seg = -1;
            if (seg_line_intersection_2d(tri_evs_2d, evs_2d, t_seg)) {
                std::array<int, 2> tri_e = {{cut_fs[i][j], cut_fs[i][(j + 1) % 3]}};
                if (tri_e[0] > tri_e[1])
                    std::swap(tri_e[0], tri_e[1]);
                points.push_back((1 - t_seg) * mesh.tet_vertices[cut_fs[i][j]].pos
                                 + t_seg * mesh.tet_vertices[cut_fs[i][(j + 1) % 3]].pos);
                map_edge_to_intersecting_point[tri_e] = points.size() - 1;
                cnt++;
            }
        }
//        if (cnt != 2) {
        if (cnt < 2) {
//            cout<<"return 2, cnt = "<<cnt<<endl;
            return false;///has to return false here
        }
    }

    return true;
}

void floatTetWild::mark_surface_fs(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
                                   std::vector<std::array<std::vector<int>, 4>> &track_surface_fs,
                                   const std::vector<bool> &is_face_inserted,
                                   const std::vector<std::array<int, 3>>& known_surface_fs,
                                   const std::vector<std::array<int, 3>>& known_not_surface_fs,
                                   std::vector<std::array<int, 2>>& b_edges,
                                   Mesh &mesh, AABBWrapper &tree) {

    auto is_on_bounded_side = [&](const std::array<Vector2, 3> &ps_2d, const Vector2 &c) {
        int cnt_pos = 0;
        int cnt_neg = 0;
        for (int j = 0; j < 3; j++) {
            int ori = Predicates::orient_2d(ps_2d[j], ps_2d[(j + 1) % 3], c);
            if (ori == Predicates::ORI_POSITIVE)
                cnt_pos++;
            else if (ori == Predicates::ORI_NEGATIVE)
                cnt_neg++;
        }

        if (cnt_neg > 0 && cnt_pos > 0)
            return false;
        return true;
    };

    std::vector<std::array<bool, 4>> is_visited(track_surface_fs.size(), {{false, false, false, false}});
    for (int t_id = 0; t_id < track_surface_fs.size(); t_id++) {
        for (int j = 0; j < 4; j++) {
//            //fortest
//            if (track_surface_fs[t_id][j].size() > 0)
////            if (std::find(track_surface_fs[t_id][j].begin(), track_surface_fs[t_id][j].end(), III) != track_surface_fs[t_id][j].end())
//                mesh.tets[t_id].is_surface_fs[j] = -1;
//            continue;
//            //fortest

            //
//            if (mesh.tets[t_id].is_surface_fs[j] != NOT_SURFACE || is_visited[t_id][j])
//                continue;
//            is_visited[t_id][j] = true;
//            if (track_surface_fs[t_id][j].empty())
//                continue;
//            //
//            int opp_t_id = get_opp_t_id(t_id, j, mesh);
//            if (opp_t_id < 0)
//                continue;
//            int k = get_local_f_id(opp_t_id, mesh.tets[t_id][(j + 1) % 4], mesh.tets[t_id][(j + 2) % 4],
//                                   mesh.tets[t_id][(j + 3) % 4], mesh);
//            is_visited[opp_t_id][k] = true;
//            //
//            std::sort(track_surface_fs[t_id][j].begin(), track_surface_fs[t_id][j].end());
//            std::sort(track_surface_fs[opp_t_id][k].begin(), track_surface_fs[opp_t_id][k].end());
//            std::vector<int> f_ids;
//            if (track_surface_fs[t_id][j] != track_surface_fs[opp_t_id][k])
//                std::set_union(track_surface_fs[t_id][j].begin(), track_surface_fs[t_id][j].end(),
//                               track_surface_fs[opp_t_id][k].begin(), track_surface_fs[opp_t_id][k].end(),
//                               std::back_inserter(f_ids));
//            else
//                f_ids = track_surface_fs[t_id][j];

//            //fortest
//            if (mesh.tets[t_id].is_surface_fs[j] != NOT_SURFACE) {
//                auto &tp1_3d = mesh.tet_vertices[mesh.tets[t_id][(j + 1) % 4]].pos;
//                auto &tp2_3d = mesh.tet_vertices[mesh.tets[t_id][(j + 2) % 4]].pos;
//                auto &tp3_3d = mesh.tet_vertices[mesh.tets[t_id][(j + 3) % 4]].pos;
//                std::vector<GEO::vec3> ps;
//                sample_triangle({{tp1_3d, tp2_3d, tp3_3d}}, ps, mesh.params.dd);
//                if (tree.is_out_sf_envelope(ps, mesh.params.eps_2)) {
//                    cout << "tree.is_out_sf_envelope(ps, mesh.params.eps_2)" << endl;
//                    cout << "mesh.tets[t_id].is_surface_fs[j] = " << (int) mesh.tets[t_id].is_surface_fs[j] << endl;
//                    cout << "t_id = " << t_id << endl;
//                    cout << "j = " << j << endl;
//                    cout<<"is_visited[t_id][j] = "<<is_visited[t_id][j]<<endl;
//                    pausee();
//                }
//            }
//            //fortest

            int ff_id = -1;
            int opp_t_id = -1;
            int k = -1;
            int case_id = 0;
            if (mesh.tets[t_id].is_surface_fs[j] == KNOWN_SURFACE
                || mesh.tets[t_id].is_surface_fs[j] == KNOWN_NOT_SURFACE) {
                opp_t_id = get_opp_t_id(t_id, j, mesh);
                if (opp_t_id < 0) {
                    mesh.tets[t_id].is_surface_fs[j] = NOT_SURFACE;
                    continue;
                }
                k = get_local_f_id(opp_t_id, mesh.tets[t_id][(j + 1) % 4], mesh.tets[t_id][(j + 2) % 4],
                                   mesh.tets[t_id][(j + 3) % 4], mesh);
                is_visited[t_id][j] = true;
                is_visited[opp_t_id][k] = true;
                if (mesh.tets[t_id].is_surface_fs[j] == KNOWN_NOT_SURFACE || track_surface_fs[t_id][j].empty()) {
                    mesh.tets[t_id].is_surface_fs[j] = NOT_SURFACE;
                    mesh.tets[opp_t_id].is_surface_fs[k] = NOT_SURFACE;
                    continue;
                } else
                    ff_id = track_surface_fs[t_id][j].front();
                case_id = 1;
            } else {
                if (mesh.tets[t_id].is_surface_fs[j] != NOT_SURFACE || is_visited[t_id][j])
                    continue;
                is_visited[t_id][j] = true;
                if (track_surface_fs[t_id][j].empty())
                    continue;

                auto &f_ids = track_surface_fs[t_id][j];

                auto &tp1_3d = mesh.tet_vertices[mesh.tets[t_id][(j + 1) % 4]].pos;
                auto &tp2_3d = mesh.tet_vertices[mesh.tets[t_id][(j + 2) % 4]].pos;
                auto &tp3_3d = mesh.tet_vertices[mesh.tets[t_id][(j + 3) % 4]].pos;

//                std::vector<GEO::vec3> ps;
//                sample_triangle({{tp1_3d, tp2_3d, tp3_3d}}, ps, mesh.params.dd);
//                double eps = (mesh.params.eps + mesh.params.eps_simplification) / 2;
//                eps *= eps;//fortest: use a smaller eps here
//                if (tree.is_out_sf_envelope(ps, eps))
//                    continue;
//                else
//                    ff_id = track_surface_fs[t_id][j].front();

                int t = get_t(tp1_3d, tp2_3d, tp3_3d);
                std::array<Vector2, 3> tps_2d = {{to_2d(tp1_3d, t), to_2d(tp2_3d, t), to_2d(tp3_3d, t)}};
                Vector2 c = (tps_2d[0] + tps_2d[1] + tps_2d[2]) / 3;

                std::array<Vector2, 3> ps_2d;
                for (int f_id: f_ids) {
                    if (!is_face_inserted[f_id])
                        continue;

                    ps_2d = {{to_2d(input_vertices[input_faces[f_id][0]], t),
                                     to_2d(input_vertices[input_faces[f_id][1]], t),
                                     to_2d(input_vertices[input_faces[f_id][2]], t)}};
                    if (!is_on_bounded_side(ps_2d, c))
                        continue;

                    ff_id = f_id;
                    break;
                }
                if (ff_id < 0) {
                    continue;
//                    std::vector<GEO::vec3> ps;
//                    sample_triangle({{tp1_3d, tp2_3d, tp3_3d}}, ps, mesh.params.dd);
//                    if (tree.is_out_sf_envelope(ps, mesh.params.eps_2))
//                        continue;
//                    else
//                        ff_id = track_surface_fs[t_id][j].front();
                }

                opp_t_id = get_opp_t_id(t_id, j, mesh);
                if (opp_t_id < 0) {
                    mesh.tets[t_id].is_surface_fs[j] = NOT_SURFACE;
                    continue;
                }
                k = get_local_f_id(opp_t_id, mesh.tets[t_id][(j + 1) % 4], mesh.tets[t_id][(j + 2) % 4],
                                   mesh.tets[t_id][(j + 3) % 4], mesh);
                is_visited[opp_t_id][k] = true;

                case_id = 2;
            }

//            //fortest
//            auto &tp1_3d = mesh.tet_vertices[mesh.tets[opp_t_id][(k + 1) % 4]].pos;
//            auto &tp2_3d = mesh.tet_vertices[mesh.tets[opp_t_id][(k + 2) % 4]].pos;
//            auto &tp3_3d = mesh.tet_vertices[mesh.tets[opp_t_id][(k + 3) % 4]].pos;
//            std::vector<GEO::vec3> ps;
//            sample_triangle({{tp1_3d, tp2_3d, tp3_3d}}, ps, mesh.params.dd);
//            if (tree.is_out_sf_envelope(ps, mesh.params.eps_2)) {
//                cout << "OPP tree.is_out_sf_envelope(ps, mesh.params.eps_2)" << endl;
//                cout << "mesh.tets[t_id].is_surface_fs[j] = " << (int) mesh.tets[t_id].is_surface_fs[j] << endl;
//                cout << "t_id = " << t_id << endl;
//                cout << "j = " << j << endl;
//                cout << mesh.tets[opp_t_id][(k + 1) % 4] << " "
//                     << mesh.tets[opp_t_id][(k + 2) % 4] << " "
//                     << mesh.tets[opp_t_id][(k + 3) % 4] << endl;
//                cout << mesh.tets[t_id][(j + 1) % 4] << " "
//                     << mesh.tets[t_id][(j + 2) % 4] << " "
//                     << mesh.tets[t_id][(j + 3) % 4] << endl;
//                pausee();
//            }
//            //fortest

            myassert(ff_id>=0, "ff_id<0!!!");//fortest

            auto &fv1 = input_vertices[input_faces[ff_id][0]];
            auto &fv2 = input_vertices[input_faces[ff_id][1]];
            auto &fv3 = input_vertices[input_faces[ff_id][2]];
            //
            int ori = Predicates::orient_3d(fv1, fv2, fv3, mesh.tet_vertices[mesh.tets[t_id][j]].pos);
            int opp_ori = Predicates::orient_3d(fv1, fv2, fv3, mesh.tet_vertices[mesh.tets[opp_t_id][k]].pos);
            //
            if (ori == Predicates::ORI_POSITIVE && opp_ori == Predicates::ORI_NEGATIVE
                || ori == Predicates::ORI_NEGATIVE && opp_ori == Predicates::ORI_POSITIVE) {
                mesh.tets[t_id].is_surface_fs[j] = ori;
                mesh.tets[opp_t_id].is_surface_fs[k] = opp_ori;
                continue;
            }
            //
            if (ori == Predicates::ORI_ZERO && opp_ori != Predicates::ORI_ZERO) {
                mesh.tets[t_id].is_surface_fs[j] = -opp_ori;
                mesh.tets[opp_t_id].is_surface_fs[k] = opp_ori;
                continue;
            }
            if (opp_ori == Predicates::ORI_ZERO && ori != Predicates::ORI_ZERO) {
                mesh.tets[t_id].is_surface_fs[j] = ori;
                mesh.tets[opp_t_id].is_surface_fs[k] = -ori;
                continue;
            }
            //
            if (ori == opp_ori) {
                Vector3 n = (fv2 - fv1).cross(fv3 - fv1);
                n.normalize();
                Scalar dist = n.dot(mesh.tet_vertices[mesh.tets[t_id][j]].pos - fv1);
                Scalar opp_dist = n.dot(mesh.tet_vertices[mesh.tets[opp_t_id][k]].pos - fv1);
                if (ori == Predicates::ORI_ZERO) {
//                    cout << "impossible!! " << dist << " " << opp_dist << endl;
//                    cout << (mesh.tet_vertices[mesh.tets[t_id][j]].pos == mesh.tet_vertices[mesh.tets[opp_t_id][k]].pos)
//                         << endl;
//                    cout << t_id << " " << j << " " << mesh.tets[t_id][j] << endl;
//                    cout << opp_t_id << " " << k << " " << mesh.tets[opp_t_id][k] << endl;
//                    cout << "n = " << n << endl;

                    auto &tv1 = mesh.tet_vertices[mesh.tets[t_id][(j + 1) % 4]].pos;
                    auto &tv2 = mesh.tet_vertices[mesh.tets[t_id][(j + 2) % 4]].pos;
                    auto &tv3 = mesh.tet_vertices[mesh.tets[t_id][(j + 3) % 4]].pos;
                    Vector3 nt;
                    if (Predicates::orient_3d(tv1, tv2, tv3, mesh.tet_vertices[mesh.tets[t_id][j]].pos)
                        == Predicates::ORI_POSITIVE)
                        nt = (tv2 - tv1).cross(tv3 - tv1);
                    else
                        nt = (tv3 - tv1).cross(tv2 - tv1);
                    nt.normalize();
                    if (n.dot(nt) > 0) {
                        mesh.tets[opp_t_id].is_surface_fs[k] = 1;
                        mesh.tets[t_id].is_surface_fs[j] = -1;
                    } else {
                        mesh.tets[opp_t_id].is_surface_fs[k] = -1;
                        mesh.tets[t_id].is_surface_fs[j] = 1;
                    }
                } else {
                    if (dist < opp_dist) {
                        mesh.tets[opp_t_id].is_surface_fs[k] = opp_ori;
                        mesh.tets[t_id].is_surface_fs[j] = -ori;
                    } else {
                        mesh.tets[opp_t_id].is_surface_fs[k] = -opp_ori;
                        mesh.tets[t_id].is_surface_fs[j] = ori;
                    }
                }
            }
            //
        }
    }

//    fortest: output and check
    output_surface(mesh, mesh.params.output_path+"surface.stl");

    cout<<"known_surface_fs.size = "<<known_surface_fs.size()<<endl;
    cout<<"known_not_surface_fs.size = "<<known_not_surface_fs.size()<<endl;
    if(known_surface_fs.empty() && known_not_surface_fs.empty())
        return;

    //b_edges
    std::vector<std::array<int, 2>> tmp_edges;
    for (const auto &f: known_surface_fs) {
        for (int j = 0; j < 3; j++) {
            if (f[j] < f[(j + 1) % 3])
                tmp_edges.push_back({{f[j], f[(j + 1) % 3]}});
            else
                tmp_edges.push_back({{f[(j + 1) % 3], f[j]}});
        }
    }
    for (const auto &f: known_not_surface_fs) {
        for (int j = 0; j < 3; j++) {
            if (f[j] < f[(j + 1) % 3])
                tmp_edges.push_back({{f[j], f[(j + 1) % 3]}});
            else
                tmp_edges.push_back({{f[(j + 1) % 3], f[j]}});
        }
    }
    vector_unique(tmp_edges);

    for (auto &e: tmp_edges) {
        int cnt = 0;
        std::vector<int> n_t_ids;
        set_intersection(mesh.tet_vertices[e[0]].conn_tets, mesh.tet_vertices[e[1]].conn_tets, n_t_ids);
        for (int t_id: n_t_ids) {
            for (int j = 0; j < 4; j++) {
                if (mesh.tets[t_id][j] != e[0] && mesh.tets[t_id][j] != e[1]) {
                    if (mesh.tets[t_id].is_surface_fs[j] != NOT_SURFACE) {
                        cnt++;
                        if (cnt > 2)
                            break;
                    }
                }
                if (cnt > 2)
                    break;
            }
        }
//        cout<<"cnt = "<<cnt<<endl;
        if (cnt == 2) {
            b_edges.push_back(e);
            mesh.tet_vertices[e[0]].is_on_boundary = true;
//            cout<<"b_edges.push_back(e);"<<endl;
        }
    }
}

void floatTetWild::mark_boundary_vs(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
                      std::vector<std::array<std::vector<int>, 4>> &track_surface_fs,
                      const std::vector<std::pair<std::array<int, 2>, std::vector<int>>>& b_edge_infos,
                      Mesh &mesh, std::vector<std::array<int, 2>> &b_edges) {
    std::vector<std::vector<std::pair<int, int>>> covered_fs_infos(input_faces.size());
    for (int i = 0; i < track_surface_fs.size(); i++) {
        for (int j = 0; j < 4; j++) {
            for (int f_id: track_surface_fs[i][j])
                covered_fs_infos[f_id].push_back(std::make_pair(i, j));
        }
    }

    b_edges.resize(b_edge_infos.size());
    for (int i = 0; i < b_edge_infos.size(); i++) {
        const auto &e = b_edge_infos[i].first;
        auto &n_f_ids = b_edge_infos[i].second;///it is sorted

        b_edges[i] = e;

        std::vector<int> v_ids;
        for (int f_id: n_f_ids) {
            for (const auto &info: covered_fs_infos[f_id]) {
                int t_id = info.first;
                int j = info.second;
                for (int k = 0; k < 3; k++)
                    v_ids.push_back(mesh.tets[t_id][(j + 1 + k) % 4]);
            }
        }
        vector_unique(v_ids);

        for (int v_id: v_ids) {
            for (int j = 0; j < 3; j++) {
                double dist_2 = p_line_squared_dist_3d(mesh.tet_vertices[v_id].pos,
                                                       input_vertices[e[0]], input_vertices[e[1]]);
                if (dist_2 <= mesh.params.eps_2_coplanar * 1.1) {
                    mesh.tet_vertices[v_id].is_on_boundary = true;
                    mesh.tet_vertices[v_id].on_boundary_e_id = i;
                }
            }
        }
    }
}

int floatTetWild::get_opp_t_id(int t_id, int j, const Mesh &mesh){
    std::vector<int> tmp;
    set_intersection(mesh.tet_vertices[mesh.tets[t_id][(j + 1) % 4]].conn_tets,
                     mesh.tet_vertices[mesh.tets[t_id][(j + 2) % 4]].conn_tets,
                     mesh.tet_vertices[mesh.tets[t_id][(j + 3) % 4]].conn_tets,
                     tmp);
//    //fortest
//    if(tmp.size() != 1 && tmp.size() != 2) {
//        cout << "tmp.size() = " << tmp.size() << endl;
//        cout << "t_id = " << t_id << ", j = " << j << endl;
//        vector_print(mesh.tet_vertices[mesh.tets[t_id][(j + 1) % 4]].conn_tets);
//        vector_print(mesh.tet_vertices[mesh.tets[t_id][(j + 2) % 4]].conn_tets);
//        vector_print(mesh.tet_vertices[mesh.tets[t_id][(j + 3) % 4]].conn_tets);
//        vector_print(tmp);
//        for(int i: tmp){
//            mesh.tets[i].print();
//        }
//        pausee();
//    }
//    //fortest
    if (tmp.size() == 2)
        return tmp[0] == t_id ? tmp[1] : tmp[0];
    else
        return -1;
}

void floatTetWild::maintain_covered_region(const std::vector<std::pair<int, int>> &old_covered_tet_fs,
                                           const std::vector<int> &cut_t_ids,
                                           std::vector<std::pair<int, int>> &new_covered_tet_fs) {
    std::unordered_set<int> set_cut_t_ids(cut_t_ids.begin(), cut_t_ids.end());
    for (const auto &f: old_covered_tet_fs) {
        if (set_cut_t_ids.find(f.first) != set_cut_t_ids.end())
            continue;
        new_covered_tet_fs.push_back(f);
    }
}

bool floatTetWild::preserve_edges(Mesh& mesh, std::vector<std::array<std::vector<int>, 4>> &track_surface_fs,
        int insert_f_id, const std::array<Vector3, 3>& f_vs, const Vector3& n, int t,
        int e_id,
        std::vector<std::pair<int, int>> &covered_tet_fs) {
    //
    std::vector<std::array<int, 2>> edges;
    for (const auto &f: covered_tet_fs) {
        int t_id = f.first;
        int j = f.second;
        std::array<int, 3> v_ids(
                {{mesh.tets[t_id][(j + 1) % 4], mesh.tets[t_id][(j + 2) % 4], mesh.tets[t_id][(j + 3) % 4]}});
        for (int k = 0; k < 3; k++) {
            if (v_ids[k] < v_ids[(k + 1) % 3])
                edges.push_back({{v_ids[k], v_ids[(k + 1) % 3]}});
            else
                edges.push_back({{v_ids[(k + 1) % 3], v_ids[k]}});
        }
    }
    vector_unique(edges);
    //
    std::vector<int> v_ids;
    for (const auto &e: edges) {
        v_ids.push_back(e[0]);
        v_ids.push_back(e[1]);
    }
    vector_unique(v_ids);
    //
    std::map<int, int> map_v_ids;
    std::vector<Vector2> vs_2d(v_ids.size());
    std::vector<int> oris(v_ids.size());
    Vector2 v1 = to_2d(f_vs[e_id], t);
    Vector2 v2 = to_2d(f_vs[(e_id + 1) % 3], t);
    for (int i = 0; i < v_ids.size(); i++) {
        map_v_ids[v_ids[i]] = i;
        vs_2d[i] = to_2d(mesh.tet_vertices[v_ids[i]].pos, t);//todo: project to plane first
        oris[i] = Predicates::orient_2d(v1, v2, vs_2d[i]);
    }
    //
    std::vector<int> subdivide_t_ids;
    std::vector<Vector3> points;
    std::map<std::array<int, 2>, int> map_edge_to_intersecting_point;
    for(const auto& e: edges) {
        if (!is_crossing(oris[map_v_ids[e[0]]], oris[map_v_ids[e[1]]]))
            continue;
        bool is_snapped = false;
        for (int j = 0; j < 2; j++) {
            Scalar dis_2 = p_line_squared_dist_3d(mesh.tet_vertices[e[0]].pos, f_vs[e_id], f_vs[(e_id + 1) % 3]);
            if (dis_2 < mesh.params.eps_2_coplanar) {
                is_snapped = true;
                break;
            }
        }
        if (is_snapped)
            continue;
        double t2;
        if (!seg_seg_intersection_2d({{v1, v2}}, {{vs_2d[map_v_ids[e[0]]], vs_2d[map_v_ids[e[1]]]}}, t2))
            continue;
        points.push_back((1 - t2) * mesh.tet_vertices[e[0]].pos + t2 * mesh.tet_vertices[e[1]].pos);
        map_edge_to_intersecting_point[e] = points.size() - 1;
        std::vector<int> n_t_ids;
        set_intersection(mesh.tet_vertices[e[0]].conn_tets, mesh.tet_vertices[e[1]].conn_tets, n_t_ids);
        subdivide_t_ids.insert(subdivide_t_ids.begin(), n_t_ids.begin(), n_t_ids.end());
    }
    vector_unique(subdivide_t_ids);
    //
    std::vector<MeshTet> new_tets;
    std::vector<std::array<std::vector<int>, 4>> new_track_surface_fs;
    std::vector<int> modified_t_ids;
    std::vector<std::pair<int, int>> new_covered_tet_fs;
    if(!subdivide_tets_edge(insert_f_id, mesh, points, map_edge_to_intersecting_point, track_surface_fs,
                   subdivide_t_ids, new_tets, new_track_surface_fs, modified_t_ids, new_covered_tet_fs))
        return false;
    //
    std::map<int, int> map_t_ids = push_new_tets(mesh, track_surface_fs, points, new_tets, new_track_surface_fs,
            modified_t_ids, false);
    for(auto& f: new_covered_tet_fs)
        f.first = map_t_ids[f.first];
    //
    maintain_covered_region(covered_tet_fs, subdivide_t_ids, new_covered_tet_fs);
    covered_tet_fs = new_covered_tet_fs;

    return true;
}

void floatTetWild::myassert(bool b, const std::string& s) {
    if (b == false) {
        cout << "myassert fail: " << s << endl;
        pausee();
    }
}

void floatTetWild::check_track_surface_fs(Mesh &mesh, std::vector<std::array<std::vector<int>, 4>> &track_surface_fs){
    return;

//    //check connection
//    std::vector<std::vector<int>> conn_tets(mesh.tet_vertices.size());
//    for (int i = 0; i < mesh.tets.size(); i++) {
//        for (int j = 0; j < 4; j++)
//            conn_tets[mesh.tets[i][j]].push_back(i);
//    }
//    for (int i = 0; i < mesh.tet_vertices.size(); i++) {
//        std::sort(mesh.tet_vertices[i].conn_tets.begin(), mesh.tet_vertices[i].conn_tets.end());
//        if (mesh.tet_vertices[i].conn_tets != conn_tets[i]) {
//            cout << "mesh.tet_vertices[i].conn_tets!=conn_tets[i]" << endl;
//            pausee();
//        }
//    }
//    cout<<"check 1 done"<<endl;
//
//    for (int i = 0; i < mesh.tets.size(); i++) {
//        for (int j = 0; j < 4; j++) {
//            int opp_t_id = get_opp_t_id(i, j, mesh);
//        }
//    }
//    cout<<"check 2 done"<<endl;

    for (int t_id = 0; t_id < track_surface_fs.size(); t_id++) {
        for (int j = 0; j < 4; j++) {
            std::sort(track_surface_fs[t_id][j].begin(), track_surface_fs[t_id][j].end());
            int opp_t_id = get_opp_t_id(t_id, j, mesh);
            if (opp_t_id < 0) {
                if (!track_surface_fs[t_id][j].empty()) {
                    cout << "bbox face but !track_surface_fs[t_id][j].empty()" << endl;
                    pausee();
                }
                continue;
            }
            int k = 0;
            for (k = 0; k < 4; k++) {
                if (mesh.tets[opp_t_id][k] != mesh.tets[t_id][(j + 1) % 4]
                    && mesh.tets[opp_t_id][k] != mesh.tets[t_id][(j + 2) % 4]
                    && mesh.tets[opp_t_id][k] != mesh.tets[t_id][(j + 3) % 4])
                    break;
            }
            std::sort(track_surface_fs[opp_t_id][k].begin(), track_surface_fs[opp_t_id][k].end());
            if (track_surface_fs[t_id][j] != track_surface_fs[opp_t_id][k]) {
                cout << "track_surface_fs[t_id][j]!=track_surface_fs[opp_t_id][k]" << endl;
                cout << "t " << t_id << ": ";
                vector_print(track_surface_fs[t_id][j]);
                cout << "opp_t " << opp_t_id << ": ";
                vector_print(track_surface_fs[opp_t_id][k]);
                pausee();
            }
        }
    }
    cout<<"check 3 done"<<endl;
}

#include <igl/writeOFF.h>
void floatTetWild::check_is_input_face_covered(const std::array<Vector3, 3> f_vs,
        const Mesh& mesh, const std::vector<std::pair<int, int>>& t_fs) {
    return;

    //
    std::vector<GEO::vec3> ps;
    sample_triangle(f_vs, ps, mesh.params.dd);
    //
    double eps = mesh.params.eps_coplanar * mesh.params.eps_coplanar;
    //
    if (ps.size() < 4) {
        cout << "ps.size() < 4!!" << endl;
        ps.clear();
        sample_triangle(f_vs, ps, mesh.params.dd / 4);
    }
    for (const auto &p: ps) {
        bool is_inside = false;
        for (const auto &info: t_fs) {
            int t_id = info.first;
            int j = info.second;
            std::array<int, 3> f = {{mesh.tets[t_id][(j + 1) % 4], mesh.tets[t_id][(j + 2) % 4],
                                     mesh.tets[t_id][(j + 3) % 4]}};
            std::array<GEO::vec3, 3> vs;
            for (int k = 0; k < 3; k++) {
                vs[k] = GEO::vec3(mesh.tet_vertices[f[k]].pos[0], mesh.tet_vertices[f[k]].pos[1],
                                  mesh.tet_vertices[f[k]].pos[2]);
            }
            double dis_2 = GEO::Geom::point_triangle_squared_distance(p, vs[0], vs[1], vs[2]);
            if (dis_2 <= eps) {
                is_inside = true;
                break;
            }
        }
        if (!is_inside) {
            cout << "check is_input_face_covered fail!!" << endl;
            pausee();
            break;
        }
    }

//    if (global_i > 2700) {
//        cout << global_i << " covered" << endl;
//        //output input/tet triangles in different colors
//        {
//            Eigen::MatrixXd V(t_fs.size() * 3, 3), C(t_fs.size() * 3, 3);
//            Eigen::MatrixXi F(t_fs.size(), 3);
//            for (int i = 0; i < t_fs.size(); i++) {
//                for (int j = 0; j < 3; j++) {
//                    V.row(i * 3 + j) = mesh.tet_vertices[t_fs[i][j]].pos;
//                    C.row(i * 3 + j) << 0, 0, 255;
//                }
//                F.row(i) << i * 3, i * 3 + 1, i * 3 + 2;
//            }
//            igl::writeOFF("covered_tet_fs_" + std::to_string(global_i) + ".off", V, F, C);
//        }
//        {
//            Eigen::MatrixXd V(1 * 3, 3), C(1 * 3, 3);
//            Eigen::MatrixXi F(1, 3);
//            int i = 0;
//            for (int j = 0; j < 3; j++) {
//                V.row(i * 3 + j) = f_vs[j];
//                C.row(i * 3 + j) << 255, 0, 0;
//            }
//            F.row(i) << i * 3, i * 3 + 1, i * 3 + 2;
//            igl::writeOFF("covered_input_f_" + std::to_string(global_i) + ".off", V, F, C);
//        }
//    }
}

void floatTetWild::check_edge_preservation(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
                             const std::vector<int> &sorted_f_ids, const Mesh &mesh,
                             std::vector<std::array<std::vector<int>, 4>> &track_surface_fs) {
    for (int i = sorted_f_ids.size() * 0.8; i < sorted_f_ids.size(); i++) {
        std::array<Vector3, 3> f_vs({{input_vertices[input_faces[sorted_f_ids[i]][0]],
                                             input_vertices[input_faces[sorted_f_ids[i]][1]],
                                             input_vertices[input_faces[sorted_f_ids[i]][2]]}});
        std::vector<std::array<int, 3>> t_fs;
        //todo
    }
}