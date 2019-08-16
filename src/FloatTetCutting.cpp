#include <floattetwild/FloatTetCutting.h>
#include <floattetwild/FloatTetCuttingCheck.h>
#include <floattetwild/FloatTetCuttingParallel.h>
#include <floattetwild/LocalOperations.h>
#include <floattetwild/Predicates.hpp>
#include <floattetwild/MeshIO.hpp>
#include <floattetwild/auto_table.hpp>
#include <floattetwild/Logger.hpp>
#include <floattetwild/intersections.h>


#include <igl/writeSTL.h>
#include <igl/Timer.h>

#ifdef USE_TBB
#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include <tbb/atomic.h>
#include <floattetwild/FloatTetCuttingParallel.h>
#include <tbb/concurrent_queue.h>
#endif

#include <bitset>
#include <numeric>
#include <unordered_map>

double find_tets_for_cut_time;
double tet_subdivision_time;
double update_opp_t_ids_time;
double insert_new_time;
double one_face_cut_other_time;

const int FAIL_CONFIG = 0;
const int FAIL_INVERSION = -1;
const int FAIL = -2;
const int SUCCESS = 1;

void floatTetWild::cutting(const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces, const std::vector<int> &input_tags,
        Mesh& mesh, std::vector<bool>& is_face_inserted, AABBWrapper& tree, bool is_again) {

    igl::Timer timer;
    find_tets_for_cut_time = 0;
    tet_subdivision_time = 0;
    update_opp_t_ids_time = 0;
    one_face_cut_other_time = 0;
    insert_new_time = 0;

    std::vector<int> current_inserted_cutting_fs;
    std::vector<std::array<std::vector<int>, 4>> cut_f_ids(mesh.tets.size());
    std::vector<bool> is_face_matched(input_faces.size(), false);
    if (!is_again) {
        timer.start();
        match_surface_fs(mesh, input_vertices, input_faces, is_face_inserted, cut_f_ids);
        is_face_matched = is_face_inserted;
        current_inserted_cutting_fs.reserve(std::count(is_face_inserted.begin(), is_face_inserted.end(), true));
        for (int i = 0; i < is_face_inserted.size(); i++) {
            if (is_face_inserted[i])
                current_inserted_cutting_fs.push_back(i);
        }
        timer.stop();
        logger().info("match_surface_fs {}s", timer.getElapsedTimeInSec());
    }
    int cnt_fail = std::count(is_face_inserted.begin(), is_face_inserted.end(), false);
    logger().info("failing before cutting {}/{}", cnt_fail, is_face_inserted.size());


    timer.start();
#ifdef USE_TBB
    int N = std::min((int)mesh.params.num_threads, 4);//tbb::task_scheduler_init::default_num_threads()
    if(cnt_fail < 5000)
        N = std::min((int)mesh.params.num_threads, 2);
    if(cnt_fail < 500)
        N = 1;

    if(N <= 1){ // || cnt_fail < 500
        for (int f_id = 0; f_id < input_faces.size(); f_id++) {
            if(is_face_inserted[f_id])
                continue;
            if(insert_one(input_vertices, input_faces, mesh, cut_f_ids, f_id, input_faces[f_id][0], false, is_again)){
                current_inserted_cutting_fs.push_back(f_id);
                is_face_inserted[f_id] = true;
            }
        }
        std::sort(current_inserted_cutting_fs.begin(), current_inserted_cutting_fs.end());
    } else {
        std::vector<std::vector<int>> partition_t_ids;
        mesh.partition(N, partition_t_ids);

        if (!is_again) {
            // mesh.partition(N, partition_t_ids);

            for (int i = 0; i < partition_t_ids.size(); i++) {
                for (int t_id: partition_t_ids[i])
                    mesh.tets[t_id].scalar = i;
            }

            int cnt = 0;
            for (int f_id = 0; f_id < input_faces.size(); f_id++) {
                if (is_face_inserted[f_id])
                    continue;

                int color = -1;
                for (int j = 0; j < 3; j++) {
                    for (int t_id:mesh.tet_vertices[input_faces[f_id][j]].conn_tets) {
                        if (color < 0)
                            color = mesh.tets[t_id].scalar;
                        else if (color != mesh.tets[t_id].scalar) {
                            color = -2;
                            break;
                        }
                    }
                }
                if (color != -2)
                    continue;

                // cnt++;
                if (insert_one(input_vertices, input_faces, mesh, cut_f_ids, f_id, input_faces[f_id][0], false, is_again)) {
                    current_inserted_cutting_fs.push_back(f_id);
                    is_face_inserted[f_id] = true;
                    cnt++;
                    if (cnt % 1000 == 0){
                        cout << cnt << endl;
                        cout<<"#v = " << mesh.tet_vertices.size()<<std::endl;
                        cout<<"#t = " << mesh.tets.size()<<std::endl;
                    }
                }
            }
            // cout<<"cnt = "<<cnt<<endl;
            logger().info("failing after 1st serial cutting {}/{}",
                    std::count(is_face_inserted.begin(), is_face_inserted.end(), false),
                    is_face_inserted.size());

            //update partition_t_ids
            int old_partition_t_ids_size = partition_t_ids.size();
            partition_t_ids.clear();
            partition_t_ids.resize(old_partition_t_ids_size);
            for (int t_id = 0; t_id < mesh.tets.size(); t_id++) {
                if (mesh.tets[t_id].scalar >= old_partition_t_ids_size || mesh.tets[t_id].scalar < 0) {
                    cout << "t.scalar>=old_partition_t_ids_size || t.scalar<0" << endl;
                    pausee();
                }
                partition_t_ids[mesh.tets[t_id].scalar].push_back(t_id);
            }

            //clear t scalars
            for (auto &t: mesh.tets)
                t.scalar = 0;
        }

        check(mesh);////
        cout << "check done" << endl;
        check_cut_f_ids(input_vertices, input_faces, mesh, cut_f_ids);////
        cout << "check_cut_f_ids done" << endl;

        // partition_t_ids.clear();
        // mesh.partition(N, partition_t_ids);
        std::vector<int> map_input_to_merged_mesh;

        {

        std::vector<Mesh> sub_meshes;
        std::vector<std::vector<std::array<std::vector<int>, 4>>> sub_cut_f_ids;
        // tbb::concurrent_vector<Mesh> sub_meshes;
        // tbb::concurrent_vector<std::vector<std::array<std::vector<int>, 4>>> sub_cut_f_ids;
        std::vector<std::array<int, 2>> map_input_to_partition(input_vertices.size(), {{-1, -1}});
        std::vector<std::vector<std::array<int, 2>>> v_p_ids;
        partition_mesh(partition_t_ids, mesh, sub_meshes,
                       cut_f_ids, sub_cut_f_ids,
                       map_input_to_partition, v_p_ids, is_again);

        std::vector<std::vector<int>> inserted_f_ids(sub_meshes.size());
        // for(int n=0;n<sub_meshes.size();n++) {
        tbb::parallel_for(size_t(0), size_t(sub_meshes.size()), [&](size_t n) {
            int cnt_insert = 0;
            int cnt_tried = 0;
            logger().info("partition {}", n);
            logger().info("#t = {}", sub_meshes[n].tets.size());

            for (int f_id = 0; f_id < input_faces.size(); f_id++) {
                if (is_face_inserted[f_id])
                    continue;
                if (map_input_to_partition[input_faces[f_id][0]][0] != n
                    || map_input_to_partition[input_faces[f_id][1]][0] != n
                    || map_input_to_partition[input_faces[f_id][2]][0] != n)
                    continue;
                if (insert_one(input_vertices, input_faces, sub_meshes[n], sub_cut_f_ids[n],
                               f_id, map_input_to_partition[input_faces[f_id][0]][1], true, is_again)) {
                    inserted_f_ids[n].push_back(f_id);
                    cnt_insert++;
                }
                cnt_tried++;
            }

            logger().info("insert #f: {}", cnt_insert);
            logger().info("tried #f: {}", cnt_tried);
        // }
        });

        for (const auto &f_ids:inserted_f_ids) {
            current_inserted_cutting_fs.insert(current_inserted_cutting_fs.end(), f_ids.begin(), f_ids.end());
            for (int f_id:f_ids)
                is_face_inserted[f_id] = true;
        }

        logger().info("failing after parallel cutting {}/{}",
                      std::count(is_face_inserted.begin(), is_face_inserted.end(), false),
                      is_face_inserted.size());

        map_input_to_merged_mesh.clear();
        map_input_to_merged_mesh.resize(input_vertices.size(), -1);
        //std::vector<int> map_input_to_merged_mesh(input_vertices.size(), -1);
        merge_meshes(sub_meshes, mesh, sub_cut_f_ids, cut_f_ids, map_input_to_partition, map_input_to_merged_mesh,
                     v_p_ids);
        }

        int cnt = 0;
        for (int f_id = 0; f_id < input_faces.size(); f_id++) {
            if (is_face_inserted[f_id])
                continue;
            if (insert_one(input_vertices, input_faces, mesh, cut_f_ids,
                           f_id, map_input_to_merged_mesh[input_faces[f_id][0]], false, is_again)) {
                current_inserted_cutting_fs.push_back(f_id);
                is_face_inserted[f_id] = true;
                cnt++;
                if(cnt%1000 == 0){
                    cout<<cnt<<endl;
                    cout<<"#v = " << mesh.tet_vertices.size()<<std::endl;
                    cout<<"#t = " << mesh.tets.size()<<std::endl;
                }
            }
        }

        vector_unique(current_inserted_cutting_fs);
    }
#else
    for (int f_id = 0; f_id < input_faces.size(); f_id++) {
        if(is_face_inserted[f_id])
            continue;
        if(insert_one(input_vertices, input_faces, mesh, cut_f_ids, f_id, input_faces[f_id][0], false, is_again)){
            current_inserted_cutting_fs.push_back(f_id);
            is_face_inserted[f_id] = true;
        }
    }
    std::sort(current_inserted_cutting_fs.begin(), current_inserted_cutting_fs.end());
#endif

    timer.stop();
//    logger().info("find tets for cut {}s", find_time);
    logger().info("cutting tets {}s", timer.getElapsedTimeInSec());
    logger().info("find_tets_for_cut_time {}s", find_tets_for_cut_time);
    logger().info("tet_subdivision_time {}s", tet_subdivision_time);
    logger().info("insert_new_time {}s", insert_new_time);
    logger().info("update_opp_t_ids_time {}s", update_opp_t_ids_time);
    logger().info("one_face_cut_other_time {}s", one_face_cut_other_time);

    cnt_fail = std::count(is_face_inserted.begin(), is_face_inserted.end(), false);
    logger().info("failing after {}/{}", cnt_fail, is_face_inserted.size());
    logger().info("#v = {}", mesh.tet_vertices.size());
    logger().info("#t = {}", mesh.tets.size());

//    std::ofstream fout;
//    fout.open(mesh.params.log_filename);
//    fout << is_face_inserted.size() << " " << cnt_fail << " " << cnt_fail_config << endl;
//    fout.close();

    check(mesh);////
    cout << "check done" << endl;
    check_cut_f_ids(input_vertices, input_faces, mesh, cut_f_ids);////
    cout << "check_cut_f_ids done" << endl;

    find_tets_for_cut_time = 0;
    tet_subdivision_time = 0;
    insert_new_time = 0;
    update_opp_t_ids_time = 0;
    one_face_cut_other_time = 0;

//    plot_cover_for_trif(COMMON_INPUT_FOR_CHECK, 930);
//    plot_cover_for_trif(COMMON_INPUT_FOR_CHECK, 960);

    timer.start();
    std::vector<bool> is_boundary_preserved(input_faces.size(), true);

//    for(int i=0;i<is_face_inserted.size();i++) {//todo: write better
//        if (is_face_inserted[i])
//            current_inserted_cutting_fs.push_back(i);
//    }
//    vector_unique(current_inserted_cutting_fs);

    preserve_cutting_open_boundary(input_vertices, input_faces, current_inserted_cutting_fs, cut_f_ids, mesh, tree,
            is_face_matched, is_boundary_preserved, is_again);
    for(int i=0;i<is_face_inserted.size();i++) {
        if (!is_boundary_preserved[i])
            is_face_inserted[i] = false;
    }
    logger().info("preserving open boundary {}s", timer.getElapsedTimeInSec());
    logger().info("find_tets_for_cut_time {}s", find_tets_for_cut_time);
    logger().info("tet_subdivision_time {}s", tet_subdivision_time);
    logger().info("insert_new_time {}s", insert_new_time);
    logger().info("update_opp_t_ids_time {}s", update_opp_t_ids_time);
    logger().info("one_face_cut_other_time {}s", one_face_cut_other_time);
    cnt_fail = std::count(is_face_inserted.begin(), is_face_inserted.end(), false);
    if(cnt_fail == 0)
        mesh.is_input_all_inserted = true;
    logger().info("failing after preserving boundary {}/{}", cnt_fail, is_face_inserted.size());

    update_tmp_b_tree(input_vertices, input_faces, is_face_inserted, tree);

    check(mesh);////
    cout << "check_cut_f_ids 1" << endl;
    check_cut_f_ids(input_vertices, input_faces, mesh, cut_f_ids);////
    track_surface(mesh, input_vertices, input_faces, input_tags, cut_f_ids, tree, is_boundary_preserved);
    check_is_surface_fs(mesh);////

    ///////test
//    plot_cover_for_trif(COMMON_INPUT_FOR_CHECK, 930);
//    plot_cover_for_trif(COMMON_INPUT_FOR_CHECK, 960);
//    plot_cover_for_trif(COMMON_INPUT_FOR_CHECK, 962);
//
//    plot_cover_for_tetf(COMMON_INPUT_FOR_CHECK, 5248, 0);
//    plot_cover_for_tetf(COMMON_INPUT_FOR_CHECK, 5248, 1);
//    plot_cover_for_tetf(COMMON_INPUT_FOR_CHECK, 5248, 2);
//    plot_cover_for_tetf(COMMON_INPUT_FOR_CHECK, 5248, 3);
//
//    plot_cover_for_tetf(COMMON_INPUT_FOR_CHECK, 119574, 0);
//    plot_cover_for_tetf(COMMON_INPUT_FOR_CHECK, 119574, 1);
//    plot_cover_for_tetf(COMMON_INPUT_FOR_CHECK, 119574, 2);
//    plot_cover_for_tetf(COMMON_INPUT_FOR_CHECK, 119574, 3);
    ///////test


    ///////test
//    int f_id = 930;
//    std::vector<std::unordered_set<int>> conn_fs(input_vertices.size());
//    for(int i=0;i<input_faces.size();i++) {
//        for (int j = 0; j < 3; j++)
//            conn_fs[input_faces[i][j]].insert(i);
//    }
//    std::array<int, 4> conn_f;
//    conn_f[0] = f_id;
//    for(int j=0;j<3;j++) {
//        std::vector<int> pair;
//        set_intersection(conn_fs[input_faces[f_id][j]], conn_fs[input_faces[f_id][(j + 1) % 3]], pair);
//        cout << pair.size() << endl;
//        if (pair.size() < 2)
//            conn_f[j + 1] = -1;
//        else
//            conn_f[j + 1] = pair[0] == f_id ? pair[1] : pair[0];
//    }
//    for(int f_id:conn_f) {
//        if (f_id < 0)
//            continue;
//        plot_cover_for_trif(COMMON_INPUT_FOR_CHECK, f_id);
//    }
    ///////test

    if (cnt_fail > 0)
        mesh.is_input_all_inserted = false;
    else
        mesh.is_input_all_inserted = true;

    timer.start();
    if (!is_again) {
        for (auto &t: mesh.tets) {
            if (t.is_removed)
                continue;
            t.quality = get_quality(mesh, t);
        }
    }
    logger().info("recomputing quality {}s", timer.getElapsedTimeInSec());


//    ////////////////////
//
//    Eigen::Matrix<Scalar, Eigen::Dynamic, 3> V(input_vertices.size(), 3);
//    Eigen::MatrixXi F(/*input_faces.size() - */cnt_fail, 3);
//    for (int i = 0; i < input_vertices.size(); ++i)
//        V.row(i) = input_vertices[i];
//    int cnt = 0;
//    for (int i = 0; i < input_faces.size(); ++i) {
//        if (!is_face_inserted[i]) {
//            F.row(cnt) = input_faces[i];
//            cnt++;
//        }
//    }
//    igl::writeSTL("/Users/yixinhu/Downloads/inserted.stl", V, F);
////    igl::writeSTL("/Users/HuYixin/Downloads/inserted.stl", V, F);
}

bool floatTetWild::insert_one(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
                Mesh &sub_mesh, std::vector<std::array<std::vector<int>, 4>>& sub_cut_f_ids,
                int f_id, int seed_v_id, bool is_parallel, bool is_again){
    igl::Timer igl_timer;
//    cout<<"inserting "<<f_id<<endl;

    //find tets for cut
    std::vector<int> oris(sub_mesh.tet_vertices.size(), Predicates::ORI_UNKNOWN);
    std::vector<int> intersection_results;
    igl_timer.start();
    if (is_again)
        find_tets_for_cut_again(input_vertices, input_faces, f_id, sub_mesh, intersection_results, oris);
    else
        find_tets_for_cut(input_vertices, input_faces, f_id, seed_v_id, sub_mesh, intersection_results, oris);
    find_tets_for_cut_time += igl_timer.getElapsedTime();
    if (intersection_results.empty()) {
        if(!is_parallel) {
            cout << "intersection_results.empty() " << f_id << " ";//could happen when parallel
            cout << get_area(input_vertices[input_faces[f_id][0]], input_vertices[input_faces[f_id][1]],
                             input_vertices[input_faces[f_id][2]]) << endl;
//            check(sub_mesh);
//            cout<<"check done"<<endl;

//            std::vector<std::array<int, 3>> fs;
//            for(int t_id: sub_mesh.tet_vertices[7397].conn_tets){
//                cout<<Predicates::orient_3d_volume(sub_mesh.tet_vertices[sub_mesh.tets[t_id][0]].pos,
//                                                   sub_mesh.tet_vertices[sub_mesh.tets[t_id][1]].pos,
//                                                   sub_mesh.tet_vertices[sub_mesh.tets[t_id][2]].pos,
//                                                   sub_mesh.tet_vertices[sub_mesh.tets[t_id][3]].pos)<<endl;
//                for(int j=0;j<4;j++) {
//                    if(sub_mesh.tets[t_id][j]==7397){
//                        fs.push_back({{sub_mesh.tets[t_id][(j+1)%4],
//                                              sub_mesh.tets[t_id][(j+2)%4],
//                                              sub_mesh.tets[t_id][(j+3)%4]}});
//                        break;
//                    }
////                    int opp_t_id = get_opp_t_id(sub_mesh, t_id, j);
////                    if(sub_mesh.tet_vertices[7397].conn_tets.find(opp_t_id)==sub_mesh.tet_vertices[7397].conn_tets.end()){
////                        b_v_ids.push_back(sub_mesh.tets[t_id][(j+1)%4]);
////                        b_v_ids.push_back(sub_mesh.tets[t_id][(j+2)%4]);
////                        b_v_ids.push_back(sub_mesh.tets[t_id][(j+3)%4]);
////                    }
//                }
//            }
////            vector_unique(b_v_ids);
////            for(auto& i: b_v_ids){
////                cout<<(sub_mesh.tet_vertices[7397].pos-sub_mesh.tet_vertices[i].pos).norm()<<endl;
////            }
////            cout<<(std::find(b_v_ids.begin(), b_v_ids.end(), 7397) == b_v_ids.end())<<endl;
//            std::vector<std::unordered_set<int>> conn_fs(sub_mesh.tet_vertices.size());
//            for(int i=0;i<fs.size();i++){
//                for(int j=0;j<3;j++) {
//                    conn_fs[fs[i][j]].insert(i);
//                }
//                cout<<"f "<<i<<": "<<fs[i][0]<<" "<<fs[i][1]<<" "<<fs[i][2]<<": ";
//                cout<<get_area(sub_mesh.tet_vertices[fs[i][0]].pos,
//                               sub_mesh.tet_vertices[fs[i][1]].pos,
//                               sub_mesh.tet_vertices[fs[i][2]].pos)<<endl;
//            }
//            for(auto& f: fs){
//                for(int j=0;j<3;j++) {
//                    std::vector<int> tmp;
//                    set_intersection(conn_fs[f[j]], conn_fs[f[(j+1)%3]], tmp);
//                    if(tmp.size()!=2){
//                        cout<<tmp.size()<<endl;
//                        pausee();
//                    } else {
//                        cout<<f[j]<<" "<<f[(j+1)%3]<<": "<<tmp[0]<<" "<<tmp[1]<<endl;
//                    }
//                }
//            }
//            cout<<"check done"<<endl;

        }
//        pausee();
        return false;
    }

    int result = one_face_cut(input_vertices, input_faces, f_id, sub_mesh, intersection_results, sub_cut_f_ids, oris, is_parallel,
                              is_again);

    if (result == SUCCESS) {
        return true;
    } else {
//        if(is_again) {
//            if (result == FAIL_CONFIG)
//                cout << "FAIL_CONFIG" << endl;
//            else if (result == FAIL)
//                cout << "FAIL" << endl;
//            else
//                cout << "FAIL_INVERSION" << endl;
//        }
        return false;
    }
}

void floatTetWild::match_surface_fs(const Mesh &mesh, const std::vector<Vector3> &input_vertices,
                                    const std::vector<Vector3i> &input_faces,
                                    std::vector<bool> &is_face_inserted,
                                    std::vector<std::array<std::vector<int>, 4>>& cut_f_ids) {
    auto comp = [](const std::array<int, 4> &a, const std::array<int, 4> &b) {
        return std::tuple<int, int, int>(a[0], a[1], a[2]) < std::tuple<int, int, int>(b[0], b[1], b[2]);
    };

    std::vector<std::array<int, 4>> input_fs(input_faces.size());
    for (int i = 0; i < input_faces.size(); i++) {
        input_fs[i] = {{input_faces[i][0], input_faces[i][1], input_faces[i][2], i}};
        std::sort(input_fs[i].begin(), input_fs[i].begin() + 3);
    }
    std::sort(input_fs.begin(), input_fs.end(), comp);

    for (int i = 0; i < mesh.tets.size(); i++) {
        auto &t = mesh.tets[i];
        for (int j = 0; j < 4; j++) {
            std::array<int, 3> f = {{t[mod4(j + 1)], t[mod4(j + 2)], t[mod4(j + 3)]}};
            std::sort(f.begin(), f.end());
            auto bounds = std::equal_range(input_fs.begin(), input_fs.end(),
                                           std::array<int, 4>({{f[0], f[1], f[2], -1}}),
                                           comp);
            for (auto it = bounds.first; it != bounds.second; ++it) {
                int f_id = (*it)[3];
                is_face_inserted[f_id] = true;
                cut_f_ids[i][j].push_back(f_id);
            }
        }
    }
}

void floatTetWild::find_tets_for_cut_again(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
                                    int f_id, Mesh &mesh, std::vector<int> &intersection_results, std::vector<int>& oris) {
    auto &tet_vertices = mesh.tet_vertices;
    auto &tets = mesh.tets;

    const Vector3 &p1 = input_vertices[input_faces[f_id][0]];
    const Vector3 &p2 = input_vertices[input_faces[f_id][1]];
    const Vector3 &p3 = input_vertices[input_faces[f_id][2]];

    Vector3 min_f, max_f;
    get_bbox_face(p1, p2, p3, min_f, max_f);

//    if(f_id == 6214) {
//        cout << std::setprecision(17);
//        cout << p1.transpose() << endl;
//        cout << p2.transpose() << endl;
//        cout << p3.transpose() << endl;
//    }

//    std::vector<bool> is_visited(tets.size(), false);

#ifdef USE_TBB
    tbb::concurrent_vector<int> t_id_queue;
    tbb::parallel_for(size_t(0), size_t(tets.size()), [&](size_t t_id) {
        if (tets[t_id].is_removed)
            return;

        Vector3 min_t, max_t;
        get_bbox_tet(tet_vertices[tets[t_id][0]].pos, tet_vertices[tets[t_id][1]].pos, tet_vertices[tets[t_id][2]].pos,
                     tet_vertices[tets[t_id][3]].pos, min_t, max_t);
        if (!is_bbox_intersected(min_f, max_f, min_t, max_t))
            return;

        t_id_queue.push_back(t_id);
    });
#else
    std::vector<int> t_id_queue;
    for (int t_id = 0; t_id < tets.size(); t_id++) {
        if (tets[t_id].is_removed)
            continue;

//        if(is_point_inside_tet(p1, tet_vertices[tets[t_id][0]].pos, tet_vertices[tets[t_id][1]].pos,
//                tet_vertices[tets[t_id][2]].pos, tet_vertices[tets[t_id][3]].pos)) {
//            t_id_queue.push(t_id);
//            cout<<"t_id = "<<t_id<<endl;
//            break;
//        }

        Vector3 min_t, max_t;
        get_bbox_tet(tet_vertices[tets[t_id][0]].pos, tet_vertices[tets[t_id][1]].pos, tet_vertices[tets[t_id][2]].pos,
                     tet_vertices[tets[t_id][3]].pos, min_t, max_t);
        if (!is_bbox_intersected(min_f, max_f, min_t, max_t))
            continue;

        t_id_queue.push_back(t_id);
    }
#endif

//    bool is_seed_found = false;
    for(int t_id:t_id_queue){
//        int t_id = t_id_queue.front();
//        t_id_queue.pop();
//        if (is_visited[t_id])
//            continue;
//        is_visited[t_id] = true;

        //compare bbox
//        Vector3 min_t, max_t;
//        get_bbox_tet(tet_vertices[tets[t_id][0]].pos, tet_vertices[tets[t_id][1]].pos, tet_vertices[tets[t_id][2]].pos,
//                     tet_vertices[tets[t_id][3]].pos, min_t, max_t);
//        if (!is_bbox_intersected(min_f, max_f, min_t, max_t))
//            continue;

        //check if the tet is on the same side of the plane
        int cnt_pos = 0;
        int cnt_neg = 0;
        int cnt_zero = 0;
        for (int j = 0; j < 4; j++) {
            if (oris[tets[t_id][j]] == Predicates::ORI_UNKNOWN) {
                oris[tets[t_id][j]] = Predicates::orient_3d(p1, p2, p3, tet_vertices[tets[t_id][j]].pos);
            }
            if (oris[tets[t_id][j]] == Predicates::ORI_POSITIVE)
                cnt_pos++;
            else if (oris[tets[t_id][j]] == Predicates::ORI_NEGATIVE)
                cnt_neg++;
            else
                cnt_zero++;
        }
        if (cnt_neg == 4 || cnt_pos == 4 || ((cnt_neg == 0 || cnt_pos == 0) && cnt_zero == 1))
            continue;

        //check if the triangle is contained inside the tet
        if (is_tri_inside_tet({{p1, p2, p3}},
                              tet_vertices[tets[t_id][0]].pos, tet_vertices[tets[t_id][1]].pos,
                              tet_vertices[tets[t_id][2]].pos, tet_vertices[tets[t_id][3]].pos)) {
            intersection_results = {t_id};
            break;
        }

        //check triangle intersection
        for (int j = 0; j < 4; j++) {
            int opp_t_id = tets[t_id].opp_t_ids[j];
            if(opp_t_id == OPP_T_ID_UNKNOWN){
                opp_t_id = get_opp_t_id(mesh, t_id, j);
                if(opp_t_id!=OPP_T_ID_BOUNDARY) {
                    int k = mesh.tets[opp_t_id].find_opp(tets[t_id][mod4(j + 1)], tets[t_id][mod4(j + 2)],
                                                         tets[t_id][mod4(j + 3)]);
                    mesh.tets[opp_t_id].opp_t_ids[k] = t_id;
                }
            }

//            if (opp_t_id < 0 || is_visited[opp_t_id])
//            if (opp_t_id < 0)
//                continue;

            int cnt_pos_f = 0;
            int cnt_neg_f = 0;
            int cnt_zero_f = 0;
            for (int k = 0; k < 3; k++) {
                if (oris[tets[t_id][mod4(j + 1 + k)]] == Predicates::ORI_POSITIVE)
                    cnt_pos_f++;
                else if (oris[tets[t_id][mod4(j + 1 + k)]] == Predicates::ORI_NEGATIVE)
                    cnt_neg_f++;
                else
                    cnt_zero_f++;
            }

            if (cnt_zero_f == 3) {
                int result = is_tri_tri_cutted_hint(p1, p2, p3, tet_vertices[tets[t_id][mod4(j + 1)]].pos,
                                                    tet_vertices[tets[t_id][mod4(j + 2)]].pos,
                                                    tet_vertices[tets[t_id][mod4(j + 3)]].pos, CUT_COPLANAR);
                if (result == CUT_COPLANAR) {
                    intersection_results.push_back(t_id);
//                    if (opp_t_id >= 0)
//                        intersection_results.push_back(opp_t_id);
                }
            } else if (oris[tets[t_id][mod4(j + 1)]] == Predicates::ORI_ZERO
                       && oris[tets[t_id][mod4(j + 2)]] == Predicates::ORI_ZERO) {
//                if (is_tri_tri_cutted_hint(p1, p2, p3, tet_vertices[tets[t_id][(j + 1) % 4]].pos,
//                                           tet_vertices[tets[t_id][(j + 2) % 4]].pos,
//                                           tet_vertices[tets[t_id][(j + 3) % 4]].pos, CUT_EDGE_0) == CUT_EDGE_0) {
//                    std::vector<int> tmp;
//                    set_intersection(tet_vertices[tets[t_id][(j + 1) % 4]].conn_tets,
//                                     tet_vertices[tets[t_id][(j + 2) % 4]].conn_tets, tmp);
//                    for (int n_t_id: tmp)
//                        if (!is_visited[n_t_id])
//                            t_id_queue.push(n_t_id);
//                }
            } else if (oris[tets[t_id][mod4(j + 2)]] == Predicates::ORI_ZERO
                       && oris[tets[t_id][mod4(j + 3)]] == Predicates::ORI_ZERO) {
//                if (is_tri_tri_cutted_hint(p1, p2, p3, tet_vertices[tets[t_id][(j + 1) % 4]].pos,
//                                           tet_vertices[tets[t_id][(j + 2) % 4]].pos,
//                                           tet_vertices[tets[t_id][(j + 3) % 4]].pos, CUT_EDGE_1) == CUT_EDGE_1) {
//                    std::vector<int> tmp;
//                    set_intersection(tet_vertices[tets[t_id][(j + 2) % 4]].conn_tets,
//                                     tet_vertices[tets[t_id][(j + 3) % 4]].conn_tets, tmp);
//                    for (int n_t_id: tmp)
//                        if (!is_visited[n_t_id])
//                            t_id_queue.push(n_t_id);
//                }
            } else if (oris[tets[t_id][mod4(j + 3)]] == Predicates::ORI_ZERO
                       && oris[tets[t_id][mod4(j + 1)]] == Predicates::ORI_ZERO) {
//                if (is_tri_tri_cutted_hint(p1, p2, p3, tet_vertices[tets[t_id][(j + 1) % 4]].pos,
//                                           tet_vertices[tets[t_id][(j + 2) % 4]].pos,
//                                           tet_vertices[tets[t_id][(j + 3) % 4]].pos, CUT_EDGE_2) == CUT_EDGE_2) {
//                    std::vector<int> tmp;
//                    set_intersection(tet_vertices[tets[t_id][(j + 3) % 4]].conn_tets,
//                                     tet_vertices[tets[t_id][(j + 1) % 4]].conn_tets, tmp);
//                    for (int n_t_id: tmp)
//                        if (!is_visited[n_t_id])
//                            t_id_queue.push(n_t_id);
//                }
            } else if (cnt_pos_f > 0 && cnt_neg_f > 0) {
                int result;
                result = is_tri_tri_cutted_hint(p1, p2, p3, tet_vertices[tets[t_id][mod4(j + 1)]].pos,
                                                tet_vertices[tets[t_id][mod4(j + 2)]].pos,
                                                tet_vertices[tets[t_id][mod4(j + 3)]].pos, CUT_FACE);
                if (result == CUT_FACE) {
                    intersection_results.push_back(t_id);
//                    if (opp_t_id >= 0)
//                        intersection_results.push_back(opp_t_id);
                }
            }
        }
    }

    vector_unique(intersection_results);

//    if(f_id == 6214) {
//        for (int t_id:intersection_results) {
//            for (int j = 0; j < 4; j++) {
//                cout << tet_vertices[tets[t_id][(j + 1) % 4]].pos.transpose() << endl;
//                cout << tet_vertices[tets[t_id][(j + 2) % 4]].pos.transpose() << endl;
//                cout << tet_vertices[tets[t_id][(j + 3) % 4]].pos.transpose() << endl;
//            }
//        }
//        pausee();
//    }
}

void floatTetWild::find_tets_for_cut(const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces, int f_id,
                                    int seed_v_id,
                                    Mesh& mesh, std::vector<int>& intersection_results, std::vector<int>& oris) {
    auto &tet_vertices = mesh.tet_vertices;
    auto &tets = mesh.tets;

    std::vector<bool> is_visited(tets.size(), false);

    const Vector3 &p1 = input_vertices[input_faces[f_id][0]];
    const Vector3 &p2 = input_vertices[input_faces[f_id][1]];
    const Vector3 &p3 = input_vertices[input_faces[f_id][2]];

    int test_f_id = -1;
    if(f_id == test_f_id) {
        cout << std::setprecision(17);
        cout << p1.transpose() << endl;
        cout << p2.transpose() << endl;
        cout << p3.transpose() << endl;
    }

    std::queue<int> t_id_queue;
//    for (int t_id:tet_vertices[input_faces[f_id][0]].conn_tets)
    for (int t_id:tet_vertices[seed_v_id].conn_tets)
        t_id_queue.push(t_id);
    while (!t_id_queue.empty()) {
        int t_id = t_id_queue.front();
        t_id_queue.pop();

        if (is_visited[t_id])
            continue;
        is_visited[t_id] = true;

        //check if the tet is on the same side of the plane
        int cnt_pos = 0;
        int cnt_neg = 0;
        int cnt_zero = 0;
        for (int j = 0; j < 4; j++) {
            if (oris[tets[t_id][j]] == Predicates::ORI_UNKNOWN) {
                oris[tets[t_id][j]] = Predicates::orient_3d(p1, p2, p3, tet_vertices[tets[t_id][j]].pos);
            }
            if (oris[tets[t_id][j]] == Predicates::ORI_POSITIVE)
                cnt_pos++;
            else if (oris[tets[t_id][j]] == Predicates::ORI_NEGATIVE)
                cnt_neg++;
            else
                cnt_zero++;
        }
        if (cnt_neg == 4 || cnt_pos == 4 || ((cnt_neg == 0 || cnt_pos == 0) && cnt_zero == 1))
            continue;

        //check triangle intersection
        for (int j = 0; j < 4; j++) {
            int opp_t_id = tets[t_id].opp_t_ids[j];
            if(opp_t_id == OPP_T_ID_UNKNOWN){
                opp_t_id = get_opp_t_id(mesh, t_id, j);
                if(opp_t_id!=OPP_T_ID_BOUNDARY) {
                    int k = mesh.tets[opp_t_id].find_opp(tets[t_id][mod4(j + 1)], tets[t_id][mod4(j + 2)],
                                                         tets[t_id][mod4(j + 3)]);
                    mesh.tets[opp_t_id].opp_t_ids[k] = t_id;
                }
            }

            if (opp_t_id >=0 && is_visited[opp_t_id])
                continue;

            int cnt_pos_f = 0;
            int cnt_neg_f = 0;
            int cnt_zero_f = 0;
            for (int k = 0; k < 3; k++) {
                if (oris[tets[t_id][mod4(j + 1 + k)]] == Predicates::ORI_POSITIVE)
                    cnt_pos_f++;
                else if (oris[tets[t_id][mod4(j + 1 + k)]] == Predicates::ORI_NEGATIVE)
                    cnt_neg_f++;
                else
                    cnt_zero_f++;
            }

            const auto& tp1 = tet_vertices[tets[t_id][mod4(j + 1)]].pos;
            const auto& tp2 = tet_vertices[tets[t_id][mod4(j + 2)]].pos;
            const auto& tp3 = tet_vertices[tets[t_id][mod4(j + 3)]].pos;

            if (cnt_zero_f == 3) {
                int result = is_tri_tri_cutted_hint(p1, p2, p3, tp1, tp2, tp3, CUT_COPLANAR);
                if (result == CUT_COPLANAR) {
                    intersection_results.push_back(t_id);
                    if(opp_t_id >= 0)
                        intersection_results.push_back(opp_t_id);
//                    t_id_queue.push(opp_t_id);

                    std::vector<int> new_t_ids;
                    for(int k=0;k<3;k++){
                        for(int n_t_id: tet_vertices[tets[t_id][mod4(j + 1+k)]].conn_tets)
                            if (!is_visited[n_t_id])
                                new_t_ids.push_back(n_t_id);
                    }
                    vector_unique(new_t_ids);
                    for(int n_t_id:new_t_ids)
                        t_id_queue.push(n_t_id);
                }
            } else if (oris[tets[t_id][mod4(j + 1)]] == Predicates::ORI_ZERO
                       && oris[tets[t_id][mod4(j + 2)]] == Predicates::ORI_ZERO) {
                if (is_tri_tri_cutted_hint(p1, p2, p3, tp1, tp2, tp3, CUT_EDGE_0) == CUT_EDGE_0) {
//                    std::vector<int> tmp;
//                    set_intersection(tet_vertices[tets[t_id][mod4(j + 1)]].conn_tets,
//                                     tet_vertices[tets[t_id][mod4(j + 2)]].conn_tets, tmp);
//                    for (int n_t_id: tmp)
//                        if (!is_visited[n_t_id])
//                            t_id_queue.push(n_t_id);

                    std::vector<int> new_t_ids;
                    for(int k=0;k<2;k++){
                        for(int n_t_id: tet_vertices[tets[t_id][mod4(j + 1+k)]].conn_tets)
                            if (!is_visited[n_t_id])
                                new_t_ids.push_back(n_t_id);
                    }
                    vector_unique(new_t_ids);
                    for(int n_t_id:new_t_ids)
                        t_id_queue.push(n_t_id);
                }
            } else if (oris[tets[t_id][mod4(j + 2)]] == Predicates::ORI_ZERO
                       && oris[tets[t_id][mod4(j + 3)]] == Predicates::ORI_ZERO) {
                if (is_tri_tri_cutted_hint(p1, p2, p3, tp1, tp2, tp3, CUT_EDGE_1) == CUT_EDGE_1) {
//                    std::vector<int> tmp;
//                    set_intersection(tet_vertices[tets[t_id][mod4(j + 2)]].conn_tets,
//                                     tet_vertices[tets[t_id][mod4(j + 3)]].conn_tets, tmp);
//                    for (int n_t_id: tmp)
//                        if (!is_visited[n_t_id])
//                            t_id_queue.push(n_t_id);

                    std::vector<int> new_t_ids;
                    for(int k=0;k<2;k++){
                        for(int n_t_id: tet_vertices[tets[t_id][mod4(j + 2+k)]].conn_tets)
                            if (!is_visited[n_t_id])
                                new_t_ids.push_back(n_t_id);
                    }
                    vector_unique(new_t_ids);
                    for(int n_t_id:new_t_ids)
                        t_id_queue.push(n_t_id);
                }
            } else if (oris[tets[t_id][mod4(j + 3)]] == Predicates::ORI_ZERO
                       && oris[tets[t_id][mod4(j + 1)]] == Predicates::ORI_ZERO) {
                if (is_tri_tri_cutted_hint(p1, p2, p3, tp1, tp2, tp3, CUT_EDGE_2) == CUT_EDGE_2) {
//                    std::vector<int> tmp;
//                    set_intersection(tet_vertices[tets[t_id][mod4(j + 3)]].conn_tets,
//                                     tet_vertices[tets[t_id][mod4(j + 1)]].conn_tets, tmp);
//                    for (int n_t_id: tmp)
//                        if (!is_visited[n_t_id])
//                            t_id_queue.push(n_t_id);

                    std::vector<int> new_t_ids;
                    for(int k=0;k<2;k++){
                        for(int n_t_id: tet_vertices[tets[t_id][mod4(j + (3+k)%3)]].conn_tets)
                            if (!is_visited[n_t_id])
                                new_t_ids.push_back(n_t_id);
                    }
                    vector_unique(new_t_ids);
                    for(int n_t_id:new_t_ids)
                        t_id_queue.push(n_t_id);
                }
            } else if (cnt_pos_f > 0 && cnt_neg_f > 0) {
                if (is_tri_tri_cutted_hint(p1, p2, p3, tp1, tp2, tp3, CUT_FACE) == CUT_FACE) {
                    intersection_results.push_back(t_id);
                    if(opp_t_id >= 0)
                        intersection_results.push_back(opp_t_id);
//                    t_id_queue.push(opp_t_id);

                    std::vector<int> new_t_ids;//NOTE: you have to traverse like this!! otherwise cannot find all tets!
                    for(int k=0;k<3;k++){
                        for(int n_t_id: tet_vertices[tets[t_id][mod4(j + 1+k)]].conn_tets)
                            if (!is_visited[n_t_id])
                                new_t_ids.push_back(n_t_id);
                    }
                    vector_unique(new_t_ids);
                    for(int n_t_id:new_t_ids)
                        t_id_queue.push(n_t_id);
                }
            }
        }
    }

    vector_unique(intersection_results);

    if(f_id == test_f_id) {
        std::vector<int> check_t_ids(tet_vertices[seed_v_id].conn_tets.begin(),
                                     tet_vertices[seed_v_id].conn_tets.end());// = intersection_results;
        for (int t_id:check_t_ids) {
            for (int j = 0; j < 4; j++) {
                cout << tet_vertices[tets[t_id][(j + 1) % 4]].pos.transpose() << endl;
                cout << tet_vertices[tets[t_id][(j + 2) % 4]].pos.transpose() << endl;
                cout << tet_vertices[tets[t_id][(j + 3) % 4]].pos.transpose() << endl;
            }
        }

        cout<<"is_cutting_cross_partitions = "<<is_cutting_cross_partitions(mesh, check_t_ids)<<endl;

//        for (int t_id:check_t_ids) {
//            cout << t_id << " "
//                 << Predicates::orient_3d_volume(tet_vertices[tets[t_id][0]].pos, tet_vertices[tets[t_id][1]].pos,
//                                                 tet_vertices[tets[t_id][2]].pos, tet_vertices[tets[t_id][3]].pos)<<endl;
//            for(int j=0;j<4;j++)
//                cout<<tets[t_id][j]<<" ";
//            cout<<endl;
//
//            for(int j=0;j<4;j++)
//                cout<<oris[tets[t_id][j]]<<" ";
//            cout<<endl;
//
////            for (int j = 0; j < 4; j++) {
////                auto& tp1 = tet_vertices[tets[t_id][mod4(j + 1)]].pos;
////                auto& tp2 = tet_vertices[tets[t_id][mod4(j + 2)]].pos;
////                auto& tp3 = tet_vertices[tets[t_id][mod4(j + 3)]].pos;
////                cout<<is_tri_tri_cutted_hint(p1, p2, p3, tp1, tp2, tp3, CUT_FACE)<<endl;
////            }
//        }
        pausee();
    }
}


void floatTetWild::snapping(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces, int f_id,
              Mesh &mesh,
              std::vector<int> &intersection_results, std::vector<int> &oris, std::vector<Scalar>& signed_dist) {
    auto &tet_vertices = mesh.tet_vertices;
    auto &tets = mesh.tets;

    const Vector3 &a = input_vertices[input_faces[f_id][0]];
    const Vector3 &b = input_vertices[input_faces[f_id][1]];
    const Vector3 &c = input_vertices[input_faces[f_id][2]];
    Vector3 n = ((b - c).cross(a - c)).normalized();

//    Vector3 min_f, max_f;
//    get_bbox_face(a, b, c, min_f, max_f);

    ////get all oris for vs
    std::vector<int> v_ids;
    v_ids.reserve(intersection_results.size() * 4);
    for (auto &t_id:intersection_results) {
        for (int j = 0; j < 4; j++)
            v_ids.push_back(tets[t_id][j]);
    }
    vector_unique(v_ids);

    signed_dist.resize(tet_vertices.size(), INT_MAX);
    std::fill(oris.begin(), oris.end(), Predicates::ORI_UNKNOWN);
    for (int i:v_ids) {
        if (signed_dist[i] != INT_MAX)
            continue;
        oris[i] = Predicates::orient_3d(a, b, c, tet_vertices[i].pos);
        if (oris[i] == Predicates::ORI_ZERO) {
            signed_dist[i] = 0;
            continue;
        }
        signed_dist[i] = n.dot(tet_vertices[i].pos - a);
    }

    ////compute edges intersections
    std::vector<std::array<int, 2>> edges;
    edges.reserve(intersection_results.size() * 6);
    for (unsigned int i = 0; i < intersection_results.size(); i++) {
        auto &t = mesh.tets[intersection_results[i]];
        for (int j = 0; j < 3; j++) {
            std::array<int, 2> e = {{t[0], t[j + 1]}};
            if (e[0] > e[1])
                std::swap(e[0], e[1]);
            edges.push_back(e);
            e = {{t[j + 1], t[mod3(j + 1) + 1]}};
            if (e[0] > e[1])
                std::swap(e[0], e[1]);
            edges.push_back(e);
        }
    }
    vector_unique(edges);

    std::unordered_set<int> set_intersection_results(intersection_results.begin(), intersection_results.end());
    int M = 2;
    for (int m = 0; m < M; m++) {
        std::vector<bool> is_crossed(edges.size(), false);
        for (int i = 0; i < edges.size(); i++)
            is_crossed[i] = is_crossing(oris[edges[i][0]], oris[edges[i][1]]);
        int old_edges_size = edges.size();
        int cnt_insert = 0;
        for (int i = 0; i < old_edges_size; i++) {
            if (!is_crossed[i])
                continue;
            for (int j = 0; j < 2; j++) {
                const int v_id = edges[i][j];
                assert(oris[v_id] != Predicates::ORI_UNKNOWN);
                if (oris[v_id] != Predicates::ORI_ZERO && abs(signed_dist[v_id]) < mesh.params.eps_coplanar) {
                    if (m == M - 1) {
                        oris[v_id] = Predicates::ORI_ZERO;
                        continue;
                    }

                    for (int t_id: tet_vertices[v_id].conn_tets) {
                        //insert new tets/edges
//                        if (std::find(intersection_results.begin(), intersection_results.end(), t_id) ==
//                            intersection_results.end()) {
//                            intersection_results.push_back(t_id);
                        if (set_intersection_results.find(t_id) != set_intersection_results.end())
                            continue;

                        //you CANNOT use bbox here!
//                        Vector3 min_t, max_t;
//                        get_bbox_tet(tet_vertices[tets[t_id][0]].pos, tet_vertices[tets[t_id][1]].pos,
//                                     tet_vertices[tets[t_id][2]].pos, tet_vertices[tets[t_id][3]].pos,
//                                     min_t, max_t);
//                        if (is_bbox_intersected(min_f, max_f, min_t, max_t))
//                            continue;

                        //update oris/signed_dists
                        for (int k = 0; k < 4; k++) {
                            if (oris[tets[t_id][k]] != Predicates::ORI_UNKNOWN)
                                continue;
                            oris[tets[t_id][k]] = Predicates::orient_3d(a, b, c, tet_vertices[tets[t_id][k]].pos);
                            if (oris[tets[t_id][k]] == Predicates::ORI_ZERO) {
                                signed_dist[tets[t_id][k]] = 0;
                                continue;
                            }
                            signed_dist[tets[t_id][k]] = n.dot(tet_vertices[tets[t_id][k]].pos - a);
                        }

                        for (int k = 0; k < 3; k++) {
                            std::array<int, 2> ee = {{tets[t_id][0], tets[t_id][k + 1]}};
                            if (ee[0] > ee[1])
                                std::swap(ee[0], ee[1]);
                            edges.push_back(ee);
                            ee = {{tets[t_id][k + 1], tets[t_id][mod3(k + 1) + 1]}};
                            if (ee[0] > ee[1])
                                std::swap(ee[0], ee[1]);
                            edges.push_back(ee);
                        }

                        set_intersection_results.insert(t_id);
                        cnt_insert++;
                    }
                }
            }
        }
//        vector_unique(intersection_results);
        vector_unique(edges);
    }
//    vector_unique(intersection_results);
    intersection_results.clear();
    if(!set_intersection_results.empty()) {
        intersection_results = std::vector<int>(set_intersection_results.begin(), set_intersection_results.end());
        std::sort(intersection_results.begin(), intersection_results.end());//has to be sorted
    } else {
        cout<<"intersection_results empty!!!"<<endl;
        pausee();
    }

//    ///////test
//    for (int t_id:intersection_results) {
//        for (int j = 0; j < 4; j++) {
//            if (oris[tets[t_id][j]] == Predicates::ORI_UNKNOWN) {
//                cout << "oris[tets[t_id][j]] == Predicates::ORI_UNKNOWN" << endl;
//                pausee();
//            }
//        }
//    }
////    cout << "edges.size = " << edges.size() << endl;
//    for (int i = 0; i < edges.size(); i++) {
//        auto &e = edges[i];
//        if (is_crossing(oris[e[0]], oris[e[1]])) {
//            if (abs(signed_dist[e[0]]) < mesh.params.eps_coplanar ||
//                abs(signed_dist[e[1]]) < mesh.params.eps_coplanar) {
//                cout << "cross edge un-snapped" << endl;
//                pausee();
//            }
//        }
//    }
//    ///////test
}

int floatTetWild::one_face_cut(const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces, int f_id,
        Mesh& mesh, std::vector<int>& intersection_results,
        std::vector<std::array<std::vector<int>, 4>>& cut_f_ids, std::vector<int>& oris, bool is_parallel, bool is_again) {
    igl::Timer timer;
//    timer.start();

    const Vector3 &a = input_vertices[input_faces[f_id][0]];
    const Vector3 &b = input_vertices[input_faces[f_id][1]];
    const Vector3 &c = input_vertices[input_faces[f_id][2]];
    Vector3 n = get_normal(a, b, c);
    auto &tet_vertices = mesh.tet_vertices;
    auto &tets = mesh.tets;

    if(is_parallel && is_cutting_cross_partitions(mesh, intersection_results))
        return FAIL;

    timer.start();
    std::vector<Scalar> signed_dist;
    snapping(input_vertices, input_faces, f_id, mesh, intersection_results, oris, signed_dist);
    one_face_cut_other_time += timer.getElapsedTimeInSec();

    if(is_parallel && is_cutting_cross_partitions(mesh, intersection_results))
        return FAIL;

//    if(intersection_results.size()>50)
//        return FAIL;

    std::vector<std::array<int, 2>> edges;
    edges.reserve(intersection_results.size() * 6);
    for (unsigned int i = 0; i < intersection_results.size(); i++) {
        auto &t = mesh.tets[intersection_results[i]];
        for (int j = 0; j < 3; j++) {
            std::array<int, 2> e = {{t[0], t[j + 1]}};
            if (e[0] > e[1])
                std::swap(e[0], e[1]);
            edges.push_back(e);
            e = {{t[j + 1], t[mod3(j + 1) + 1]}};
            if (e[0] > e[1])
                std::swap(e[0], e[1]);
            edges.push_back(e);
        }
    }
    vector_unique(edges);

    std::map<std::array<int, 2>, int> e_p_map;
    std::vector<MeshVertex> tmp_tet_vertices;
    tmp_tet_vertices.reserve(edges.size());
    for (int i = 0; i < edges.size(); i++) {
        auto &e = edges[i];
        if (!is_crossing(oris[e[0]], oris[e[1]]))
            continue;

        Vector3 p;
        Scalar d1;
        if (!seg_plane_intersection(tet_vertices[e[0]].pos, tet_vertices[e[1]].pos, a, n, p, d1)) {
            cout << "!seg_plane_intersection !!!" << endl;
            cout<<signed_dist[e[0]]<<endl;
            cout<<signed_dist[e[1]]<<endl;
            cout<<oris[e[0]]<<endl;
            cout<<oris[e[1]]<<endl;
            cout<<e[0]<<" "<<e[1]<<endl;
            cout<<mesh.params.eps_coplanar<<endl;
            pausee();
        } else {
            tmp_tet_vertices.emplace_back(p);
            e_p_map[e] = tet_vertices.size() + tmp_tet_vertices.size();
        }
    }

    std::vector<int> intersection_results_wn;
    intersection_results_wn.reserve(intersection_results.size() * 4 * 6);
    for (int t_id: intersection_results) {
        for (int j = 0; j < 4; j++) {
            intersection_results_wn.insert(intersection_results_wn.end(),
                                           tet_vertices[tets[t_id][j]].conn_tets.begin(),
                                           tet_vertices[tets[t_id][j]].conn_tets.end());
        }
    }
    vector_unique(intersection_results_wn);

//    one_face_cut_other_time += timer.getElapsedTimeInSec();

    if(is_parallel) {
        std::vector<int> tmp_wn;
        tmp_wn.reserve(intersection_results_wn.size() - intersection_results.size());
        std::set_difference(intersection_results_wn.begin(), intersection_results_wn.end(),
                            intersection_results.begin(), intersection_results.end(),
                            std::back_inserter(tmp_wn));
        if (is_cutting_cross_partitions(mesh, tmp_wn))
            return FAIL;
    }

    int result = tet_subdivision(f_id, mesh, intersection_results, intersection_results_wn, cut_f_ids, signed_dist, oris, e_p_map,
                           tmp_tet_vertices, is_again);

    ///////test
    if(false && result == 1){
        cout << "find_bad_cover_for_trif " << f_id << endl;
        find_bad_cover_for_trif(COMMON_INPUT_FOR_CHECK, f_id);
    }
    ///////test

    return result;
}

int floatTetWild::tet_subdivision(int f_id,
        Mesh& mesh, std::vector<int>& intersection_results, std::vector<int>& intersection_results_wn,
        std::vector<std::array<std::vector<int>, 4>>& cut_f_ids, std::vector<Scalar>& signed_dist,std::vector<int>& oris,
        std::map<std::array<int, 2>, int>& e_p_map, std::vector<MeshVertex>& tmp_tet_vertices,
        bool is_again){
    igl::Timer timer;
    igl::Timer timer1;
    timer.start();
    auto &tet_vertices = mesh.tet_vertices;
    auto &tets = mesh.tets;

    ////categorize ALL tets
    std::vector<MeshTet> tmp_tets;
    std::vector<std::array<std::vector<int>, 4>> tmp_cut_f_ids;
    static const std::array<std::array<int, 2>, 6> t_es = {{{{0, 1}}, {{1, 2}}, {{2, 0}}, {{0, 3}}, {{1, 3}}, {{2, 3}}}};
    static const std::array<std::array<int, 3>, 4> t_f_es = {{{{1, 5, 4}}, {{5, 3, 2}}, {{3, 0, 4}}, {{0, 1, 2}}}};
    static const std::array<std::array<int, 3>, 4> t_f_vs = {{{{3, 1, 2}}, {{0, 2, 3}}, {{1, 3, 0}}, {{2, 0, 1}}}};
//    ///clean up invalid cases
//    for (int t_id: intersection_results_wn) {
//        int cnt_zero = 0;
//        int j_non_zero;
//        for (int j = 0; j < 4; j++) {
//            if (oris[tets[t_id][j]] == Predicates::ORI_ZERO)
//                cnt_zero++;
//            else
//                j_non_zero = j;
//        }
//        if (cnt_zero == 4) {
//            Scalar max_dist = 0;
//            for (int j = 0; j < 4; j++) {
//                if (abs(signed_dist[tets[t_id][j]]) > max_dist) {
//                    max_dist = abs(signed_dist[tets[t_id][j]]);
//                    j_non_zero = j;
//                }
//            }
//            if (signed_dist[tets[t_id][j_non_zero]] < 0)
//                oris[tets[t_id][j_non_zero]] = Predicates::ORI_NEGATIVE;
//            else
//                oris[tets[t_id][j_non_zero]] = Predicates::ORI_POSITIVE;
//        }
//    }

//    if(f_id == -1){
//        for(int i=0;i<intersection_results.size();i++){
//            int t_id = intersection_results[i];
//            cout<<t_id<<": "<<i<<" "<<(1+i*4+1)<<endl;
//            for(int j=0;j<4;j++)
//                cout<<tets[t_id][j]<<" ";
//            cout<<endl;
//            for(int j=0;j<4;j++)
//                cout<<tets[t_id].opp_t_ids[j]<<" ";
//            cout<<endl;
//            for(int j=0;j<4;j++)
//                cout<<oris[tets[t_id][j]]<<" ";
//            cout<<endl;
//            for(int j=0;j<4;j++)
//                cout<<signed_dist[tets[t_id][j]]<<" ";
//            cout<<endl;
//        }
////        for (int t_id: intersection_results_wn)
////            cout<<t_id<<" ";
////        cout<<endl;
//        pausee();
//    }

    int i_intersected = 0;
    std::vector<bool> is_updated(tets.size(), false);
    std::vector<int> updated_t_ids;
    updated_t_ids.reserve(intersection_results_wn.size());
    std::vector<int> is_coplanar(tets.size(), -1);
    //reuse containers
    std::vector<Vector2i> my_diags;
    std::vector<int> global_v_ids;
    for (int t_id: intersection_results_wn) {
//        if(tets[t_id].is_removed)
//            cout<<"tets[t_id].is_removed"<<endl;

        bool is_neighbor = true;
        if(!intersection_results.empty()
            && i_intersected < intersection_results.size()
            && t_id == intersection_results[i_intersected]){
            i_intersected++;
            is_neighbor = false;
        }

        //check if coplanar faces
        int cnt_zero = 0;
        int j_non_zero;
        for (int j = 0; j < 4; j++) {
            if (oris[tets[t_id][j]] == Predicates::ORI_ZERO)
                cnt_zero++;
            else
                j_non_zero = j;
        }
        if (cnt_zero == 3) {
            is_coplanar[t_id] = j_non_zero;
            continue;
        } else if (cnt_zero == 4) {
            is_coplanar[t_id] = 4;
            continue;
        }

        //get idx
        std::bitset<6> idx;
        std::array<int, 6> global_p_ids;
        for (int i = 0; i < t_es.size(); i++) {
            std::array<int, 2> e = {{tets[t_id][t_es[i][0]], tets[t_id][t_es[i][1]]}};
            if (e[0] > e[1])
                std::swap(e[0], e[1]);
            const auto it = e_p_map.find(e);
            global_p_ids[i] = (it == e_p_map.end()) ? -1 : (it->second - 1); //[e] - 1;
            if (global_p_ids[i] >= 0)
                idx.set(i);
        }
        if (idx.none()) //not cutted
            continue;

        //match diags
        my_diags.clear();
        for (int i = 0; i < t_f_es.size(); i++) {
            std::array<int, 3> g_p_id = {{-1, -1, -1}};
            int cnt = 0;
            for (int j = 0; j < 3; j++) {
                if (idx[t_f_es[i][j]]) {
                    g_p_id[j] = global_p_ids[t_f_es[i][j]];
                    cnt++;
                }
            }
            if (cnt < 2)
                continue;

            auto it = std::max_element(g_p_id.begin(), g_p_id.end());
            cnt = 0;
            for (int ii:global_p_ids) {
                if (ii == *it)
                    break;
                if (ii >= 0)
                    cnt++;
            }
            int v_id = cnt + 4;
            int opp_v_id = t_f_vs[i][it - g_p_id.begin()];
            if (v_id < opp_v_id)
                my_diags.emplace_back(v_id, opp_v_id);
            else
                my_diags.emplace_back(opp_v_id, v_id);
        }
        std::sort(my_diags.begin(), my_diags.end(), [](const Vector2i& a, const Vector2i& b) {
            return std::tuple<int, int>(a[0], a[1]) < std::tuple<int, int>(b[0], b[1]);
        });

        int idx_int = idx.to_ulong();
        const auto &all_diags = CutTable::get_diag_confs(idx_int);
        if(all_diags.empty()) {
//            cout<<"invalid conf"<<endl;
//
//            cout<<"tet "<<t_id<<endl;
//            cout<<oris[tets[t_id][0]]<<" "<<oris[tets[t_id][1]]<<" "<<oris[tets[t_id][2]]<<" "<<oris[tets[t_id][3]]<<endl;
//            cout<<tets[t_id][0]<<" "<<tets[t_id][1]<<" "<<tets[t_id][2]<<" "<<tets[t_id][3]<<endl;
//            cout<<idx_int<<endl;
//            for(auto& ii:global_p_ids)
//                cout<<ii<<" ";
//            cout<<endl;
            tet_subdivision_time += timer.getElapsedTimeInSec();
            return FAIL_CONFIG;
        }

        //select confs
        global_v_ids.clear();
        global_v_ids.reserve(4+global_p_ids.size());
        for(int j=0;j<4;j++)
            global_v_ids.push_back(tets[t_id][j]);
        for(int i=0;i<global_p_ids.size();i++) {
            if (global_p_ids[i] >= 0)
                global_v_ids.push_back(global_p_ids[i]);
        }

        int config_id = 0;
        if(my_diags.empty()) {
            const auto &new_tets = CutTable::get_tet_conf(idx_int, config_id);
            std::vector<Vector3> centroids;
            get_centroids(mesh, new_tets, global_v_ids, tmp_tet_vertices, centroids);
            Scalar min_volume = get_min_volume(mesh, new_tets, centroids, global_v_ids, tmp_tet_vertices);
            if (min_volume <= SCALAR_ZERO_3) {
//                cout << "inverted 1 " << min_volume << endl;
                timer.stop();
                tet_subdivision_time += timer.getElapsedTimeInSec();
                return FAIL_INVERSION;
            } else {
                for (auto &c:centroids) {
                    tmp_tet_vertices.push_back(MeshVertex(c));
                    global_v_ids.push_back(tmp_tet_vertices.size() + tet_vertices.size() - 1);
                }
            }
        } else {
            Scalar max_min_volume = -1;
            std::vector<Vector3> centroids;

//            int i = std::find(all_diags.begin(), all_diags.end(), my_diags) - all_diags.begin();
//            for (; i < all_diags.size(); i++) {
//                if (my_diags != all_diags[i])
//                    break;
//                const auto &new_tets = CutTable::get_tet_conf(idx_int, i);
//                std::vector<Vector3> tmp_centroids;
//                get_centroids(mesh, new_tets, global_v_ids, tmp_tet_vertices, tmp_centroids);
//                Scalar min_volume = get_min_volume(mesh, new_tets, tmp_centroids, global_v_ids, tmp_tet_vertices);
//                if (min_volume > max_min_volume) {
//                    max_min_volume = min_volume;
//                    centroids = tmp_centroids;
//                    config_id = i;
//                }
//            }

            for (int i = 0; i < all_diags.size(); i++) {
                if (my_diags == all_diags[i]) {
                    const auto &new_tets = CutTable::get_tet_conf(idx_int, i);
                    std::vector<Vector3> tmp_centroids;
                    get_centroids(mesh, new_tets, global_v_ids, tmp_tet_vertices, tmp_centroids);
                    Scalar min_volume = get_min_volume(mesh, new_tets, tmp_centroids, global_v_ids, tmp_tet_vertices);
                    if(min_volume > max_min_volume) {
                        max_min_volume = min_volume;
                        centroids = tmp_centroids;
                        config_id = i;
                    }
                }
            }

            if (max_min_volume <= SCALAR_ZERO_3) {
//                cout << "inverted 2" << endl;
//                cout << "idx_int = " << idx_int << endl;
//                cout << "my_diags = " << endl;
//                for (auto &ii:my_diags)
//                    cout << ii[0] << " " << ii[1] << endl;
                tet_subdivision_time += timer.getElapsedTimeInSec();
                return FAIL_INVERSION;
            } else {
                for (auto &c:centroids) {
                    tmp_tet_vertices.push_back(MeshVertex(c));
                    global_v_ids.push_back(tmp_tet_vertices.size() + tet_vertices.size() - 1);
                }
            }
        }

        const auto &new_tets = CutTable::get_tet_conf(idx_int, config_id);
        const auto &new_is_surface_fs = CutTable::get_surface_conf(idx_int, config_id);
        const auto &new_local_f_ids = CutTable::get_face_id_conf(idx_int, config_id);
        //add new_tets
        for (int i = 0; i < new_tets.size(); i++) {
            tmp_tets.emplace_back(global_v_ids[new_tets[i][0]], global_v_ids[new_tets[i][1]],
                                  global_v_ids[new_tets[i][2]], global_v_ids[new_tets[i][3]]);
            tmp_tets.back().scalar = tets[t_id].scalar;///for tracking partition
            tmp_cut_f_ids.emplace_back();
            for (int j = 0; j < 4; j++) {
                if (new_is_surface_fs[i][j] && !is_neighbor) {
                    (tmp_cut_f_ids.back())[j].push_back(f_id);
//                    if(f_id == 38 && (tmp_tets.back())[mod4(j + 1)] == 550 && (tmp_tets.back())[mod4(j + 2)] == 544 && (tmp_tets.back())[mod4(j + 3)] == 523) {
//                        auto& v0 = (tmp_tets.back())[mod4(j + 1)]>=tet_vertices.size()?
//                                tmp_tet_vertices[(tmp_tets.back())[mod4(j + 1)] - tet_vertices.size()]:tet_vertices[(tmp_tets.back())[mod4(j + 1)]];
//                        auto& v1 = (tmp_tets.back())[mod4(j + 2)]>=tet_vertices.size()?
//                                   tmp_tet_vertices[(tmp_tets.back())[mod4(j + 2)] - tet_vertices.size()]:tet_vertices[(tmp_tets.back())[mod4(j + 2)]];
//                        auto& v2 = (tmp_tets.back())[mod4(j + 3)]>=tet_vertices.size()?
//                                   tmp_tet_vertices[(tmp_tets.back())[mod4(j + 3)] - tet_vertices.size()]:tet_vertices[(tmp_tets.back())[mod4(j + 3)]];
//                        Vector3 n = ((v0.pos-v2.pos).cross(v1.pos-v2.pos)).normalized();
//                        if(abs(cutting_n.dot(n))<0.1) {
//                            cout << "abs(cutting_n.dot(n))<0.1" << endl;
//                            cout << abs(cutting_n.dot(n)) << endl;
//                            cout << (tmp_tets.back())[mod4(j + 1)] << " " << (tmp_tets.back())[mod4(j + 2)] << " "
//                                 << (tmp_tets.back())[mod4(j + 3)] << endl;
//                            cout << v0.pos.transpose() << endl;
//                            cout << v1.pos.transpose() << endl;
//                            cout << v2.pos.transpose() << endl;
//                            cout<<idx_int<<endl;
//                            cout<<config_id<<endl;
//                            cout<<new_tets[i][0]<<" "<<new_tets[i][1]<<" "<<new_tets[i][2]<<" "<<new_tets[i][3]<<endl;
//                            cout<<new_is_surface_fs[i][0]<<" "<<new_is_surface_fs[i][1]<<" "<<new_is_surface_fs[i][2]<<" "<<new_is_surface_fs[i][3]<<endl;
//                            cout<<"t_id = "<<t_id<<endl;
//                            cout<<signed_dist[tets[t_id][0]]<<" "<<signed_dist[tets[t_id][1]]<<" "<<signed_dist[tets[t_id][2]]<<" "<<signed_dist[tets[t_id][3]]<<endl;
//                            cout<<oris[tets[t_id][0]]<<" "<<oris[tets[t_id][1]]<<" "<<oris[tets[t_id][2]]<<" "<<oris[tets[t_id][3]]<<endl;
//                            cout<<"is_neighbor "<<is_neighbor<<endl;
//                            cout<<"is found "<<(std::find(intersection_results.begin(), intersection_results.end(), t_id)!=intersection_results.end())<<endl;
//                            pausee();
//                        }
//                    }
                }

                int old_local_f_id = new_local_f_ids[i][j];
                if (old_local_f_id < 0)
                    continue;
                (tmp_cut_f_ids.back())[j].insert((tmp_cut_f_ids.back())[j].end(),
                                                 cut_f_ids[t_id][old_local_f_id].begin(),
                                                 cut_f_ids[t_id][old_local_f_id].end());
                (tmp_tets.back()).is_bbox_fs[j] = tets[t_id].is_bbox_fs[old_local_f_id];
                (tmp_tets.back()).is_surface_fs[j] = tets[t_id].is_surface_fs[old_local_f_id];
                (tmp_tets.back()).surface_tags[j] = tets[t_id].surface_tags[old_local_f_id];
            }
        }

        is_updated[t_id] = true;
        updated_t_ids.push_back(t_id);
    }

    for(int i=0;i<is_coplanar.size();i++) {
        if (is_coplanar[i] >= 0) {
            if (is_coplanar[i] < 4) {
                cut_f_ids[i][is_coplanar[i]].push_back(f_id);
            } else {
                for (int k = 0; k < 4; k++)
                    cut_f_ids[i][k].push_back(f_id);
            }
        }
    }

    if(tmp_tet_vertices.empty()){
        tet_subdivision_time += timer.getElapsedTimeInSec();
        return SUCCESS;
    }
    tet_subdivision_time += timer.getElapsedTimeInSec();

    timer.start();
    //push back new elements
    tet_vertices.insert(tet_vertices.end(), tmp_tet_vertices.begin(), tmp_tet_vertices.end());
    std::vector<bool> is_dirty(tet_vertices.size(), false);
    std::vector<int> dirty_t_ids;
    dirty_t_ids.reserve(tmp_tets.size());
    for (auto &t: tmp_tets) {
        for (int j = 0; j < 4; j++) {
            is_dirty[t[j]] = true;
            dirty_t_ids.insert(dirty_t_ids.end(),
                               tet_vertices[t[j]].conn_tets.begin(), tet_vertices[t[j]].conn_tets.end());
            tet_vertices[t[j]].conn_tets.clear();
        }
        if(is_again)
            t.quality = get_quality(mesh, t);
    }

    int cnt = 0;
//    for (int i = 0; i < tets.size(); i++) {
//        if (is_updated[i]) {
//            tets[i] = tmp_tets[cnt];
//            cut_f_ids[i] = tmp_cut_f_ids[cnt];
//            cnt++;
//        }
//    }

    for (int t_id: updated_t_ids) {
        tets[t_id] = tmp_tets[cnt];
        cut_f_ids[t_id] = tmp_cut_f_ids[cnt];
        dirty_t_ids.push_back(t_id);
        cnt++;
    }
    if(is_again) {
        int i = mesh.t_empty_start;
        for (; i < tets.size(); i++) {
            if (tets[i].is_removed) {
                tets[i] = tmp_tets[cnt];
                cut_f_ids[i] = tmp_cut_f_ids[cnt];
                dirty_t_ids.push_back(i);
                cnt++;

                if (cnt >= tmp_tets.size())
                    break;
            }
        }
        mesh.t_empty_start = i + 1;
    }
    if(cnt < tmp_tets.size()) {
        const int tets_size = tets.size();
        tets.insert(tets.end(), tmp_tets.begin() + cnt, tmp_tets.end());
        cut_f_ids.insert(cut_f_ids.end(), tmp_cut_f_ids.begin() + cnt, tmp_cut_f_ids.end());

        const int old_dirty_t_ids_size = dirty_t_ids.size();
        dirty_t_ids.resize(old_dirty_t_ids_size + tmp_tets.size() - cnt);
        std::iota(dirty_t_ids.begin() + old_dirty_t_ids_size, dirty_t_ids.end(), tets_size);
    }
    vector_unique(dirty_t_ids);

    //update conn_tets
    for(int t_id: dirty_t_ids){
        bool is_reset = false;
        for (int j = 0; j < 4; j++) {
            if (is_dirty[tets[t_id][j]]) {
//                tet_vertices[tets[t_id][j]].conn_tets.insert(t_id);
                tet_vertices[tets[t_id][j]].conn_tets.push_back(t_id);
                if (!is_reset) {
                    tets[t_id].opp_t_ids = {{OPP_T_ID_UNKNOWN, OPP_T_ID_UNKNOWN, OPP_T_ID_UNKNOWN, OPP_T_ID_UNKNOWN}};
                    is_reset = true;
                }
            }
        }
    }

//    for (int i = 0; i < tets.size(); i++) {
//        if(tets[i].is_removed)
//            continue;
//
//        bool is_reset = false;
//        for (int j = 0; j < 4; j++) {
//            if (is_dirty[tets[i][j]]) {
//                tet_vertices[tets[i][j]].conn_tets.insert(i);
//                if (!is_reset) {
//                    tets[i].opp_t_ids = {{OPP_T_ID_UNKNOWN, OPP_T_ID_UNKNOWN, OPP_T_ID_UNKNOWN, OPP_T_ID_UNKNOWN}};
//                    is_reset = true;
//                }
//            }
//        }
//    }

    //update opp_t_ids
//    for (int i = 0; i < tets.size(); i++) {
//        if (tets[i].is_removed)
//            continue;
//
//        if (is_dirty[tets[i][0]] || is_dirty[tets[i][1]] || is_dirty[tets[i][2]] || is_dirty[tets[i][3]]) {
//            for (int j = 0; j < 4; j++) {
//                if (tets[i].opp_t_ids[j] >= 0)
//                    continue;
//                std::vector<int> pair;
//                set_intersection(mesh.tet_vertices[tets[i][(j + 1) % 4]].conn_tets,
//                                 mesh.tet_vertices[tets[i][(j + 2) % 4]].conn_tets,
//                                 mesh.tet_vertices[tets[i][(j + 3) % 4]].conn_tets, pair);
//                if (pair.size() == 2) {
//                    int opp_t_id = pair[0] == i ? pair[1] : pair[0];
//                    tets[i].opp_t_ids[j] = opp_t_id;
//                    for (int k = 0; k < 4; k++) {
//                        if (tets[opp_t_id][k] != tets[i][mod4(j + 1)]
//                            && tets[opp_t_id][k] != tets[i][mod4(j + 2)]
//                            && tets[opp_t_id][k] != tets[i][mod4(j + 3)]) {
//                            tets[opp_t_id].opp_t_ids[k] = i;
//                            break;
//                        }
//                    }
//                }
//            }
//        }
//    }

    update_opp_t_ids_time += timer.getElapsedTimeInSec();

    return SUCCESS;
}

floatTetWild::Scalar floatTetWild::get_min_volume(const Mesh &mesh, const std::vector<Vector4i>& new_tets, const std::vector<Vector3>& centroids,
        const std::vector<int>& global_v_ids, const std::vector<MeshVertex>& tmp_tet_vertices){
    auto &tet_vertices = mesh.tet_vertices;
    auto &tets = mesh.tets;

    Scalar min = INT_MAX;
    for (const auto &t: new_tets) {
        auto &v0 = t[0] < global_v_ids.size() ? (global_v_ids[t[0]] < tet_vertices.size() ? tet_vertices[global_v_ids[t[0]]].pos : tmp_tet_vertices[global_v_ids[t[0]] - tet_vertices.size()].pos)
                                              : centroids[t[0] - global_v_ids.size()];
        auto &v1 = t[1] < global_v_ids.size() ? (global_v_ids[t[1]] < tet_vertices.size() ? tet_vertices[global_v_ids[t[1]]].pos : tmp_tet_vertices[global_v_ids[t[1]] - tet_vertices.size()].pos)
                                              : centroids[t[1] - global_v_ids.size()];
        auto &v2 = t[2] < global_v_ids.size() ? (global_v_ids[t[2]] < tet_vertices.size() ? tet_vertices[global_v_ids[t[2]]].pos : tmp_tet_vertices[global_v_ids[t[2]] - tet_vertices.size()].pos)
                                              : centroids[t[2] - global_v_ids.size()];
        auto &v3 = t[3] < global_v_ids.size() ? (global_v_ids[t[3]] < tet_vertices.size() ? tet_vertices[global_v_ids[t[3]]].pos : tmp_tet_vertices[global_v_ids[t[3]] - tet_vertices.size()].pos)
                                              : centroids[t[3] - global_v_ids.size()];

        Scalar volume = Predicates::orient_3d_volume(v0, v1, v2, v3);
        if (volume < min)
            min = volume;
    }

    return min;
}

void floatTetWild::get_centroids(const Mesh &mesh, const std::vector<Vector4i>& new_tets, const std::vector<int>& global_v_ids, const std::vector<MeshVertex>& tmp_tet_vertices,  std::vector<Vector3>& cs) {
    const auto &tet_vertices = mesh.tet_vertices;
    const auto &tets = mesh.tets;

    std::array<std::vector<int>, 2> centroid_n_v_ids;
    for (int i = 0; i < new_tets.size(); i++) {
        for (int j = 0; j < 4; j++) {
            if (new_tets[i][j] >= global_v_ids.size()) {
                int slot = new_tets[i][j] - global_v_ids.size();
                centroid_n_v_ids[slot].push_back(global_v_ids[new_tets[i][mod4(j + 1)]]);
                centroid_n_v_ids[slot].push_back(global_v_ids[new_tets[i][mod4(j + 2)]]);
                centroid_n_v_ids[slot].push_back(global_v_ids[new_tets[i][mod4(j + 3)]]);
                break;
            }
        }
    }
    for (auto &c_n_v_ids: centroid_n_v_ids) {
        if (c_n_v_ids.empty())
            break;
        vector_unique(c_n_v_ids);
        Vector3 p(0, 0, 0);
        for (int v_id: c_n_v_ids) {
            auto &tmp_p = v_id < tet_vertices.size() ? tet_vertices[v_id].pos
            : tmp_tet_vertices[v_id - tet_vertices.size()].pos;
            p += tmp_p;
        }
        p /= c_n_v_ids.size();
        cs.push_back(p);
    }
}

void floatTetWild::track_surface(Mesh &mesh,
                                 const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces, const std::vector<int> &input_tags,
                                 std::vector<std::array<std::vector<int>, 4>> &cut_f_ids,
                                 const AABBWrapper& tree, const std::vector<bool>& is_boundary_preserved) {
    auto &tet_vertices = mesh.tet_vertices;
    auto &tets = mesh.tets;

    for (int t_id = 0; t_id < cut_f_ids.size(); t_id++) {
        if (tets[t_id].is_removed)
            continue;

        for (int j = 0; j < 4; j++) {
            if (cut_f_ids[t_id][j].empty())
                continue;
            if (tets[t_id].is_surface_fs[j] != NOT_SURFACE)
                continue;

//            ///////test
//            int f_id = cut_f_ids[t_id][j].front();
//            std::array<int, 3> f_tet = {{tets[t_id][mod4(j + 1)], tets[t_id][mod4(j + 2)], tets[t_id][mod4(j + 3)]}};
//            tets[t_id].is_surface_fs[j] = Predicates::orient_3d(
//                    input_vertices[input_faces[f_id][0]], input_vertices[input_faces[f_id][1]],
//                    input_vertices[input_faces[f_id][2]], mesh.tet_vertices[tets[t_id][j]].pos);
//            int opp_t_id = tets[t_id].opp_t_ids[j];
//            for (int k = 0; k < 4; k++) {
//                if (tets[opp_t_id][k] != f_tet[0] && tets[opp_t_id][k] != f_tet[1]
//                    && tets[opp_t_id][k] != f_tet[2]) {
//                    tets[opp_t_id].is_surface_fs[k] = -tets[t_id].is_surface_fs[j];
//                    break;
//                }
//            }
//            continue;
//            ///////test

            Vector3 &v0 = tet_vertices[tets[t_id][mod4(j + 1)]].pos;
            Vector3 &v1 = tet_vertices[tets[t_id][mod4(j + 2)]].pos;
            Vector3 &v2 = tet_vertices[tets[t_id][mod4(j + 3)]].pos;
            Vector3 n_tf = (v1 - v0).cross(v2 - v0);
//            Vector3 n_tf = get_normal(v0, v1, v2);

            bool is_found = false;

            int t = get_t(v0, v1, v2);
            std::array<Vector2, 3> vs_tet = {{to_2d(v0, t), to_2d(v1, t), to_2d(v2, t)}};
            for (int f_id:cut_f_ids[t_id][j]) {
                if (!is_boundary_preserved[f_id])
                    continue;

                std::array<int, 3> f_in = {{input_faces[f_id][0], input_faces[f_id][1], input_faces[f_id][2]}};
                std::sort(f_in.begin(), f_in.end());
                std::array<int, 3> f_tet = {{tets[t_id][mod4(j + 1)], tets[t_id][mod4(j + 2)],
                                                    tets[t_id][mod4(j + 3)]}};
                std::sort(f_tet.begin(), f_tet.end());
                if (f_in == f_tet) {
                    tets[t_id].is_surface_fs[j] = Predicates::orient_3d(
                            input_vertices[input_faces[f_id][0]], input_vertices[input_faces[f_id][1]],
                            input_vertices[input_faces[f_id][2]], mesh.tet_vertices[tets[t_id][j]].pos);
                    tets[t_id].surface_tags[j] = input_tags[f_id];
                    int opp_t_id = tets[t_id].opp_t_ids[j];
//                    for (int k = 0; k < 4; k++) {
//                        if (tets[opp_t_id][k] != f_tet[0] && tets[opp_t_id][k] != f_tet[1]
//                            && tets[opp_t_id][k] != f_tet[2]) {
//                            tets[opp_t_id].is_surface_fs[k] = -tets[t_id].is_surface_fs[j];
//                            break;
//                        }
//                    }

                    if (opp_t_id == OPP_T_ID_UNKNOWN)
                        opp_t_id = get_opp_t_id(mesh, t_id, j);
                    if (opp_t_id != OPP_T_ID_BOUNDARY) {
                        int k = mesh.tets[opp_t_id].find_opp(tets[t_id][mod4(j + 1)], tets[t_id][mod4(j + 2)],
                                                             tets[t_id][mod4(j + 3)]);
                        mesh.tets[opp_t_id].opp_t_ids[k] = t_id;
                        tets[opp_t_id].is_surface_fs[k] = -tets[t_id].is_surface_fs[j];
                        tets[opp_t_id].surface_tags[k] = input_tags[f_id];
                    }

                    is_found = true;
                    break;
                }

                //check
                std::array<Vector2, 3> vs_tri = {{to_2d(input_vertices[input_faces[f_id][0]], t),
                                                         to_2d(input_vertices[input_faces[f_id][1]], t),
                                                         to_2d(input_vertices[input_faces[f_id][2]], t)}};

//                Vector2 c = (vs_tet[0] + vs_tet[1] + vs_tet[2]) / 3;
//                int cnt_pos = 0;
//                int cnt_neg = 0;
//                for (int i = 0; i < 3; i++) {
//                    int ori = Predicates::orient_2d(vs_tri[i], vs_tri[mod3(i + 1)], c);
//                    if (ori == Predicates::ORI_POSITIVE)
//                        cnt_pos++;
//                    else if (ori == Predicates::ORI_NEGATIVE)
//                        cnt_neg++;
//                }
//                if (!(cnt_pos == 0 || cnt_neg == 0)) {
//                    continue;
//                }

                Vector2 c0 = (vs_tet[0] + vs_tet[1] + vs_tet[2]) / 3;
                std::vector<Vector2> cs;
                cs.push_back(c0);
                cs.push_back((vs_tet[0] + c0) / 2);
                cs.push_back((vs_tet[1] + c0) / 2);
                cs.push_back((vs_tet[2] + c0) / 2);

                //add more samples.
                //this bug only happens on planar regions
                //a flat tet goes across two coplanar input triangles
                //but the tet faces are marked as coplanar to different triangles
                //the centroids could be a little off
//                bool find_one_outside = false;
//                for(auto& c: cs) {
//                    int cnt_pos = 0;
//                    int cnt_neg = 0;
//                    for (int i = 0; i < 3; i++) {
//                        int ori = Predicates::orient_2d(vs_tri[i], vs_tri[mod3(i + 1)], c);
//                        if (ori == Predicates::ORI_POSITIVE)
//                            cnt_pos++;
//                        else if (ori == Predicates::ORI_NEGATIVE)
//                            cnt_neg++;
//                    }
//                    if (!(cnt_pos == 0 || cnt_neg == 0)) {
////                        continue;
//                        find_one_outside = true;
//                        break;
//                    }
////                    find_one_inside = true;
////                    break;
//                }
//                if(find_one_outside)
//                    continue;

                bool find_one_inside = false;
                for (auto &c: cs) {
                    int cnt_pos = 0;
                    int cnt_neg = 0;
                    for (int i = 0; i < 3; i++) {
                        int ori = Predicates::orient_2d(vs_tri[i], vs_tri[mod3(i + 1)], c);
                        if (ori == Predicates::ORI_POSITIVE)
                            cnt_pos++;
                        else if (ori == Predicates::ORI_NEGATIVE)
                            cnt_neg++;
                    }
                    if (!(cnt_pos == 0 || cnt_neg == 0)) {
                        continue;
                    }
                    find_one_inside = true;
                    break;
                }
                if (!find_one_inside)
                    continue;

                //double check
                std::vector<GEO::vec3> ps;
                sample_triangle({{v0, v1, v2}}, ps, mesh.params.dd);
                if (tree.is_out_sf_envelope(ps, mesh.params.eps_2))
                    continue;

                //mark
                int opp_t_id = tets[t_id].opp_t_ids[j];
                if (opp_t_id == OPP_T_ID_UNKNOWN)
                    opp_t_id = get_opp_t_id(mesh, t_id, j);
                int k = opp_t_id >= 0 ? mesh.tets[opp_t_id].find_opp(tets[t_id][mod4(j + 1)], tets[t_id][mod4(j + 2)],
                                                                     tets[t_id][mod4(j + 3)]) : -1;

//                int k;
//                for (k = 0; k < 4; k++) {
//                    if (tets[opp_t_id][k] != f_tet[0] && tets[opp_t_id][k] != f_tet[1]
//                        && tets[opp_t_id][k] != f_tet[2])
//                        break;
//                }

                int t_ori = Predicates::orient_3d(
                        input_vertices[input_faces[f_id][0]], input_vertices[input_faces[f_id][1]],
                        input_vertices[input_faces[f_id][2]], mesh.tet_vertices[tets[t_id][j]].pos);
                int opp_t_ori = opp_t_id >= 0 ? Predicates::orient_3d(
                        input_vertices[input_faces[f_id][0]], input_vertices[input_faces[f_id][1]],
                        input_vertices[input_faces[f_id][2]], mesh.tet_vertices[tets[opp_t_id][k]].pos) : -1;
                if (t_ori == Predicates::ORI_POSITIVE && opp_t_ori == Predicates::ORI_NEGATIVE
                    || opp_t_ori == Predicates::ORI_POSITIVE && t_ori == Predicates::ORI_NEGATIVE) {
                    tets[t_id].is_surface_fs[j] = t_ori;
                    if (opp_t_id >= 0)
                        tets[opp_t_id].is_surface_fs[k] = opp_t_ori;
                } else {
                    Vector3 n_f = (input_vertices[input_faces[f_id][1]] - input_vertices[input_faces[f_id][0]]).cross(
                            input_vertices[input_faces[f_id][2]] - input_vertices[input_faces[f_id][0]]);
//                    Vector3 n_f = get_normal(input_vertices[input_faces[f_id][0]], input_vertices[input_faces[f_id][1]],
//                            input_vertices[input_faces[f_id][2]]);
                    Scalar cos_a = n_f.dot(n_tf);
                    int ori;
                    if (cos_a >= 0)
                        ori = Predicates::orient_3d(v0, v1, v2, mesh.tet_vertices[tets[t_id][j]].pos);
                    else
                        ori = Predicates::orient_3d(v0, v2, v1, mesh.tet_vertices[tets[t_id][j]].pos);
                    tets[t_id].is_surface_fs[j] = ori;
                    if (opp_t_id >= 0)
                        tets[opp_t_id].is_surface_fs[k] = -ori;
                }

                tets[t_id].surface_tags[j] = input_tags[f_id];
                if (opp_t_id >= 0)
                    tets[opp_t_id].surface_tags[k] = input_tags[f_id];

                is_found = true;
                break;
            }

//            if(is_found){
//                int t_is_surface_fs = 0;
//
//                //mark
//                int opp_t_id = tets[t_id].opp_t_ids[j];
//
//                if (opp_t_id == OPP_T_ID_UNKNOWN)
//                    opp_t_id = get_opp_t_id(mesh, t_id, j);
//                int k = mesh.tets[opp_t_id].find_opp(tets[t_id][mod4(j + 1)], tets[t_id][mod4(j + 2)],
//                                                     tets[t_id][mod4(j + 3)]);
//
//                for (int f_id:cut_f_ids[t_id][j]) {
////                int k;
////                for (k = 0; k < 4; k++) {
////                    if (tets[opp_t_id][k] != f_tet[0] && tets[opp_t_id][k] != f_tet[1]
////                        && tets[opp_t_id][k] != f_tet[2])
////                        break;
////                }
//
//                    int t_ori = Predicates::orient_3d(
//                            input_vertices[input_faces[f_id][0]], input_vertices[input_faces[f_id][1]],
//                            input_vertices[input_faces[f_id][2]], mesh.tet_vertices[tets[t_id][j]].pos);
//                    int opp_t_ori = Predicates::orient_3d(
//                            input_vertices[input_faces[f_id][0]], input_vertices[input_faces[f_id][1]],
//                            input_vertices[input_faces[f_id][2]], mesh.tet_vertices[tets[opp_t_id][k]].pos);
//                    if (t_ori == Predicates::ORI_POSITIVE && opp_t_ori == Predicates::ORI_NEGATIVE
//                        || opp_t_ori == Predicates::ORI_POSITIVE && t_ori == Predicates::ORI_NEGATIVE) {
//                        t_is_surface_fs += t_ori;
////                        tets[t_id].is_surface_fs[j] = t_ori;
////                        tets[opp_t_id].is_surface_fs[k] = opp_t_ori;
//                    } else {
//                        Vector3 n_f = (input_vertices[input_faces[f_id][1]] -
//                                       input_vertices[input_faces[f_id][0]]).cross(
//                                input_vertices[input_faces[f_id][2]] - input_vertices[input_faces[f_id][0]]);
////                    Vector3 n_f = get_normal(input_vertices[input_faces[f_id][0]], input_vertices[input_faces[f_id][1]],
////                            input_vertices[input_faces[f_id][2]]);
//                        Scalar cos_a = n_f.dot(n_tf);
//                        int ori;
//                        if (cos_a >= 0)
//                            ori = Predicates::orient_3d(v0, v1, v2, mesh.tet_vertices[tets[t_id][j]].pos);
//                        else
//                            ori = Predicates::orient_3d(v0, v2, v1, mesh.tet_vertices[tets[t_id][j]].pos);
//                        t_is_surface_fs += ori;
////                        tets[t_id].is_surface_fs[j] = ori;
////                        tets[opp_t_id].is_surface_fs[k] = -ori;
//                    }
//                }
//                if(t_is_surface_fs > 0){
//                    tets[t_id].is_surface_fs[j] = 1;
//                    tets[opp_t_id].is_surface_fs[k] = -1;
//                } else {
//                    tets[t_id].is_surface_fs[j] = -1;
//                    tets[opp_t_id].is_surface_fs[k] = 1;
//                }
//            }

            ///////test
            continue;
//            if (!is_found)
//                find_bad_cover_for_tetf(COMMON_INPUT_FOR_CHECK, tree, t_id, j, false);
            if (is_found)
                find_bad_cover_for_tetf(COMMON_INPUT_FOR_CHECK, tree, t_id, j, true);
            ///////test
        }
    }

//    std::ofstream fout("bp.obj");
//    for(auto& v:mesh.tet_vertices){
//        if(v.is_removed)
//            continue;
//        if(v.is_on_boundary)
//            fout<<"v "<<v.pos.transpose()<<endl;
//    }
//    fout.close();
}


void floatTetWild::preserve_cutting_open_boundary(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
        const std::vector<int>& current_inserted_cutting_fs,
        std::vector<std::array<std::vector<int>, 4>> &cut_f_ids, Mesh& mesh, AABBWrapper& tree,
        const std::vector<bool>& is_face_matched, std::vector<bool>& is_boundary_preserved, bool is_again){
//    std::ofstream fout("boundayr_p.obj");

    auto &tet_vertices = mesh.tet_vertices;
    auto &tets = mesh.tets;
    std::vector<std::array<int, 2>> edges;
    edges.reserve(current_inserted_cutting_fs.size()*3);
    std::vector<std::vector<int>> conn_tris(input_vertices.size());
    for (int f_id: current_inserted_cutting_fs) {
        for (int j = 0; j < 3; j++) {
            conn_tris[input_faces[f_id][j]].push_back(f_id);
            edges.push_back({{input_faces[f_id][j], input_faces[f_id][mod3(j + 1)]}});
        }
    }
    vector_unique(edges);

    std::vector<std::array<int, 2>> f_id_le_id;
    get_current_open_boundary(input_vertices, input_faces, edges, conn_tris, f_id_le_id);

    std::vector<std::vector<int>> hint_cut_t_ids(input_faces.size());
    for (int t_id = 0; t_id < cut_f_ids.size(); t_id++) {
        if (mesh.tets[t_id].is_removed)
            continue;
        for (int j = 0; j < 4; j++) {
            for (int f_id:cut_f_ids[t_id][j])
                hint_cut_t_ids[f_id].push_back(t_id);
        }
    }

//    if(!mesh.is_input_all_inserted) {
//        std::vector<std::array<int, 2>> b_edges;
//        b_edges.reserve(f_id_le_id.size());
//        for (auto &info: f_id_le_id) {
//            int f_id = info[0];
//            int le_id = info[1];
//            b_edges.push_back({{input_faces[f_id][le_id], input_faces[f_id][mod3(le_id + 1)]}});
//        }
//
//        std::ofstream fout("boundary.obj");
//        for(auto& v:input_vertices){
//            fout<<"v "<<v.transpose()<<endl;
//        }
//        for(auto& e: b_edges){
//            fout<<"l "<<e[0]+1<<" "<<e[1]+1<<endl;
//        }
//        fout.close();
//
//        tree.init_tmp_b_mesh_and_tree(input_vertices, input_faces, b_edges);
//    }


//    cout<<"f_id_le_id.size = "<<f_id_le_id.size()<<endl;
//    cout<<"current_inserted_cutting_fs.size = "<<current_inserted_cutting_fs.size()<<endl;
//    cout<<cnt_1<<endl;
//    cout<<cnt_2<<endl;
    for (auto &info: f_id_le_id) {
        int cur_f_id = info[0];
        int le_id = info[1];

//        if(is_face_matched[cur_f_id]) {
//            continue;
//        }

        int t = get_t(input_vertices[input_faces[cur_f_id][0]], input_vertices[input_faces[cur_f_id][1]],
                      input_vertices[input_faces[cur_f_id][2]]);
        std::array<Vector2, 2> seg = {{to_2d(input_vertices[input_faces[cur_f_id][le_id]], t),
                                              to_2d(input_vertices[input_faces[cur_f_id][mod3(le_id + 1)]], t)}};

        std::vector<MeshVertex> tmp_tet_vertices;
        std::map<std::array<int, 2>, int> e_p_map;
        std::vector<int> intersection_results_wn;

        ////new
        std::vector<bool> is_visited(tets.size(), false);
        std::array<int, 2> seed = {{-1, -1}};
        for (int t_id:hint_cut_t_ids[cur_f_id]) {//hint first
            is_visited[t_id] = true;

            for (int j = 0; j < 4; j++) {
                if (std::find(cut_f_ids[t_id][j].begin(), cut_f_ids[t_id][j].end(), cur_f_id) !=
                    cut_f_ids[t_id][j].end()) {
                    seed[0] = t_id;
                    seed[1] = j;
                    break;
                }
            }
            if (seed[0] >= 0)
                break;
        }
        if (seed[0] < 0) {
            for (int t_id = 0; t_id < cut_f_ids.size(); t_id++) {
                is_visited[t_id] = true;

                if (mesh.tets[t_id].is_removed)
                    continue;

                for (int j = 0; j < 4; j++) {
                    if(cut_f_ids[t_id][j].empty())
                        continue;
                    if (std::find(cut_f_ids[t_id][j].begin(), cut_f_ids[t_id][j].end(), cur_f_id) !=
                        cut_f_ids[t_id][j].end()) {
                        seed[0] = t_id;
                        seed[1] = j;
                        break;
                    }
                }
                if (seed[0] >= 0)
                    break;
            }
        }
        if (seed[0] < 0)
            continue;

        std::vector<std::array<int, 2>> arr_tf_ids;

        std::queue<std::array<int, 2>> tet_queue;
        tet_queue.push(seed);
        while (!tet_queue.empty()) {
            int t_id = tet_queue.front()[0];
            int j = tet_queue.front()[1];
            tet_queue.pop();

            for (int k = 0; k < 3; k++) {
                for (int n_t_id:mesh.tet_vertices[mesh.tets[t_id][mod4(j + 1 + k)]].conn_tets) {
                    if (!is_visited[n_t_id]) {
                        is_visited[n_t_id] = true;
                        for (int kk = 0; kk < 4; kk++) {
                            if (std::find(cut_f_ids[n_t_id][kk].begin(), cut_f_ids[n_t_id][kk].end(), cur_f_id) !=
                                cut_f_ids[n_t_id][kk].end()) {
                                tet_queue.push({{n_t_id, kk}});
                                break;
                            }
                        }
                    }
                }
            }

            //mark close points//todo: maybe just record them in a list and mark them after cutting
            std::array<bool, 3> is_cross = {{true, true, true}};
            for (int i = 0; i < 3; i++) {
                auto &tet_v = mesh.tet_vertices[mesh.tets[t_id][mod4(j + 1 + i)]];
                Scalar dist = p_seg_squared_dist_3d(tet_v.pos,
                                                    input_vertices[input_faces[cur_f_id][le_id]],
                                                    input_vertices[input_faces[cur_f_id][mod3(le_id+1)]]);
                if (dist < mesh.params.eps_2_coplanar) {
//                if (dist < SCALAR_ZERO_2) {
                    is_cross[i] = false;
                    is_cross[mod3(i + 2)] = false;
                    tet_v.is_on_boundary = true;
//                    fout<<"v "<<tet_v.pos.transpose()<<endl;
                }
            }

            if(is_face_matched[cur_f_id]) // no need to do 2d arragement later
                continue;

            //compute insersection points
            std::array<Vector2, 3> f_tet = {{to_2d(mesh.tet_vertices[mesh.tets[t_id][mod4(j + 1)]].pos, t),
                                                    to_2d(mesh.tet_vertices[mesh.tets[t_id][mod4(j + 2)]].pos, t),
                                                    to_2d(mesh.tet_vertices[mesh.tets[t_id][mod4(j + 3)]].pos, t)}};

            std::array<Vector2, 3> ps;
            std::array<Scalar, 3> t2s = {{-1, -1, -1}};
            for (int i = 0; i < 3; i++) {
                std::array<int, 2> e = {{tets[t_id][mod4(j + 1 + i)], tets[t_id][mod4(j + 1 + mod3(i + 1))]}};
                if (e[0] > e[1])
                    std::swap(e[0], e[1]);
//                if (e_p_map[e] - 1 < 0) {
                if (e_p_map.find(e) == e_p_map.end()) {
                    Scalar t2 = -1;
                    if (seg_seg_intersection_2d(seg, {{f_tet[i], f_tet[mod3(i + 1)]}}, t2)) {
                        ps[i] = (1 - t2) * f_tet[i] + t2 * f_tet[mod3(i + 1)];
                        t2s[i] = t2;
                    }
                }
            }
            bool is_cutted = false;
            for (int i = 0; i < 3; i++) {
                if (t2s[i] > 1 || t2s[i] < 0)
                    continue;

                Vector2 v = ps[i] - f_tet[i];
//                if (v.squaredNorm() < SCALAR_ZERO_2) {
                if (v.squaredNorm() < mesh.params.eps_2_coplanar) {
                    mesh.tet_vertices[mesh.tets[t_id][mod4(j + 1 + i)]].is_on_boundary = true;
                    t2s[i] = -1;
                    continue;
                }
                v = ps[i] - f_tet[mod3(i + 1)];
//                if (v.squaredNorm() < SCALAR_ZERO_2) {
                if (v.squaredNorm() < mesh.params.eps_2_coplanar) {
                    mesh.tet_vertices[mesh.tets[t_id][mod4(j + 1 + mod3(i + 1))]].is_on_boundary = true;
                    t2s[i] = -1;
                    continue;
                }

                const Scalar &t2 = t2s[i];
                Vector3 p = (1 - t2) * tet_vertices[tets[t_id][mod4(j + 1 + i)]].pos
                            + t2 * tet_vertices[tets[t_id][mod4(j + 1 + mod3(i + 1))]].pos;
                tmp_tet_vertices.push_back(MeshVertex(p));
                tmp_tet_vertices.back().is_on_boundary = true;//mark intersecting points
//                fout<<"v "<<tmp_tet_vertices.back().pos.transpose()<<endl;
                std::array<int, 2> e = {{tets[t_id][mod4(j + 1 + i)], tets[t_id][mod4(j + 1 + mod3(i + 1))]}};
                if (e[0] > e[1])
                    std::swap(e[0], e[1]);
                e_p_map[e] = tet_vertices.size() + tmp_tet_vertices.size();

                intersection_results_wn.insert(intersection_results_wn.end(), tet_vertices[e[0]].conn_tets.begin(),
                                               tet_vertices[e[0]].conn_tets.end());
                intersection_results_wn.insert(intersection_results_wn.end(), tet_vertices[e[1]].conn_tets.begin(),
                                               tet_vertices[e[1]].conn_tets.end());
                is_cutted = true;
            }
            if(is_cutted)
                arr_tf_ids.push_back({{t_id, j}});
        }
        if (tmp_tet_vertices.empty())
            continue;

        //cutting tets
        int _f_id = 0;
        std::vector<int> _intersection_results;
        std::vector<int> _oris(mesh.tet_vertices.size(), Predicates::ORI_UNKNOWN);
        std::vector<Scalar> _signed_dist(mesh.tet_vertices.size(), INT_MAX);
        vector_unique(intersection_results_wn);
        int result = tet_subdivision(_f_id, mesh, _intersection_results, intersection_results_wn, cut_f_ids,
                _signed_dist, _oris, e_p_map, tmp_tet_vertices, is_again);
        if (result != SUCCESS) {
            //todo: this is wrong, you need to do real snapping!!!

//            for (const auto &arr_tf: arr_tf_ids) {
//                int t_id = arr_tf[0];
//                int j = arr_tf[1];
//                std::vector<GEO::vec3> ps;
//                sample_triangle({{tet_vertices[tets[t_id][mod4(j + 1)]].pos, tet_vertices[tets[t_id][mod4(j + 2)]].pos,
//                                         tet_vertices[tets[t_id][mod4(j + 3)]].pos}}, ps, mesh.params.dd);
//                if (tree.is_out_sf_envelope(ps, mesh.params.eps_2)) {
//                    is_boundary_preserved[cur_f_id] = false;
//                    break;
//                }
//            }
        }
    }

//    fout.close();
}

void floatTetWild::update_tmp_b_tree(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
        const std::vector<bool> &is_face_inserted, AABBWrapper& tree) {
    std::vector<std::vector<int>> conn_tris(input_vertices.size());
    std::vector<std::array<int, 2>> edges;
    edges.reserve(input_faces.size()*3);
    for (int f_id = 0; f_id < input_faces.size(); f_id++) {
        if (!is_face_inserted[f_id])
            continue;
        for (int j = 0; j < 3; j++) {
            conn_tris[input_faces[f_id][j]].push_back(f_id);
            edges.push_back({{input_faces[f_id][j], input_faces[f_id][mod3(j + 1)]}});
        }
    }
    vector_unique(edges);

    std::vector<std::array<int, 2>> f_id_le_id;
    get_current_open_boundary(input_vertices, input_faces, edges, conn_tris, f_id_le_id);

    std::vector<std::array<int, 2>> b_edges;
    b_edges.reserve(f_id_le_id.size());
    for (auto &info: f_id_le_id) {
        int f_id = info[0];
        int le_id = info[1];
        b_edges.push_back({{input_faces[f_id][le_id], input_faces[f_id][mod3(le_id + 1)]}});
    }
    tree.init_tmp_b_mesh_and_tree(input_vertices, input_faces, b_edges);


//    /////test
//    std::ofstream fout("boundary.obj");
//    for (auto &v:input_vertices) {
//        fout << "v " << v.transpose() << endl;
//    }
//    for (auto &e: b_edges) {
//        fout << "l " << e[0] + 1 << " " << e[1] + 1 << endl;
//    }
//    fout.close();
//    /////test
}

void floatTetWild::get_current_open_boundary(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
        const std::vector<std::array<int, 2>>& edges, const std::vector<std::vector<int>>& conn_tris,
        std::vector<std::array<int, 2>>& f_id_le_id){
    int cnt1 = 0;
    int cnt2 = 0;
    for (const auto &e: edges) {
        int v1_id = e[0];
        int v2_id = e[1];
        std::vector<int> n12_f_ids;
        std::set_intersection(conn_tris[v1_id].begin(), conn_tris[v1_id].end(),
                              conn_tris[v2_id].begin(), conn_tris[v2_id].end(), std::back_inserter(n12_f_ids));
        int f_id = n12_f_ids[0];
        int j = 0;
        for (; j < 3; j++) {
            if (input_faces[f_id][j] == v1_id && input_faces[f_id][mod3(j + 1)] == v2_id
                || input_faces[f_id][j] == v2_id && input_faces[f_id][mod3(j + 1)] == v1_id)
                break;
        }

        if (n12_f_ids.size() == 1) {//open boundary
            f_id_le_id.push_back({{f_id, j}});
            cnt1++;
        } else if (n12_f_ids.size() > 1) {
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
                        if (new_ori != ori) {
                            is_fine = true;
                        }
                        break;
                    }
                }
                if (is_fine)
                    break;
            }
            if (is_fine)
                continue;

            cnt2++;
            for (int ff_id: n12_f_ids) {
                int k = 0;
                for (; k < 3; k++) {
                    if (input_faces[ff_id][k] == v1_id && input_faces[ff_id][mod3(k + 1)] == v2_id
                        || input_faces[ff_id][k] == v2_id && input_faces[ff_id][mod3(k + 1)] == v1_id)
                        break;
                }
                f_id_le_id.push_back({{ff_id, k}});
            }
        }
    }

    cout<<"cnt1 = "<<cnt1<<endl;
    cout<<"cnt2 = "<<cnt2<<endl;
}