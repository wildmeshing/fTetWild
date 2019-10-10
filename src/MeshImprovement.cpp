#include <floattetwild/MeshImprovement.h>

#include <floattetwild/LocalOperations.h>
#include <floattetwild/EdgeSplitting.h>
#include <floattetwild/EdgeCollapsing.h>
#include <floattetwild/EdgeSwapping.h>
#include <floattetwild/VertexSmoothing.h>
#include <floattetwild/Parameters.h>
#include <floattetwild/MeshIO.hpp>

#include <floattetwild/FloatTetCutting.h>
#include <floattetwild/TriangleInsertion.h>
#include <floattetwild/Statistics.h>

#include <floattetwild/Logger.hpp>

#include <igl/Timer.h>

void floatTetWild::init(Mesh &mesh, AABBWrapper& tree) {
    cout << "initializing..." << endl;
//    for (auto &t: mesh.tets) {
//        if (t.is_removed)
//            continue;
//        t.quality = get_quality(mesh, t);
//    }

    for (auto &v: mesh.tet_vertices) {
        if (v.is_removed)
            continue;
        v.is_on_surface = false;
        v.is_on_bbox = false;
//        v.is_on_boundary = false;
    }
    int cnt = 0;
    for (auto &t: mesh.tets) {
        if (t.is_removed)
            continue;
        for (int j = 0; j < 4; j++) {
//            if (t.is_surface_fs[j] != NOT_SURFACE) {
            if (t.is_surface_fs[j] <= 0) {
                mesh.tet_vertices[t[mod4(j + 1)]].is_on_surface = true;
                mesh.tet_vertices[t[mod4(j + 2)]].is_on_surface = true;
                mesh.tet_vertices[t[mod4(j + 3)]].is_on_surface = true;
                cnt++;
            }
            if (t.is_bbox_fs[j] != NOT_BBOX) {
                mesh.tet_vertices[t[mod4(j + 1)]].is_on_bbox = true;
                mesh.tet_vertices[t[mod4(j + 2)]].is_on_bbox = true;
                mesh.tet_vertices[t[mod4(j + 3)]].is_on_bbox = true;
            }
        }
    }

    //todo: after all faces are inserted, mark true is_on_boundary

//    //vertices: mark is_on_boundary
//    std::vector<std::array<int, 2>> b_edges;
//    b_edges.reserve(cnt*3);
//    for (auto &t: mesh.tets) {
//        if (t.is_removed)
//            continue;
//        for (int j = 0; j < 4; j++) {
//            if (t.is_surface_fs[j] < 0) {
//                for (int k = 0; k < 3; k++) {
//                    std::array<int, 2> e = {{t[mod4(j + 1 + k)], t[mod4(j + 1 + mod3(k + 1))]}};
//                    if (e[0] > e[1])
//                        std::swap(e[0], e[1]);
//                    b_edges.push_back(e);
//                }
//            }
//        }
//    }
//    std::sort(b_edges.begin(), b_edges.end());
//
//    bool is_unique = true;
//    for (int i = 0; i < b_edges.size() - 1; i++) {
//        if (b_edges[i] == b_edges[i + 1]) {
//            is_unique = false;
//        } else {
//            if (is_unique) {
//                mesh.tet_vertices[b_edges[i][0]].is_on_boundary = true;
//                mesh.tet_vertices[b_edges[i][1]].is_on_boundary = true;
//            }
//            is_unique = true;
//        }
//    }
//    if (is_unique) {
//        mesh.tet_vertices[b_edges.back()[0]].is_on_boundary = true;
//        mesh.tet_vertices[b_edges.back()[1]].is_on_boundary = true;
//    }
//
//    //rebuild b_tree
//    tree.init_tmp_b_mesh_and_tree(mesh, b_edges);


//    std::ofstream fout("not_bp.obj");
//    for(auto& v:mesh.tet_vertices){
//        if(v.is_removed)
//            continue;
//        if(v.is_on_surface && !v.is_on_boundary)
//            fout<<"v "<<v.pos.transpose()<<endl;
//    }
//    fout.close();

    if (mesh.params.log_level<3) {
        output_surface(mesh, mesh.params.output_path + "_" + mesh.params.postfix + "_cutting");
        output_info(mesh, tree);
        //pausee();
        int v_num, t_num;
        double max_energy, avg_energy;
        v_num = mesh.get_v_num();
        t_num = mesh.get_t_num();
        get_max_avg_energy(mesh, max_energy, avg_energy);
        cout << "#v = " << v_num << endl;
        cout << "#t = " << t_num << endl;
        cout << "max_energy = " << max_energy << endl;
        cout << "avg_energy = " << avg_energy << endl;
    }
}

void floatTetWild::optimization(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces, const std::vector<int> &input_tags, std::vector<bool> &is_face_inserted,
        Mesh &mesh, AABBWrapper& tree, const std::array<int, 4> &ops) {
    init(mesh, tree);

    ////pre-processing
    mesh.is_limit_length = false;
    operation(input_vertices, input_faces, input_tags, is_face_inserted, mesh, tree, std::array<int, 5>({{0, 1, 0, 0, 1}}));
    cleanup_empty_slots(mesh);
    mesh.is_limit_length = true;

    const int M = 5;
    const int N = 5;

    ////optimization
    int it_after_al_inserted = 0;
    bool is_just_after_update = false;
    bool is_hit_min_edge_length = false;
    std::vector<std::array<Scalar, 2>> quality_queue;
    int cnt_increase_epsilon = mesh.params.stage - 1;
    for (int it = 0; it < mesh.params.max_its; it++) {
        if (mesh.is_input_all_inserted)
            it_after_al_inserted++;

        Scalar max_energy, avg_energy;
        get_max_avg_energy(mesh, max_energy, avg_energy);
        if (max_energy <= mesh.params.stop_energy && it_after_al_inserted > M)
            break;

        if (mesh.params.stop_p > 0) {
            int p = get_max_p(mesh);
            cout << "p = " << p << endl;
            if (p <= mesh.params.stop_p && it_after_al_inserted > M)
                break;
        }

        cout << "//////////////// pass " << it << " ////////////////" << endl;
        std::array<int, 5> it_ops;
        if (it % 3 == 2)
            it_ops = {{ops[0], ops[1], ops[2], ops[3], 1}};
        else
            it_ops = {{ops[0], ops[1], ops[2], ops[3], 0}};
        operation(input_vertices, input_faces, input_tags, is_face_inserted, mesh, tree, it_ops);

        if (it > mesh.params.max_its / 4 && max_energy > 1e3) {//Scalar check
            if (cnt_increase_epsilon > 0 && cnt_increase_epsilon == mesh.params.stage - 1) {
//                mesh.params.eps += mesh.params.eps_input / mesh.params.stage;
                mesh.params.eps += mesh.params.eps_delta;
                mesh.params.eps_2 = mesh.params.eps * mesh.params.eps;
                cnt_increase_epsilon--;
                cout << "enlarge envelope, eps = " << mesh.params.eps << endl;
//                pausee();
            }
        }

        Scalar new_max_energy, new_avg_energy;
        get_max_avg_energy(mesh, new_max_energy, new_avg_energy);
        if (!is_just_after_update) {
            if (max_energy - new_max_energy < 5e-1 && (avg_energy - new_avg_energy) / avg_energy < 0.1) {
                is_hit_min_edge_length = update_scaling_field(mesh, new_max_energy) || is_hit_min_edge_length;
                is_just_after_update = true;
                if (cnt_increase_epsilon > 0) {
//                    mesh.params.eps += mesh.params.eps_input / mesh.params.stage;
                    mesh.params.eps += mesh.params.eps_delta;
                    mesh.params.eps_2 = mesh.params.eps * mesh.params.eps;
                    cnt_increase_epsilon--;
                    cout << "enlarge envelope, eps = " << mesh.params.eps << endl;
//                    pausee();
                }
            }
        } else
            is_just_after_update = false;

        quality_queue.push_back(std::array<Scalar, 2>({{new_max_energy, new_avg_energy}}));
        if (is_hit_min_edge_length && it_after_al_inserted > M && it > M + N) {
            bool is_break = true;
            for (int j = 0; j < N; j++) {
                if (quality_queue[it - j][0] - quality_queue[it - j - 1][0] < 0) {
                    is_break = false;
                    break;
                }
            }
            if (is_break)
                break;

//            bool is_loop = true;
//            for (int i = 0; i < M; i++) {
//                if (quality_queue[it - i][0] - quality_queue[it - i - 1][0] < -1e-4
//                    || quality_queue[it - i][1] - quality_queue[it - i - 1][1] < -1e-4) {
//                    is_loop = false;
//                    break;
//                }
//            }
//            if (is_loop)
//                break;
        }
    }

    ////postprocessing
    cout << "//////////////// postprocessing ////////////////" << endl;
    for (auto &v:mesh.tet_vertices) {
        if (v.is_removed)
            continue;
        v.sizing_scalar = 1;
    }
    operation(input_vertices, input_faces, input_tags, is_face_inserted, mesh, tree, std::array<int, 5>({{0, 1, 0, 0, 0}}));
}

void floatTetWild::cleanup_empty_slots(Mesh &mesh, double percentage) {
    if (mesh.tets.size() < 9e5)
        return;
    cout<<mesh.tets.size()<<" ==> ";
    ///
    const int v_end_id = mesh.tet_vertices.size() * percentage;
    const int t_end_id = mesh.tets.size() * percentage;
    //
    std::vector<int> map_v_ids(mesh.tet_vertices.size(), -1);
    int cnt = 0;
    for (int i = 0; i < v_end_id; i++) {
        if (mesh.tet_vertices[i].is_removed)
            continue;
        map_v_ids[i] = cnt++;
    }
    for (int i = v_end_id; i < mesh.tet_vertices.size(); i++) {
        map_v_ids[i] = cnt++;
    }
    //
    std::vector<int> map_t_ids(mesh.tets.size(), -1);
    cnt = 0;
    for (int i = 0; i < t_end_id; i++) {
        if (mesh.tets[i].is_removed)
            continue;
        map_t_ids[i] = cnt++;
    }
    for (int i = t_end_id; i < mesh.tets.size(); i++) {
        map_t_ids[i] = cnt++;
    }

    ///
    mesh.tet_vertices.erase(std::remove_if(mesh.tet_vertices.begin(), mesh.tet_vertices.begin() + v_end_id,
                                           [](const MeshVertex &v) { return v.is_removed; }),
                            mesh.tet_vertices.begin() + v_end_id);
    mesh.tets.erase(std::remove_if(mesh.tets.begin(), mesh.tets.begin() + t_end_id,
                                   [](const MeshTet &t) { return t.is_removed; }),
                    mesh.tets.begin() + t_end_id);

    ///
    for (auto &v: mesh.tet_vertices) {
        if (v.is_removed)
            continue;
        for (auto &t_id: v.conn_tets)
            t_id = map_t_ids[t_id];
    }
    for (auto &t: mesh.tets) {
        if (t.is_removed)
            continue;
        for (int j = 0; j < 4; j++)
            t[j] = map_v_ids[t[j]];
    }
    cout<<mesh.tets.size()<<endl;
}

void floatTetWild::operation(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces, const std::vector<int> &input_tags, std::vector<bool> &is_face_inserted,
        Mesh &mesh, AABBWrapper& tree, const std::array<int, 5> &ops) {
    igl::Timer igl_timer;
    int v_num, t_num;
    double max_energy, avg_energy;
    double time;

    for (int i = 0; i < ops[0]; i++) {
        igl_timer.start();
        cout << "edge splitting..." << endl;
        edge_splitting(mesh);
        time = igl_timer.getElapsedTime();
        cout << "edge splitting done!" << endl;
        cout << "time = " << time << "s" << endl;
        v_num = mesh.get_v_num();
        t_num = mesh.get_t_num();
        get_max_avg_energy(mesh, max_energy, avg_energy);
        cout << "#v = " << v_num << endl;
        cout << "#t = " << t_num << endl;
        cout << "max_energy = " << max_energy << endl;
        cout << "avg_energy = " << avg_energy << endl;
        stats().record(StateInfo::splitting_id, time, v_num, t_num, max_energy, avg_energy);
        output_info(mesh, tree);
    }

    for (int i = 0; i < ops[1]; i++) {
        igl_timer.start();
        cout << "edge collapsing..." << endl;
        edge_collapsing(mesh, tree);
        time = igl_timer.getElapsedTime();
        cout << "edge collapsing done!" << endl;
        cout << "time = " << time << "s" << endl;
        v_num = mesh.get_v_num();
        t_num = mesh.get_t_num();
        get_max_avg_energy(mesh, max_energy, avg_energy);
        cout << "#v = " << v_num << endl;
        cout << "#t = " << t_num << endl;
        cout << "max_energy = " << max_energy << endl;
        cout << "avg_energy = " << avg_energy << endl;
        stats().record(StateInfo::collapsing_id, time, v_num, t_num, max_energy, avg_energy);
        output_info(mesh, tree);
    }

    for (int i = 0; i < ops[2]; i++) {
        igl_timer.start();
        cout << "edge swapping..." << endl;
        edge_swapping(mesh);
        time = igl_timer.getElapsedTime();
        cout << "edge swapping done!" << endl;
        cout << "time = " << time << "s" << endl;
        v_num = mesh.get_v_num();
        t_num = mesh.get_t_num();
        get_max_avg_energy(mesh, max_energy, avg_energy);
        cout << "#v = " << v_num << endl;
        cout << "#t = " << t_num << endl;
        cout << "max_energy = " << max_energy << endl;
        cout << "avg_energy = " << avg_energy << endl;
        stats().record(StateInfo::swapping_id, time, v_num, t_num, max_energy, avg_energy);
        output_info(mesh, tree);
    }

    for (int i = 0; i < ops[3]; i++) {
        igl_timer.start();
        cout << "vertex smoothing..." << endl;
        vertex_smoothing(mesh, tree);
        time = igl_timer.getElapsedTime();
        cout << "vertex smoothing done!" << endl;
        cout << "time = " << time << "s" << endl;
        v_num = mesh.get_v_num();
        t_num = mesh.get_t_num();
        get_max_avg_energy(mesh, max_energy, avg_energy);
        cout << "#v = " << v_num << endl;
        cout << "#t = " << t_num << endl;
        cout << "max_energy = " << max_energy << endl;
        cout << "avg_energy = " << avg_energy << endl;
        stats().record(StateInfo::smoothing_id, time, v_num, t_num, max_energy, avg_energy);
        output_info(mesh, tree);
    }

    if (!mesh.is_input_all_inserted) {
        pausee();

        for (int i = 0; i < ops[4]; i++) {
            //reset boundary points
            for (auto &v: mesh.tet_vertices) {
                if (v.is_removed)
                    continue;
                v.is_on_boundary = false;
                v.on_boundary_e_id = -1;
            }
            //
            igl_timer.start();
            insert_triangles(input_vertices, input_faces, input_tags, mesh, is_face_inserted, tree, true);
            init(mesh, tree);
            stats().record(StateInfo::cutting_id, igl_timer.getElapsedTimeInSec(),
                           mesh.get_v_num(), mesh.get_t_num(),
                           mesh.get_max_energy(), mesh.get_avg_energy(),
                           std::count(is_face_inserted.begin(), is_face_inserted.end(),
                                      false));
        }

//        for (int i = 0; i < ops[4]; i++) {
//            if(!mesh.is_input_all_inserted) {
//                //check isolate boundary points
//                for (auto &v: mesh.tet_vertices) {
//                    if (v.is_removed || !v.is_on_boundary)
//                        continue;
//                    if (is_point_out_boundary_envelope(mesh, v.pos, tree))
//                        v.is_on_boundary = false;
//                }
//            }
//
////            std::ofstream fout("bp_del.obj");
////            for(auto& v:mesh.tet_vertices){
////                if(v.is_removed)
////                    continue;
////                if(v.is_on_boundary)
////                    fout<<"v "<<v.pos.transpose()<<endl;
////            }
////            fout.close();
//
////            pausee();
//            igl_timer.start();
//
//            mesh.reset_t_empty_start();
//            mesh.reset_v_empty_start();
//
//            for (int t_id = 0; t_id < mesh.tets.size(); t_id++) {
//                if (mesh.tets[t_id].is_removed)
//                    continue;
//                mesh.tets[t_id].opp_t_ids = {{OPP_T_ID_UNKNOWN, OPP_T_ID_UNKNOWN, OPP_T_ID_UNKNOWN, OPP_T_ID_UNKNOWN}};
//            }
//
//            //todo: keep track of opp_t_ids
////        for (int t_id = 0; t_id < mesh.tets.size(); t_id++) {
////            if (mesh.tets[t_id].is_removed)
////                continue;
////            for (int j = 0; j < 4; j++)
////                mesh.tets[t_id].opp_t_ids[j] = -1;
////        }
////
////        for (int t_id = 0; t_id < mesh.tets.size(); t_id++) {
////            if (mesh.tets[t_id].is_removed)
////                continue;
////
////            for (int j = 0; j < 4; j++) {
////                if (mesh.tets[t_id].opp_t_ids[j] >= 0)
////                    continue;
////
////                int opp_t_id = get_opp_t_id(mesh, t_id, j);
////                if (opp_t_id < 0)
////                    continue;
////                mesh.tets[t_id].opp_t_ids[j] = opp_t_id;
////                for (int k = 0; k < 4; k++) {
////                    if (mesh.tets[opp_t_id][k] != mesh.tets[t_id][(j + 1) % 4]
////                        && mesh.tets[opp_t_id][k] != mesh.tets[t_id][(j + 2) % 4]
////                        && mesh.tets[opp_t_id][k] != mesh.tets[t_id][(j + 3) % 4])
////                        mesh.tets[opp_t_id].opp_t_ids[k] = t_id;
////                }
////            }
////        }
//
//            cutting(input_vertices, input_faces, input_tags, mesh, is_face_inserted, tree, true);
//            init(mesh, tree);
//
//            stats().record(StateInfo::cutting_id, igl_timer.getElapsedTimeInSec(),
//                                                           mesh.get_v_num(), mesh.get_t_num(),
//                                                           mesh.get_max_energy(), mesh.get_avg_energy(),
//                                                           std::count(is_face_inserted.begin(), is_face_inserted.end(),
//                                                                      false));
//        }
    }
}

#include <geogram/points/kd_tree.h>
bool floatTetWild::update_scaling_field(Mesh &mesh, Scalar max_energy) {
    cout << "updating sclaing field ..." << endl;
    bool is_hit_min_edge_length = false;

    Scalar radius0 = mesh.params.ideal_edge_length * 1.8;//increasing the radius would increase the #v in output
//    if(is_hit_min)
//        radius0 *= 2;

    static const Scalar stop_filter_energy = mesh.params.stop_energy * 0.8;

    Scalar filter_energy = max_energy / 100 > stop_filter_energy ? max_energy / 100 : stop_filter_energy;
//    if(filter_energy > 1e3)
//        filter_energy = 1e3;

    if (filter_energy > 100) {
//        filter_energy = get_mid_energy(mesh);
//        if (filter_energy < stop_filter_energy)
//            filter_energy = stop_filter_energy;
        filter_energy = 100;
    }

    cout << "filter_energy = " << filter_energy << endl;
    Scalar recover = 1.5;
    std::vector<Scalar> scale_multipliers(mesh.tet_vertices.size(), recover);
    Scalar refine_scale = 0.5;
//    Scalar min_refine_scale = mesh.epsilon / mesh.ideal_edge_length;
    Scalar min_refine_scale = mesh.params.min_edge_len_rel;

    const int N = -int(std::log2(min_refine_scale) - 1);
    std::vector<std::vector<int>> v_ids(N, std::vector<int>());
    for (int i = 0; i < mesh.tet_vertices.size(); i++) {
        auto &v = mesh.tet_vertices[i];
        if (v.is_removed)
            continue;

        bool is_refine = false;
        for (int t_id: v.conn_tets) {
            if (mesh.tets[t_id].quality > filter_energy)
                is_refine = true;
        }
        if (!is_refine)
            continue;

        int n = -int(std::log2(v.sizing_scalar) - 0.5);
        if (n >= N)
            n = N - 1;
        v_ids[n].push_back(i);
    }

    for (int n = 0; n < N; n++) {
        if (v_ids[n].size() == 0)
            continue;

        Scalar radius = radius0 / std::pow(2, n);

        std::unordered_set<int> is_visited;
        std::queue<int> v_queue;

        std::vector<double> pts;//geogram needs double []
        pts.reserve(v_ids[n].size() * 3);
        for (int i = 0; i < v_ids[n].size(); i++) {
            pts.push_back(mesh.tet_vertices[v_ids[n][i]].pos[0]);
            pts.push_back(mesh.tet_vertices[v_ids[n][i]].pos[1]);
            pts.push_back(mesh.tet_vertices[v_ids[n][i]].pos[2]);

            v_queue.push(v_ids[n][i]);
            is_visited.insert(v_ids[n][i]);
            scale_multipliers[v_ids[n][i]] = refine_scale;
        }
        // construct the kdtree
        GEO::NearestNeighborSearch_var nnsearch = GEO::NearestNeighborSearch::create(3, "BNN");
        nnsearch->set_points(int(v_ids[n].size()), pts.data());

        while (!v_queue.empty()) {
            int v_id = v_queue.front();
            v_queue.pop();

            for (int t_id:mesh.tet_vertices[v_id].conn_tets) {
                for (int j = 0; j < 4; j++) {
                    if (is_visited.find(mesh.tets[t_id][j]) != is_visited.end())
                        continue;
                    GEO::index_t _;
                    double sq_dist;
                    const double p[3] = {mesh.tet_vertices[mesh.tets[t_id][j]].pos[0],
                                         mesh.tet_vertices[mesh.tets[t_id][j]].pos[1],
                                         mesh.tet_vertices[mesh.tets[t_id][j]].pos[2]};
                    nnsearch->get_nearest_neighbors(1, p, &_, &sq_dist);
                    Scalar dis = sqrt(sq_dist);

                    if (dis < radius) {
                        v_queue.push(mesh.tets[t_id][j]);
                        Scalar new_ss = (dis / radius) * (1 - refine_scale) + refine_scale;
                        if (new_ss < scale_multipliers[mesh.tets[t_id][j]])
                            scale_multipliers[mesh.tets[t_id][j]] = new_ss;
                    }
                    is_visited.insert(mesh.tets[t_id][j]);
                }
            }
        }
    }

    // update scalars
    for (int i=0;i< mesh.tet_vertices.size();i++) {
        auto& v = mesh.tet_vertices[i];
        if (v.is_removed)
            continue;
        Scalar new_scale = v.sizing_scalar * scale_multipliers[i];
        if (new_scale > 1)
            v.sizing_scalar = 1;
//        if (new_scale > mesh.tri_vertices[i].max_scale)
//            mesh.tri_vertices[i].scale = mesh.tri_vertices[i].max_scale;
        else if (new_scale < min_refine_scale) {
            is_hit_min_edge_length = true;
            v.sizing_scalar = min_refine_scale;
        } else
            v.sizing_scalar = new_scale;
    }

    cout << "is_hit_min_edge_length = " << is_hit_min_edge_length << endl;
    return is_hit_min_edge_length;
}

void floatTetWild::output_info(Mesh& mesh, const AABBWrapper& tree) {
    if(mesh.params.is_quiet)
        return;

    if(mesh.params.log_level >= 2)
        return;

    auto &tets = mesh.tets;
    auto &tet_vertices = mesh.tet_vertices;

    //count
    int cnt_v = mesh.get_v_num();
    int cnt_t = mesh.get_t_num();
//    cout << "#v = " << cnt_v << "(" << tet_vertices.size() << ")" << endl;
//    cout << "#t = " << cnt_t << "(" << tets.size() << ")" << endl;

    //quality
//    Scalar max_energy, avg_energy;
//    get_max_avg_energy(mesh, max_energy, avg_energy);
//    cout << "max_energy = " << max_energy << endl;
//    cout << "avg_energy = " << avg_energy << endl;

//    for (int i = 0; i < tets.size(); i++) {
//        if (tets[i].is_removed)
//            continue;
//        Scalar q = get_quality(mesh, i);
//        if (abs(tets[i].quality - q)/tets[i].quality > 0.01) {
//            cout << "tets[i].quality != get_quality(mesh,i)" << endl;
//            cout << tets[i].quality << " - " << q << " = " << tets[i].quality - q << endl;
//            //pausee();
//        }
//    }

    if(mesh.params.log_level > 1) {
        output_surface(mesh, mesh.params.output_path+"_"+mesh.params.postfix+"_opt");
        return;
    }

    //euler
    std::vector<std::array<int, 2>> edges;
    get_all_edges(mesh, edges);
    std::vector<std::array<int, 3>> faces;
    for (auto &t: tets) {
        if(t.is_removed)
            continue;
        for (int j = 0; j < 4; j++) {
            std::array<int, 3> f = {{t[j], t[mod4(j + 1)], t[mod4(j + 2)]}};
            std::sort(f.begin(), f.end());
            faces.push_back(f);
        }
    }
    vector_unique(faces);
    int euler = cnt_v - edges.size() + faces.size() - cnt_t;
    if (euler != 1) {
        cout << "euler error " << euler << endl;
        //pausee();
    }

    //inversion
    for (int i = 0; i < tets.size(); i++) {
        if (tets[i].is_removed)
            continue;
        for (int j = 0; j < 4; j++){
            if(tet_vertices[tets[i][j]].is_removed){
                cout<<"tet_vertices[tets[i][j]].is_removed"<<endl;
//                //pausee();
            }
        }
    }
    for (int i = 0; i < tets.size(); i++) {
        if (tets[i].is_removed)
            continue;
        if (is_inverted(mesh, i)) {
            cout << "tet " << i << " inverted: " << tets[i][0] << " " << tets[i][1] << " " << tets[i][2] << " "
                 << tets[i][3] << endl;
            //pausee();
        }
    }

    //conn_tets
//    for (int i = 0; i < tet_vertices.size(); i++) {
//        if (tet_vertices[i].is_removed)
//            continue;
//        int old_size = tet_vertices[i].conn_tets.size();
//        vector_unique(tet_vertices[i].conn_tets);
//        if (tet_vertices[i].conn_tets.size() != old_size) {
//            cout << "tet_vertices[i].conn_tets.size()!=old_size()" << endl;
//            cout << tet_vertices[i].conn_tets.size() << " " << old_size << endl;
//        }
//    }
//    for (int i = 0; i < tets.size(); i++) {
//        if (tets[i].is_removed)
//            continue;
//        for (int j = 0; j < 4; j++) {
////            if (tet_vertices[tets[i][j]].conn_tets.find(i) == tet_vertices[tets[i][j]].conn_tets.end()) {
//            if (std::find(tet_vertices[tets[i][j]].conn_tets.begin(), tet_vertices[tets[i][j]].conn_tets.end(), i)
//                == tet_vertices[tets[i][j]].conn_tets.end()) {
//                cout << "conn_tets error!" << endl;
//                //pausee();
//            }
//        }
//    }
//    for (int i = 0; i < tet_vertices.size(); i++) {
//        if (tet_vertices[i].is_removed)
//            continue;
//        for(int t_id:tet_vertices[i].conn_tets){
//            if(t_id>tets.size() || t_id<0){
//                cout<<"t_id>tets.size() || t_id<0"<<endl;
//                //pausee();
//            }
//            int j = tets[t_id].find(i);
//            if(j<0){
//                cout<<"conn_tets error: j<0"<<endl;
//                //pausee();
//            }
//        }
//    }

    //check conn_tets
    std::vector<std::vector<int>> tmp_conn_tets(tet_vertices.size());
    for (int i = 0; i < tets.size(); i++) {
        if(tets[i].is_removed)
            continue;
        for (int j = 0; j < 4; j++)
            tmp_conn_tets[tets[i][j]].push_back(i);
    }
    for (int i = 0; i < tet_vertices.size(); i++) {
        if(tet_vertices[i].is_removed)
            continue;
        std::vector<int> conn_tets = tet_vertices[i].conn_tets;
        std::sort(conn_tets.begin(), conn_tets.end());
        if (conn_tets != tmp_conn_tets[i]) {
            cout << "conn_tets error" << endl;
            for (auto &ii:tet_vertices[i].conn_tets)
                cout << ii << " ";
            cout << endl;
            for (auto &ii:tmp_conn_tets[i])
                cout << ii << " ";
            cout << endl;
            //pausee();
        }
    }

//    auto get_opp_t_id = [&](int t_id, int j) {
//        std::vector<int> p;
//        set_intersection(tet_vertices[tets[t_id][(j + 1) % 4]].conn_tets,
//                         tet_vertices[tets[t_id][(j + 2) % 4]].conn_tets,
//                         tet_vertices[tets[t_id][(j + 3) % 4]].conn_tets, p);
//        if (p.size() < 2)
//            return OPP_T_ID_BOUNDARY;
//        return p[0] == t_id ? p[1] : p[0];
//    };

    //surface tracking
    for (int i = 0; i < tets.size(); i++) {
        if (tets[i].is_removed)
            continue;
        for (int j = 0; j < 4; j++) {
            if (tets[i].is_surface_fs[j] != NOT_SURFACE && tets[i].is_bbox_fs[j] != NOT_BBOX){
                cout<<"tets[i].is_surface_fs[j] != NOT_SURFACE && tets[i].is_bbox_fs[j] != NOT_BBOX"<<endl;
                cout<<i<<" "<<j<<endl;
                //pausee();
            }
            if (tets[i].is_surface_fs[j] != NOT_SURFACE) {
//                if (tets[i].is_surface_fs[j] == 0){
//                    cout << "is_surface_fs error 0" << endl;
//                    cout << "t " << i << ": " << tets[i][0] << " " << tets[i][1] << " " << tets[i][2] << " "
//                         << tets[i][3] << endl;
//                    cout << tets[i].is_surface_fs[0] << " " << tets[i].is_surface_fs[1] << " "
//                         << tets[i].is_surface_fs[2] << " " << tets[i].is_surface_fs[3] << endl;
//                    cout << tet_vertices[tets[i][0]].is_on_surface << " " << tet_vertices[tets[i][1]].is_on_surface
//                         << " "
//                         << tet_vertices[tets[i][2]].is_on_surface << " " << tet_vertices[tets[i][3]].is_on_surface
//                         << endl;
//                    //pausee();
//                }

                for (int k = 0; k < 3; k++) {
                    if (!tet_vertices[tets[i][(j + 1 + k) % 4]].is_on_surface) {
                        cout << "is_surface_fs error" << endl;
                        cout << "t " << i << ": " << tets[i][0] << " " << tets[i][1] << " " << tets[i][2] << " "
                             << tets[i][3] << endl;
                        cout << tets[i].is_surface_fs[0] << " " << tets[i].is_surface_fs[1] << " "
                             << tets[i].is_surface_fs[2] << " " << tets[i].is_surface_fs[3] << endl;
                        cout << tet_vertices[tets[i][0]].is_on_surface << " " << tet_vertices[tets[i][1]].is_on_surface
                             << " "
                             << tet_vertices[tets[i][2]].is_on_surface << " " << tet_vertices[tets[i][3]].is_on_surface
                             << endl;
                        //pausee();
                    }
                }
            }

            if (tets[i].is_bbox_fs[j] != NOT_BBOX) {
                for (int k = 0; k < 3; k++) {
                    if (!tet_vertices[tets[i][(j + 1 + k) % 4]].is_on_bbox) {
                        cout<<"is_bbox_fs error"<<endl;
                        cout << "t " << i << ": " << tets[i][0] << " " << tets[i][1] << " " << tets[i][2] << " "
                             << tets[i][3] << endl;
                        cout << (int)tets[i].is_bbox_fs[0] << " " << (int)tets[i].is_bbox_fs[1] << " "
                             << (int)tets[i].is_bbox_fs[2] << " " << (int)tets[i].is_bbox_fs[3] << endl;
                        cout << tet_vertices[tets[i][0]].is_on_bbox << " " << tet_vertices[tets[i][1]].is_on_bbox << " "
                             << tet_vertices[tets[i][2]].is_on_bbox << " " << tet_vertices[tets[i][3]].is_on_bbox
                             << endl;
                        pausee();
                    }
                }
                if(get_opp_t_id(mesh, i, j) != OPP_T_ID_BOUNDARY){
                    cout<<"wrong-marked bbox face"<<endl;
                    //pausee();
                }
            } else {
                if(get_opp_t_id(mesh, i, j) == OPP_T_ID_BOUNDARY){
                    cout<<"unmarked bbox face"<<endl;
                    //pausee();
                }
            }
        }
    }

    for(int i=0;i<tet_vertices.size();i++) {
        if (tet_vertices[i].is_removed)
            continue;

        if(tet_vertices[i].is_on_bbox && tet_vertices[i].is_on_surface){
            cout<<"error tet_vertices[i].is_on_bbox && tet_vertices[i].is_on_surface"<<endl;
            cout<<"v "<<i<<endl;
            //pausee();
        }

        if (tet_vertices[i].is_on_bbox) {
            bool is_found = false;
            for (int t_id:tet_vertices[i].conn_tets) {
                int j = tets[t_id].find(i);
                for (int k = 0; k < 3; k++) {
                    if (tets[t_id].is_bbox_fs[(j + 1 + k) % 4] != NOT_BBOX) {
                        is_found = true;
                        break;
                    }
                }
            }
            if (!is_found) {
                cout << "is_on_bbox error" << endl;
                for (int t_id:tet_vertices[i].conn_tets) {
                    cout << "t " << t_id << ": " << tets[t_id][0] << " " << tets[t_id][1] << " " << tets[t_id][2] << " "
                         << tets[t_id][3] << endl;
                    cout << tets[t_id].is_bbox_fs[0] << " " << tets[t_id].is_bbox_fs[1] << " "
                         << tets[t_id].is_bbox_fs[2] << " " << tets[t_id].is_bbox_fs[3] << endl;
                    cout << tet_vertices[tets[t_id][0]].is_on_bbox << " " << tet_vertices[tets[t_id][1]].is_on_bbox << " "
                         << tet_vertices[tets[t_id][2]].is_on_bbox << " " << tet_vertices[tets[t_id][3]].is_on_bbox
                         << endl;
                }
                //pausee();
            }
        }
        if (tet_vertices[i].is_on_surface) {
            bool is_found = false;
            for (int t_id:tet_vertices[i].conn_tets) {
                int j = tets[t_id].find(i);
                for (int k = 0; k < 3; k++) {
                    if (tets[t_id].is_bbox_fs[(j + 1 + k) % 4] != NOT_SURFACE) {
                        is_found = true;
                        break;
                    }
                }
            }
            if (!is_found) {
                cout << "is_on_surface error" << endl;
                for (int t_id:tet_vertices[i].conn_tets) {
                    cout << "t " << t_id << ": " << tets[t_id][0] << " " << tets[t_id][1] << " " << tets[t_id][2] << " "
                         << tets[t_id][3] << endl;
                    cout << tets[t_id].is_surface_fs[0] << " " << tets[t_id].is_surface_fs[1] << " "
                         << tets[t_id].is_surface_fs[2] << " " << tets[t_id].is_surface_fs[3] << endl;
                    cout << tet_vertices[tets[t_id][0]].is_on_surface << " " << tet_vertices[tets[t_id][1]].is_on_surface
                         << " "
                         << tet_vertices[tets[t_id][2]].is_on_surface << " " << tet_vertices[tets[t_id][3]].is_on_surface
                         << endl;
                }
                //pausee();
            }
        }
    }
    cout<<endl;

    check_envelope(mesh, tree);

//    MeshIO::write_mesh("test.msh", mesh);
    output_surface(mesh, mesh.params.output_path+"_"+mesh.params.postfix+"_opt");
//    //pausee();

    return;
}

void floatTetWild::check_envelope(Mesh& mesh, const AABBWrapper& tree) {
//    if (mesh.params.log_level >= 1)
//        return;

//    Scalar check_eps = mesh.params.eps_input * mesh.params.eps_input;
    Scalar check_eps = mesh.params.eps_2;

    for (auto &t: mesh.tets) {
        if (t.is_removed)
            continue;
        for (int j = 0; j < 4; j++) {
            if (t.is_surface_fs[j] <= 0) {
                std::vector<GEO::vec3> ps;
                sample_triangle({{mesh.tet_vertices[t[(j + 1) % 4]].pos, mesh.tet_vertices[t[(j + 2) % 4]].pos,
                                         mesh.tet_vertices[t[(j + 3) % 4]].pos}}, ps, mesh.params.dd);
                if(tree.is_out_sf_envelope(ps, mesh.params.eps_2)){
//                Scalar d = tree.dist_sf_envelope(ps, check_eps);
//                if (d > mesh.params.eps_2) {
                    cout << "out of envelope!" << endl;
//                    cout << d << ", eps_input = " << check_eps << endl;
//                    //pausee();
                }
            }
        }
    }

    for (auto &v: mesh.tet_vertices) {
        if (v.is_removed || !v.is_on_surface)
            continue;

        std::vector<GEO::vec3> ps = {GEO::vec3(v.pos[0], v.pos[1], v.pos[2])};
        Scalar d = tree.dist_sf_envelope(ps, check_eps);
        if (d > mesh.params.eps_2) {
            cout << "v out of envelope!" << endl;
            cout << d << ", eps_input = " << check_eps << endl;
            //pausee();
        }
    }

    cout<<"envelope check done"<<endl;
}

int floatTetWild::get_max_p(const Mesh &mesh)
{
    const Scalar scaling = 1.0 / (mesh.params.bbox_max - mesh.params.bbox_min).maxCoeff();
    const Scalar B = 3;
    const int p_ref = 1;

    std::vector<std::array<int, 2>> edges;
    get_all_edges(mesh, edges);
    Scalar h_ref = 0;
    //mesh.params.ideal_edge_length;
    for(const auto &e : edges){
        const Vector3 edge = (mesh.tet_vertices[e[0]].pos - mesh.tet_vertices[e[1]].pos)*scaling;
        h_ref += edge.norm();
    }
    h_ref /= edges.size();

    const Scalar rho_ref = sqrt(6.)/12.*h_ref;
    const Scalar sigma_ref = rho_ref / h_ref;

    int max_p = 1;

    for(const auto &t : mesh.tets)
    {
        if(t.is_removed)
            continue;

        const auto &v0 = mesh.tet_vertices[t[0]].pos;
        const auto &v1 = mesh.tet_vertices[t[1]].pos;
        const auto &v2 = mesh.tet_vertices[t[2]].pos;
        const auto &v3 = mesh.tet_vertices[t[3]].pos;

        Eigen::Matrix<Scalar, 6, 3> e;
        e.row(0) = (v0 - v1) * scaling;
        e.row(1) = (v1 - v2) * scaling;
        e.row(2) = (v2 - v0) * scaling;

        e.row(3) = (v0 - v3) * scaling;
        e.row(4) = (v1 - v3) * scaling;
        e.row(5) = (v2 - v3) * scaling;

        const Eigen::Matrix<Scalar, 6, 1> en = e.rowwise().norm();

        const Scalar S = (e.row(0).cross(e.row(1)).norm() + e.row(0).cross(e.row(4)).norm() + e.row(4).cross(e.row(1)).norm() + e.row(2).cross(e.row(5)).norm()) / 2;
        const Scalar V = std::abs(e.row(3).dot(e.row(2).cross(-e.row(0))))/6;
        const Scalar rho = 3 * V / S;
        const Scalar h = en.maxCoeff();

        const Scalar sigma = rho / h;

        const Scalar ptmp = (std::log(B*std::pow(h_ref, p_ref + 1)*sigma*sigma/sigma_ref/sigma_ref) - std::log(h))/std::log(h);
        const int p = (int)std::round(ptmp);
        max_p = std::max(p, max_p);
    }

    return max_p;
}

#include <igl/writeOBJ.h>
#include <igl/writeSTL.h>
#include <floattetwild/Predicates.hpp>

void floatTetWild::output_surface(Mesh& mesh, const std::string& filename) {
#define SF_CONDITION t.is_surface_fs[j]<=0
//#define SF_CONDITION t.is_surface_fs[j]!=NOT_SURFACE
//#define SF_CONDITION t.is_surface_fs[j]<=0&&t.surface_tags[j]==2
//#define SF_CONDITION t.is_bbox_fs[j]==2
//#define SF_CONDITION t.is_bbox_fs[j]!=NOT_BBOX
//#define SF_CONDITION t.opp_t_ids[j]<0

    auto &tets = mesh.tets;
    auto &tet_vertices = mesh.tet_vertices;

    int cnt = 0;
    for (auto &t: mesh.tets) {
        if (t.is_removed)
            continue;
        for (int j = 0; j < 4; j++) {
            if (SF_CONDITION)
                cnt++;
        }
    }

    Eigen::Matrix<Scalar, Eigen::Dynamic, 3> V_sf(cnt * 3, 3);
    Eigen::MatrixXi F_sf(cnt, 3);
    cnt = 0;
    for (int t_id = 0;t_id<mesh.tets.size();t_id++) {
        auto &t = mesh.tets[t_id];
        if (t.is_removed)
            continue;
        for (int j = 0; j < 4; j++) {
            if (SF_CONDITION) {
                for (int k = 0; k < 3; k++)
                    V_sf.row(cnt * 3 + k) = tet_vertices[t[(j + k + 1) % 4]].pos;
                if (Predicates::orient_3d(tet_vertices[t[(j + 1) % 4]].pos, tet_vertices[t[(j + 2) % 4]].pos,
                                          tet_vertices[t[(j + 3) % 4]].pos, tet_vertices[t[j]].pos) ==
                    Predicates::ORI_POSITIVE)
                    F_sf.row(cnt) << cnt * 3, cnt * 3 + 2, cnt * 3 + 1;
                else
                    F_sf.row(cnt) << cnt * 3, cnt * 3 + 1, cnt * 3 + 2;
                cnt++;
            }
        }
    }
    igl::writeSTL(filename + ".stl", Eigen::MatrixXd(V_sf), Eigen::MatrixXi(F_sf));
}

#include <floattetwild/bfs_orient.h>
#include <igl/unique_rows.h>
#include <igl/remove_duplicate_vertices.h>
#include <floattetwild/TriangleInsertion.h>
void floatTetWild::get_tracked_surface(Mesh& mesh, Eigen::Matrix<Scalar, Eigen::Dynamic, 3> &V_sf, Eigen::Matrix<int, Eigen::Dynamic, 3> &F_sf, int c_id) {
#define SF_CONDITION t.is_surface_fs[j]<=0&&t.surface_tags[j]==c_id

    auto &tets = mesh.tets;
    auto &tet_vertices = mesh.tet_vertices;

    int cnt = 0;
    for (auto &t: mesh.tets) {
        if (t.is_removed)
            continue;
        for (int j = 0; j < 4; j++) {
            if (SF_CONDITION)
                cnt++;
        }
    }

    V_sf.resize(cnt * 3, 3);
    F_sf.resize(cnt, 3);
    cnt = 0;
    for (auto &t: mesh.tets) {
        if (t.is_removed)
            continue;
        for (int j = 0; j < 4; j++) {
            if (SF_CONDITION) {
                for (int k = 0; k < 3; k++)
                    V_sf.row(cnt * 3 + k) = tet_vertices[t[mod4(j + k + 1)]].pos;
                if (Predicates::orient_3d(tet_vertices[t[mod4(j + 1)]].pos, tet_vertices[t[mod4(j + 2)]].pos,
                                          tet_vertices[t[mod4(j + 3)]].pos, tet_vertices[t[j]].pos) ==
                    Predicates::ORI_POSITIVE)
                    F_sf.row(cnt) << cnt * 3, cnt * 3 + 2, cnt * 3 + 1;
                else
                    F_sf.row(cnt) << cnt * 3, cnt * 3 + 1, cnt * 3 + 2;
                cnt++;
            }
        }
    }
//    igl::writeSTL("before_bfs.stl", V_sf, F_sf);

    if(mesh.params.correct_surface_orientation) {
        Eigen::MatrixXd V;
        Eigen::MatrixXi F;
        Eigen::VectorXi _1, _2;
        igl::remove_duplicate_vertices(V_sf, F_sf, -1, V, _1, _2, F);
        V_sf = V;
        F_sf.resize(0, 0);
        bfs_orient(F, F_sf, _1);
    }
    igl::writeSTL(mesh.params.output_path + "_tracked_surface.stl", V_sf, F_sf);
}

void floatTetWild::correct_tracked_surface_orientation(Mesh &mesh, AABBWrapper& tree){
    std::vector<std::array<bool, 4>> is_visited(mesh.tets.size(), {{false, false, false, false}});
    for (int t_id = 0; t_id < mesh.tets.size(); t_id++) {
        auto &t = mesh.tets[t_id];
        if (t.is_removed)
            continue;
        for (int j = 0; j < 4; j++) {
//            if (t.is_surface_fs[j] != 0)
            if (t.is_surface_fs[j] == NOT_SURFACE || is_visited[t_id][j])
                continue;
            is_visited[t_id][j] = true;
            int opp_t_id = get_opp_t_id(t_id, j, mesh);
            if (opp_t_id < 0) {
                t.is_surface_fs[j] = NOT_SURFACE;
                continue;
            }
            int k = get_local_f_id(opp_t_id, t[(j + 1) % 4], t[(j + 2) % 4], t[(j + 3) % 4], mesh);
            is_visited[opp_t_id][k] = true;
            //
            Vector3 c = (mesh.tet_vertices[t[(j + 1) % 4]].pos + mesh.tet_vertices[t[(j + 2) % 4]].pos +
                         mesh.tet_vertices[t[(j + 3) % 4]].pos) / 3;
            int f_id = tree.get_nearest_face_sf(c);
            const auto &fv1 = tree.sf_mesh.vertices.point(tree.sf_mesh.facets.vertex(f_id, 0));
            const auto &fv2 = tree.sf_mesh.vertices.point(tree.sf_mesh.facets.vertex(f_id, 1));
            const auto &fv3 = tree.sf_mesh.vertices.point(tree.sf_mesh.facets.vertex(f_id, 2));
            auto nf = GEO::cross((fv2 - fv1), (fv3 - fv1));
            Vector3 n, nt;
            n << nf[0], nf[1], nf[2];
            //
            auto &tv1 = mesh.tet_vertices[mesh.tets[t_id][(j + 1) % 4]].pos;
            auto &tv2 = mesh.tet_vertices[mesh.tets[t_id][(j + 2) % 4]].pos;
            auto &tv3 = mesh.tet_vertices[mesh.tets[t_id][(j + 3) % 4]].pos;
            if (Predicates::orient_3d(tv1, tv2, tv3, mesh.tet_vertices[mesh.tets[t_id][j]].pos)
                == Predicates::ORI_POSITIVE)
                nt = (tv2 - tv1).cross(tv3 - tv1);
            else
                nt = (tv3 - tv1).cross(tv2 - tv1);
            //
            if (nt.dot(n) > 0) {
                t.is_surface_fs[j] = 1;
                mesh.tets[opp_t_id].is_surface_fs[k] = -1;
            } else {
                t.is_surface_fs[j] = -1;
                mesh.tets[opp_t_id].is_surface_fs[k] = 1;
            }
        }
    }
}

#include <igl/winding_number.h>
void floatTetWild::boolean_operation(Mesh& mesh, int op){
    const int OP_UNION = 0;
    const int OP_INTERSECTION = 1;
    const int OP_DIFFERENCE = 2;

    Eigen::Matrix<Scalar, Eigen::Dynamic, 3> v1;
    Eigen::Matrix<int, Eigen::Dynamic, 3> f1;
    get_tracked_surface(mesh, v1, f1, 1);

    Eigen::Matrix<Scalar, Eigen::Dynamic, 3> v2;
    Eigen::Matrix<int, Eigen::Dynamic, 3> f2;
    get_tracked_surface(mesh, v2, f2, 2);

    Eigen::MatrixXd C(mesh.get_t_num(), 3);
    C.setZero();
    int index = 0;
    for (size_t i = 0; i < mesh.tets.size(); i++) {
        if (mesh.tets[i].is_removed)
            continue;
        for (int j = 0; j < 4; j++)
            C.row(index) += mesh.tet_vertices[mesh.tets[i][j]].pos.cast<double>();
        C.row(index) /= 4.0;
        index++;
    }

    Eigen::VectorXd w1, w2;
    igl::winding_number(Eigen::MatrixXd(v1.cast<double>()), Eigen::MatrixXi(f1), C, w1);
    igl::winding_number(Eigen::MatrixXd(v2.cast<double>()), Eigen::MatrixXi(f2), C, w2);

    int cnt = 0;
    for (int t_id = 0; t_id < mesh.tets.size(); ++t_id) {
        auto &t = mesh.tets[t_id];
        if(t.is_removed)
            continue;
        switch(op){
            case OP_UNION:
                if(w1(cnt)<=0.5 && w2(cnt)<=0.5)
                    t.is_removed = true;
                break;
            case OP_INTERSECTION:
                if(w1(cnt)<=0.5 || w2(cnt)<=0.5)
                    t.is_removed = true;
                break;
            case OP_DIFFERENCE:
                if(w1(cnt)<=0.5 || w2(cnt)>0.5)
                    t.is_removed = true;
                break;
        }
        cnt++;
    }

//    output_surface(mesh, "inner.stl");
}

void floatTetWild::filter_outside(Mesh& mesh, bool invert_faces) {
    Eigen::MatrixXd C(mesh.get_t_num(), 3);
    C.setZero();
    int index = 0;
    for (size_t i = 0; i < mesh.tets.size(); i++) {
        if (mesh.tets[i].is_removed)
            continue;
        for (int j = 0; j < 4; j++)
            C.row(index) += mesh.tet_vertices[mesh.tets[i][j]].pos.cast<double>();
        C.row(index) /= 4.0;
        index++;
    }
//    C.conservativeResize(index, 3);

    Eigen::Matrix<Scalar, Eigen::Dynamic, 3> V;
    Eigen::Matrix<int, Eigen::Dynamic, 3> F;
    get_tracked_surface(mesh, V, F);
    Eigen::VectorXd W;
    if (invert_faces) {
        Eigen::Matrix<int, Eigen::Dynamic, 1> tmp = F.col(1);
        F.col(1) = F.col(2).eval();
        F.col(2) = tmp;
    }
    igl::winding_number(Eigen::MatrixXd(V.cast<double>()), Eigen::MatrixXi(F), C, W);

    index = 0;
    int n_tets = 0;
    std::vector<bool> old_flags(mesh.tets.size());
    for (int t_id = 0; t_id < mesh.tets.size(); ++t_id) {
        auto &t = mesh.tets[t_id];
        old_flags[t_id] = t.is_removed;

        if (t.is_removed)
            continue;
        if (W(index) <= 0.5) {
            t.is_removed = true;
        } else
            n_tets++;
        index++;
    }

    if (n_tets <= 0) {
        if (invert_faces)
            logger().error("Empty mesh, problem with inverted faces");
        else {
            for (int t_id = 0; t_id < mesh.tets.size(); ++t_id) {
                auto &t = mesh.tets[t_id];
                t.is_removed = old_flags[t_id];
            }
            logger().debug("Empty mesh trying to reverse the faces");
            filter_outside(mesh, true);
        }
    }

    for (auto &v: mesh.tet_vertices) {
        if (v.is_removed)
            continue;
        bool is_remove = true;
        for (int t_id: v.conn_tets) {
            if (!mesh.tets[t_id].is_removed) {
                is_remove = false;
                break;
            }
        }
        v.is_removed = is_remove;
    }
}

void floatTetWild::mark_outside(Mesh& mesh, bool invert_faces){
    Eigen::MatrixXd C(mesh.get_t_num(), 3);
    C.setZero();
    int index = 0;
    for (size_t i = 0; i < mesh.tets.size(); i++) {
        if (mesh.tets[i].is_removed)
            continue;
        for (int j = 0; j < 4; j++)
            C.row(index) += mesh.tet_vertices[mesh.tets[i][j]].pos.cast<double>();
        C.row(index) /= 4.0;
        index++;
    }
//    C.conservativeResize(index, 3);

    Eigen::Matrix<Scalar, Eigen::Dynamic, 3> V;
    Eigen::Matrix<int, Eigen::Dynamic, 3> F;
    get_tracked_surface(mesh, V, F);
    if(invert_faces){
        Eigen::Matrix<int, Eigen::Dynamic, 1> tmp = F.col(1);
        F.col(1) = F.col(2).eval();
        F.col(2) = tmp;
    }
    Eigen::VectorXd W;
    igl::winding_number(Eigen::MatrixXd(V.cast<double>()), Eigen::MatrixXi(F), C, W);

    index = 0;
    int n_tets = 0;
    for (auto& t: mesh.tets) {
        if (t.is_removed){
            t.is_outside = true;
            continue;
        }
        if (W(index) <= 0.5)
            t.is_outside = true;
        else
            n_tets++;
        index++;
    }

    if(n_tets <= 0)
    {
        if(invert_faces)
            logger().error("Empty mesh, problem with inverted faces");
        else{
            for (auto& t: mesh.tets)
            {
                t.is_outside = false;
            }
            logger().debug("Empty mesh trying to reverse the faces");
            mark_outside(mesh, true);
        }
    }

    for (auto& v: mesh.tet_vertices) {
        if (v.is_removed){
            v.is_outside = true;
            continue;
        }
        bool is_outside = true;
        for(int t_id: v.conn_tets) {
            if (!mesh.tets[t_id].is_outside) {
                is_outside = false;
                break;
            }
        }
        v.is_outside = is_outside;
    }
}
