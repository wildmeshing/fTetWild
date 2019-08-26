#include <floattetwild/FloatTetCuttingCheck.h>
#include <floattetwild/LocalOperations.h>
#include <floattetwild/Predicates.hpp>
#include <floattetwild/MeshIO.hpp>
#include <floattetwild/intersections.h>

void floatTetWild::find_tets_for_cut_new(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
                                         int f_id, Mesh &mesh, std::vector<int> &intersection_results, std::vector<int> &oris){
    auto &tet_vertices = mesh.tet_vertices;
    auto &tets = mesh.tets;

    std::vector<bool> is_visited(tets.size(), false);

    const Vector3 &p1 = input_vertices[input_faces[f_id][0]];
    const Vector3 &p2 = input_vertices[input_faces[f_id][1]];
    const Vector3 &p3 = input_vertices[input_faces[f_id][2]];

    Vector3 tri_min, tri_max;
    get_bbox_face(p1, p2, p3, tri_min, tri_max);

    int test_f_id = 1374;
    if(f_id == test_f_id) {
        cout << std::setprecision(17);
        cout << p1.transpose() << endl;
        cout << p2.transpose() << endl;
        cout << p3.transpose() << endl;
    }

    std::queue<int> t_id_queue;
    for (int t_id:tet_vertices[input_faces[f_id][0]].conn_tets) {
        t_id_queue.push(t_id);
        is_visited[t_id] = true;
    }
    while (!t_id_queue.empty()) {
        int t_id = t_id_queue.front();
        t_id_queue.pop();

        if(f_id == test_f_id && t_id == 9421){
            cout<<"here 0"<<endl;
        }

//        if (is_visited[t_id])
//            continue;
//        is_visited[t_id] = true;

        //check if the tet is on the same side of the plane
        int cnt_pos = 0;
        int cnt_neg = 0;
        int cnt_zero = 0;
        for (int j = 0; j < 4; j++) {
            if (oris[tets[t_id][j]] == Predicates::ORI_UNKNOWN)
                oris[tets[t_id][j]] = Predicates::orient_3d(p1, p2, p3, tet_vertices[tets[t_id][j]].pos);
            if (oris[tets[t_id][j]] == Predicates::ORI_POSITIVE)
                cnt_pos++;
            else if (oris[tets[t_id][j]] == Predicates::ORI_NEGATIVE)
                cnt_neg++;
            else
                cnt_zero++;
        }
        if ((cnt_neg == 0 || cnt_pos == 0) && cnt_zero <= 2)
            continue;

//        if(f_id == test_f_id && t_id == 9421){
//            cout<<"here 1"<<endl;
//        }

        //check triangle intersection
        for (int j = 0; j < 4; j++) {
//            int opp_t_id = tets[t_id].opp_t_ids[j];
//            if (opp_t_id >= 0 && is_visited[opp_t_id])
//                continue;

            const auto &tp1 = tet_vertices[tets[t_id][mod4(j + 1)]].pos;
            const auto &tp2 = tet_vertices[tets[t_id][mod4(j + 2)]].pos;
            const auto &tp3 = tet_vertices[tets[t_id][mod4(j + 3)]].pos;

//            if(f_id == test_f_id && t_id == 9421){
//                cout<<"here 2 "<<j<<endl;
//            }

            //check bbox
            Vector3 tf_min, tf_max;
            get_bbox_face(tp1, tp2, tp3, tf_min, tf_max);

//            if(f_id == test_f_id && t_id == 9421){
//                cout<<"here 2 "<<j<<endl;
//                cout<<"tri_min "<<tri_min.transpose()<<endl;
//                cout<<"tri_max "<<tri_max.transpose()<<endl;
//                cout<<"tf_min "<<tf_min.transpose()<<endl;
//                cout<<"tf_max "<<tf_max.transpose()<<endl;
//                cout << tet_vertices[tets[t_id][(j + 1) % 4]].pos.transpose() << endl;
//                cout << tet_vertices[tets[t_id][(j + 2) % 4]].pos.transpose() << endl;
//                cout << tet_vertices[tets[t_id][(j + 3) % 4]].pos.transpose() << endl;
//            }

            if (!is_bbox_intersected(tri_min, tri_max, tf_min, tf_max))
                continue;

            std::vector<int> new_t_ids;
            for (int k = 0; k < 3; k++) {
                for (int n_t_id: tet_vertices[tets[t_id][mod4(j + 1 + k)]].conn_tets)
                    if (!is_visited[n_t_id]) {
                        t_id_queue.push(n_t_id);
                        is_visited[n_t_id] = true;
                    }
            }

//            if(f_id == test_f_id && t_id == 9421){
//                cout<<"here 3 "<<j<<endl;
//            }

            //check oris
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
            if ((cnt_neg_f == 0 || cnt_pos_f == 0) && cnt_zero_f <= 2)
                continue;

//            if(f_id == test_f_id && t_id == 9421){
//                cout<<"here 4 "<<j<<endl;
//            }

            cnt_pos_f = 0;
            cnt_neg_f = 0;
            cnt_zero_f = 0;
            for (int k = 0; k < 3; k++) {
                int ori = Predicates::orient_3d(tp1, tp2, tp3, input_vertices[input_faces[f_id][k]]);
                if (ori == Predicates::ORI_POSITIVE)
                    cnt_pos_f++;
                else if (ori == Predicates::ORI_NEGATIVE)
                    cnt_neg_f++;
                else
                    cnt_zero_f++;
            }
            if ((cnt_neg_f == 0 || cnt_pos_f == 0) && cnt_zero_f <= 2)
                continue;

//            if(f_id == test_f_id && t_id == 9421){
//                cout<<"here 5 "<<j<<endl;
//            }

//            std::vector<int> new_t_ids;
//            for (int k = 0; k < 3; k++) {
//                for (int n_t_id: tet_vertices[tets[t_id][mod4(j + 1 + k)]].conn_tets)
//                    if (!is_visited[n_t_id]) {
//                        t_id_queue.push(n_t_id);
//                        is_visited[n_t_id] = true;
//                    }
//            }

            intersection_results.push_back(t_id);
        }
    }
    vector_unique(intersection_results);

    if(f_id == test_f_id) {
        std::vector<int> check_t_ids = intersection_results;
        for (int t_id:check_t_ids) {
            for (int j = 0; j < 4; j++) {
                cout << tet_vertices[tets[t_id][(j + 1) % 4]].pos.transpose() << endl;
                cout << tet_vertices[tets[t_id][(j + 2) % 4]].pos.transpose() << endl;
                cout << tet_vertices[tets[t_id][(j + 3) % 4]].pos.transpose() << endl;
            }
        }
        for (int t_id:check_t_ids) {
            for (int j = 0; j < 4; j++) {
                cout << oris[tets[t_id][(j + 1) % 4]] << " ";
                cout << oris[tets[t_id][(j + 2) % 4]] << " ";
                cout << oris[tets[t_id][(j + 3) % 4]] << " ";
            }
            cout<<endl;
        }
        //pausee();
    }
}

bool floatTetWild::check(const Mesh &mesh) {
    /*
    auto &tet_vertices = mesh.tet_vertices;
    auto &tets = mesh.tets;

    if(mesh.params.log_level > 0)
        return true;

    //conn_tets
    for (int i = 0; i < tets.size(); i++) {
        if (tets[i].is_removed)
            continue;
        for (int j = 0; j < 4; j++) {
            if (tet_vertices[tets[i][j]].conn_tets.find(i) == tet_vertices[tets[i][j]].conn_tets.end()) {
                cout << "conn_tets error!" << endl;
                //pausee();
            }
        }
    }
    for (int i = 0; i < tet_vertices.size(); i++) {
        if (tet_vertices[i].is_removed)
            continue;
        for(int t_id:tet_vertices[i].conn_tets){
            if(t_id>tets.size() || t_id<0){
                cout<<"t_id>tets.size() || t_id<0"<<endl;
                //pausee();
            }
            int j = tets[t_id].find(i);
            if(j<0){
                cout<<"conn_tets error: j<0"<<endl;
                //pausee();
            }
        }
    }

    bool is_valid = true;

//    MeshIO::write_mesh("/Users/yixinhu/Downloads/test.mesh", mesh);
    int cnt_t = 0;
    for(auto& t:tets)
        if(!t.is_removed)
            cnt_t++;
    int cnt_v = 0;
    for(auto& v:tet_vertices)
        if(!v.is_removed)
            cnt_v++;
    cout << "#v = " << cnt_v << endl;
    cout << "#t = " << cnt_t << endl;

    //check inversion
    int cnt = 0;
    for (int i = 0; i < tets.size(); i++) {
        if(tets[i].is_removed)
            continue;

        if (is_inverted(mesh, i)) {
            for (int j = 0; j < 4; j++)
                cout << tets[i][j] << " ";
            for (int j = 0; j < 4; j++)
                cout << tet_vertices[tets[i][j]].pos << endl;
            cout << endl;
            cout << Predicates::orient_3d(tet_vertices[tets[i][0]].pos, tet_vertices[tets[i][1]].pos,
                                          tet_vertices[tets[i][2]].pos, tet_vertices[tets[i][3]].pos) << endl;
            cnt++;
        }
        Scalar area = Predicates::orient_3d_volume(tet_vertices[tets[i][0]].pos, tet_vertices[tets[i][1]].pos,
                                                   tet_vertices[tets[i][2]].pos, tet_vertices[tets[i][3]].pos);
        if (fabs(area) < SCALAR_ZERO) {
//            cout << fabs(area) << endl;
//            //pausee();
        }
    }
    if (cnt > 0) {
        is_valid = false;

        cout << "inversion error " << cnt << endl;
        //pausee();
    }

    //check Euler
    auto find_euler_error = [&](const std::vector<std::array<int, 2>>& all_edges) {
        for (int i = 0; i < tets.size(); i++) {
            if(tets[i].is_removed)
                continue;
//        for (int i = 0; i < all_edges.size(); i++) {
            std::vector<int> one_ring = {i};
            for (int j = 0; j < 4; j++) {
//            for (int j = 0; j < 2; j++) {
                for (auto &ii: tet_vertices[tets[i][j]].conn_tets)
//                for (auto &ii: tet_vertices[all_edges[i][j]].conn_tets)
                    one_ring.push_back(ii);
            }
            vector_unique(one_ring);

            //check
            std::vector<std::array<int, 2>> edges;
            get_all_edges(mesh, one_ring, edges);
            std::vector<std::array<int, 3>> faces;
            std::vector<int> vertices;
            for (auto &t_id: one_ring) {
                auto &t = tets[t_id];
                for (int j = 0; j < 4; j++) {
                    vertices.push_back(t[j]);
                    std::array<int, 3> f = {{t[j], t[mod4(j + 1)], t[mod4(j + 2)]}};
                    std::sort(f.begin(), f.end());
                    faces.push_back(f);
                }
            }
            vector_unique(faces);
            vector_unique(vertices);

            int euler = vertices.size() - edges.size() + faces.size() - one_ring.size();
            if (euler != 1) {
                cout << "find" << endl;
                cout << "euler = " << euler << endl;
                for (auto &t_id:one_ring) {
                    cout << t_id << ": " << tets[t_id][0] << " " << tets[t_id][1] << " " << tets[t_id][2] << " "
                         << tets[t_id][3] << endl;
                }
//                cout << "error t: " << i << "; " << tets[i][0] << " " << tets[i][1] << " " << tets[i][2] << " "
//                     << tets[i][3] << endl;
                cout << "error edge: " << all_edges[i][0] << " " << all_edges[i][1] << endl;
                MeshIO::write_mesh("/Users/yixinhu/Downloads/60246_euler.stl_euler.mesh", mesh, one_ring);
                MeshIO::write_mesh("/Users/yixinhu/Downloads/60246_euler.stl_euler.msh", mesh, one_ring);
                //pausee();
            }
        }
    };

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
        is_valid = false;

        cout << "euler error " << euler << endl;
//        exit(0);

//        for(auto& t: tets){
//            for(int j=0;j<4;j++){
//                if(t.opp_t_ids[j] == -1){
//                    cout<<tet_vertices[t[(j+1)%4]].pos[0]<<" "<<tet_vertices[t[(j+1)%4]].pos[1]<<" "<<tet_vertices[t[(j+1)%4]].pos[2]<<endl;
//                    cout<<tet_vertices[t[(j+2)%4]].pos[0]<<" "<<tet_vertices[t[(j+2)%4]].pos[1]<<" "<<tet_vertices[t[(j+2)%4]].pos[2]<<endl;
//                    cout<<tet_vertices[t[(j+3)%4]].pos[0]<<" "<<tet_vertices[t[(j+3)%4]].pos[1]<<" "<<tet_vertices[t[(j+3)%4]].pos[2]<<endl;
//                    cout<<endl;
//                }
//            }
//        }

//        find_euler_error(edges);

        //pausee();
    }

    //check conn_tets
    std::vector<std::unordered_set<int>> tmp_conn_tets(tet_vertices.size());
    for (int i = 0; i < tets.size(); i++) {
        if(tets[i].is_removed)
            continue;

        for (int j = 0; j < 4; j++)
            tmp_conn_tets[tets[i][j]].insert(i);
    }
    for (int i = 0; i < tet_vertices.size(); i++) {
        if(tet_vertices[i].is_removed)
            continue;

        if (tet_vertices[i].conn_tets != tmp_conn_tets[i]) {
            is_valid = false;

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

    //check opp_t_ids
//    std::vector<std::array<int, 4>> tmp_opp_t_ids(tets.size());
//    for (int i = 0; i < tets.size(); i++) {
//        if(tets[i].is_removed)
//            continue;
//        for (int j = 0; j < 4; j++) {
//            tmp_opp_t_ids[i][j] = -1;
//            std::unordered_set<int> tmp;
//            set_intersection(mesh.tet_vertices[tets[i][mod4(j + 1)]].conn_tets,
//                             mesh.tet_vertices[tets[i][mod4(j + 2)]].conn_tets, tmp);
//            std::vector<int> pair;
//            set_intersection(mesh.tet_vertices[tets[i][mod4(j + 3)]].conn_tets, tmp, pair);
//            if (pair.size() == 2) {
//                int opp_t_id = pair[0] == i ? pair[1] : pair[0];
//                tmp_opp_t_ids[i][j] = opp_t_id;
//            }
//
//        }
//        if (tmp_opp_t_ids[i] != tets[i].opp_t_ids) {
//            is_valid = false;
//
//            cout << "opp_t_ids error" << endl;
//            for (auto &ii:tmp_opp_t_ids[i])
//                cout << ii << " ";
//            cout << endl;
//            for (auto &ii:tets[i].opp_t_ids)
//                cout << ii << " ";
//            cout << endl;
//            //pausee();
//        }
//    }
//
//    return true;

    //check manifoldness
    //faces - done above
    for(auto& f: faces){
        std::unordered_set<int> tmp;
        set_intersection(mesh.tet_vertices[f[0]].conn_tets,
                         mesh.tet_vertices[f[1]].conn_tets, tmp);
        std::vector<int> pair;
        set_intersection(mesh.tet_vertices[f[2]].conn_tets, tmp, pair);
        if (pair.size() > 2) {
            cout << "face nonmanifold" << endl;
            //pausee();
        }
    }
    //edges
    for(auto& e:edges){
        std::vector<std::array<int, 3>> edge_faces;
        std::unordered_set<int> tmp;
        set_intersection(mesh.tet_vertices[e[0]].conn_tets,
                         mesh.tet_vertices[e[1]].conn_tets, tmp);
        std::vector<int> tmp_v_ids;
        for(int t_id: tmp){
            for(int j=0;j<4;j++) {
                if(tets[t_id][j]!=e[0] && tets[t_id][j]!=e[1]){
                    tmp_v_ids.push_back(tets[t_id][j]);
                }else{
                    std::array<int, 3> f = {{tets[t_id][mod4(j+1)], tets[t_id][mod4(j+2)], tets[t_id][mod4(j+3)]}};
                    std::sort(f.begin(), f.end());
                    edge_faces.push_back(f);
                }
            }
        }
        vector_unique(tmp_v_ids);
        if(tmp_v_ids.size() - tmp.size() > 1){
            cout << "edge nonmanifold" << endl;
            //pausee();
        }

//        vector_unique(edge_faces);
//        cout<<"edge_faces"<<endl;
//        for(auto& f:edge_faces){
//            for(int j=0;j<3;j++)
//                cout<<f[j]<<" ";
//            cout<<endl;
//        }
//        for(int i=0;i<tet_vertices.size();i++){
//            cout<<tet_vertices[i].pos[0]<<" "<<tet_vertices[i].pos[1]<<" "<<tet_vertices[i].pos[2]<<endl;
//        }

        std::vector<int> cnt_f_v(tet_vertices.size(), 0);
        for(auto& f: edge_faces){
            for(int j=0;j<3;j++) {
                if(f[j]!=e[0] && f[j]!=e[1])
                    cnt_f_v[f[j]]++;
            }
        }
        for(int i=0;i<cnt_f_v.size();i++)
            if(cnt_f_v[i]>4){
                cout << "edge nonmanifold 2" << endl;
                cout<<"e "<<e[0]<<" "<<e[1]<<endl;
                cout<<"v "<<i<<endl;
                for(int t_id:tmp)
                    cout<<tets[t_id][0]<<' '<<tets[t_id][1]<<' '<<tets[t_id][2]<<' '<<tets[t_id][3]<<endl;
                //pausee();
            }
    }

    //vertices
    for(int i = 0; i < mesh.tet_vertices.size(); ++i) {
        if(tet_vertices[i].is_removed)
            continue;
        const auto &vertex = mesh.tet_vertices[i];
        const auto &neighs = vertex.conn_tets;
        std::vector<Vector3i> surf_faces;
        std::unordered_set<int> vertices;
        std::set<std::array<int, 2>> edges;
        surf_faces.reserve(neighs.size());

        for(const int n : neighs) {
            int index = 0;
            surf_faces.emplace_back();
            auto &f = surf_faces.back();
            for(int j = 0; j < 4; ++j) {
                const int vid = mesh.tets[n][j];
                if(vid != i){
                    f[index] = vid;
                    vertices.insert(vid);
                    ++index;
                }
            }

            assert(index == 3);

            for(int j = 0; j < 3; ++j) {
                const int jp = mod3(j+1);
                edges.insert({{{std::min(f[j], f[jp]), std::max(f[j], f[jp])}}});
            }
        }

        int e = surf_faces.size() - edges.size() + vertices.size();
        if(!(e == 2 || e == 1)) {
            cout << " v nonmanifold "<<i << endl;
            for (int i = 0; i < tet_vertices.size(); i++) {
                cout << tet_vertices[i].pos[0] << " " << tet_vertices[i].pos[1] << " " << tet_vertices[i].pos[2]
                     << endl;
            }
            for (auto &f: surf_faces) {
                cout << f[0] << " " << f[1] << " " << f[2] << endl;
            }

            cout<<"tet_vertices[i].conn_tets.size() = "<<tet_vertices[i].conn_tets.size()<<endl;

            //pausee();
        }
    }

    //check centroid
//    for(int i=0;i<tets.size();i++){
//        Vector3 centroid = (tet_vertices[tets[i][0]].pos+tet_vertices[tets[i][1]].pos+tet_vertices[tets[i][2]].pos+tet_vertices[tets[i][3]].pos)/4;
//        GEO::MeshCellsAABB
//    }

    return is_valid;
     */
    return true;
}

void floatTetWild::check_is_surface_fs(const Mesh &mesh){
    if(mesh.params.log_level > 0)
        return;

    for (auto &t: mesh.tets) {
        if(t.is_removed)
            continue;

        for (int j = 0; j < 4; j++) {
            int opp_t_id = t.opp_t_ids[j];
            if (opp_t_id < 0)
                continue;
            int my_sf = t.is_surface_fs[j];
            int opp_sf;
            for (int k = 0; k < 4; k++) {
                if (mesh.tets[opp_t_id][k] != t[mod4(j + 1)] && mesh.tets[opp_t_id][k] != t[mod4(j + 2)] &&
                    mesh.tets[opp_t_id][k] != t[mod4(j + 3)]) {
                    opp_sf = mesh.tets[opp_t_id].is_surface_fs[k];
                    break;
                }
            }
            if (
                    my_sf == 0 || opp_sf == 0 ||
                    (my_sf == NOT_SURFACE && opp_sf != NOT_SURFACE) || (my_sf != NOT_SURFACE && my_sf + opp_sf != 0)) {
                cout << "surface error (NOT_SURFACE)" << endl;
                cout << my_sf << " " << opp_sf << endl;
                //pausee();
            }

//            if (my_sf == NOT_SURFACE) {
//                if (opp_sf != NOT_SURFACE) {
//                    cout << "surface error (NOT_SURFACE)" << endl;
//                    cout << my_sf << " " << opp_sf << endl;
//                    //pausee();
//                }
//            } else {
//                if (my_sf + opp_sf != 0) {
//                    cout << "surface error " << endl;
//                    cout << my_sf << " " << opp_sf << endl;
//                    //pausee();
//                }
//            }
        }
    }
}

void floatTetWild::check_cut_f_ids(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
                                   const Mesh &mesh, std::vector<std::array<std::vector<int>, 4>>& cut_f_ids) {
    const auto &tet_vertices = mesh.tet_vertices;
    const auto &tets = mesh.tets;

//    return;
//    for (int t_id = 0; t_id < cut_f_ids.size(); t_id++) {
//        if (tets[t_id].is_removed)
//            continue;
//
//        for (int j = 0; j < 4; j++) {
//            if (cut_f_ids[t_id][j].empty())
//                continue;
//            if (tets[t_id].is_surface_fs[j] != NOT_SURFACE)
//                continue;
//
//            //todo: avoid duplicated check on neighbors, since you know opp_t_ids info here
//
//            const Vector3 &v0 = tet_vertices[tets[t_id][mod4(j + 1)]].pos;
//            const Vector3 &v1 = tet_vertices[tets[t_id][mod4(j + 2)]].pos;
//            const Vector3 &v2 = tet_vertices[tets[t_id][mod4(j + 3)]].pos;
//
//            bool is_found = false;
//
//            int t = get_t(v0, v1, v2);
//            std::array<Vector2, 3> vs_tet = {{to_2d(v0, t), to_2d(v1, t), to_2d(v2, t)}};
//            for (int f_id:cut_f_ids[t_id][j]) {
//                std::array<int, 3> f_in = {{input_faces[f_id][0], input_faces[f_id][1], input_faces[f_id][2]}};
//                std::sort(f_in.begin(), f_in.end());
//                std::array<int, 3> f_tet = {{tets[t_id][mod4(j + 1)], tets[t_id][mod4(j + 2)], tets[t_id][mod4(
//                        j + 3)]}};
//                std::sort(f_tet.begin(), f_tet.end());
//                if (f_in == f_tet) {
//                    is_found = true;
//                    break;
//                }
//
//                //check intersection
//                std::array<Vector2, 3> vs_tri = {{to_2d(input_vertices[input_faces[f_id][0]], t),
//                                                         to_2d(input_vertices[input_faces[f_id][1]], t),
//                                                         to_2d(input_vertices[input_faces[f_id][2]], t)}};
//
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
//                if (!(cnt_pos == 0 || cnt_neg == 0))
//                    continue;
//
//                is_found = true;
//                break;
//            }
//
//            if (!is_found) {
//                if(cut_f_ids[t_id][j].size()<2)
//                    continue;
//
//                cout << tet_vertices[tets[t_id][mod4(j + 1)]].pos.transpose() << endl;
//                cout << tet_vertices[tets[t_id][mod4(j + 2)]].pos.transpose() << endl;
//                cout << tet_vertices[tets[t_id][mod4(j + 3)]].pos.transpose() << endl;
//
//                for (int f_id:cut_f_ids[t_id][j]) {
//                    for (int k = 0; k < 3; k++)
//                        cout << input_vertices[input_faces[f_id][k]].transpose() << endl;
//                }
//
//                for (int f_id:cut_f_ids[t_id][j])
//                    cout << f_id << " ";
//                cout << endl;
//                //pausee();
//            }
//        }
//    }
//
//    return;

    if (mesh.params.log_level > 0)
        return;

//    const auto &tet_vertices = mesh.tet_vertices;
//    const auto &tets = mesh.tets;

    for (auto &cut_f_id: cut_f_ids) {
        for (int j = 0; j < 4; j++)
            std::sort(cut_f_id[j].begin(), cut_f_id[j].end());
    }
    for (int i = 0; i < tets.size(); i++) {
        const auto &t = mesh.tets[i];
        if (t.is_removed)
            continue;
        for (int j = 0; j < 4; j++) {
            int opp_t_id = t.opp_t_ids[j];
            if (opp_t_id < 0)
                continue;
            auto &my_sf = cut_f_ids[i][j];
            for (int k = 0; k < 4; k++) {
                if (mesh.tets[opp_t_id][k] != t[mod4(j + 1)] && mesh.tets[opp_t_id][k] != t[mod4(j + 2)] &&
                    mesh.tets[opp_t_id][k] != t[mod4(j + 3)]) {
                    if (my_sf != cut_f_ids[opp_t_id][k]) {
                        cout << "cut_f_ids error" << endl;
                        cout << "cut_f_ids[t_id]: ";
                        for (auto &ii: my_sf)
                            cout << ii << " ";
                        cout << endl;
                        cout << "cut_f_ids[opp_t_id]: ";
                        for (auto &ii: cut_f_ids[opp_t_id][k])
                            cout << ii << " ";
                        cout << endl;
                        cout << "t_id " << i << " " << j << ": " << tets[i][0] << " " << tets[i][1] << " " << tets[i][2]
                             << " "
                             << tets[i][3] << endl;
                        cout << "opp_t_id " << opp_t_id << " " << k << ": " << tets[opp_t_id][0] << " "
                             << tets[opp_t_id][1] << " "
                             << tets[opp_t_id][2] << " " << tets[opp_t_id][3] << endl;
                        //pausee();
                    }
                    break;
                }
            }
        }
    }

    return;
//check coplanar
    for (int i = 0; i < tets.size(); i++) {
        const auto &t = mesh.tets[i];
        if (t.is_removed)
            continue;
        for (int j = 0; j < 4; j++) {
            Vector3 n1 = (mesh.tet_vertices[t[(j + 1) % 4]].pos - mesh.tet_vertices[t[(j + 3) % 4]].pos).cross(
                    mesh.tet_vertices[t[(j + 2) % 4]].pos - mesh.tet_vertices[t[(j + 3) % 4]].pos);
            n1.normalize();
            for (int f_id: cut_f_ids[i][j]) {
                Vector3 n2 = (input_vertices[input_faces[f_id][0]] - input_vertices[input_faces[f_id][2]]).cross(
                        input_vertices[input_faces[f_id][1]] - input_vertices[input_faces[f_id][2]]);
                n2.normalize();
//                Vector3 v1 = n1 - n2;
//                Vector3 v2 = n1 + n2;
//                if ((abs(v1[0]) < SCALAR_ZERO && abs(v1[1]) < SCALAR_ZERO && abs(v1[2]) < SCALAR_ZERO)
//                    || (abs(v2[0]) < SCALAR_ZERO && abs(v2[1]) < SCALAR_ZERO && abs(v2[2]) < SCALAR_ZERO))
//                    continue;
                if (abs(n1.dot(n2)) < 0.1) {
                    cout<<"abs(n1.dot(n2)) < 0.1"<<endl;
                    cout << n1.transpose() << endl;
                    cout << n2.transpose() << endl;
                    cout << "tet " << i << endl;
                    cout << "tri " << f_id << endl;
                    //pausee();
                }
            }
        }
    }
}

void floatTetWild::plot_cover_for_tetf(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
        const Mesh &mesh, const std::vector<std::array<std::vector<int>, 4>>& cut_f_ids,
        int t_id, int j) {
    const auto &tet_vertices = mesh.tet_vertices;
    const auto &tets = mesh.tets;

    const auto &v0 = tet_vertices[tets[t_id][mod4(j + 1)]].pos;
    const auto &v1 = tet_vertices[tets[t_id][mod4(j + 2)]].pos;
    const auto &v2 = tet_vertices[tets[t_id][mod4(j + 3)]].pos;
    int t = get_t(v0, v1, v2);

    cout << "plot_cover_for_tetf: " << t_id << " " << j << endl;
    cout<<tets[t_id][0]<<" "<<tets[t_id][1]<<" "<<tets[t_id][2]<<" "<<tets[t_id][3]<<endl;
    cout << "/////////3d" << endl;
    cout << v0.transpose() << endl;
    cout << v1.transpose() << endl;
    cout << v2.transpose() << endl;
    for (int f_id:cut_f_ids[t_id][j]) {
        cout << input_vertices[input_faces[f_id][0]].transpose() << endl;
        cout << input_vertices[input_faces[f_id][1]].transpose() << endl;
        cout << input_vertices[input_faces[f_id][2]].transpose() << endl;
    }

    cout << "/////////2d" << endl;
    cout << to_2d(v0, t).transpose() << endl;
    cout << to_2d(v1, t).transpose() << endl;
    cout << to_2d(v2, t).transpose() << endl;
    for (int f_id:cut_f_ids[t_id][j]) {
        cout << to_2d(input_vertices[input_faces[f_id][0]], t).transpose() << endl;
        cout << to_2d(input_vertices[input_faces[f_id][1]], t).transpose() << endl;
        cout << to_2d(input_vertices[input_faces[f_id][2]], t).transpose() << endl;
    }

    //pausee();
}

void floatTetWild::plot_cover_for_trif(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
        const Mesh &mesh, const std::vector<std::array<std::vector<int>, 4>>& cut_f_ids,
        int f_id) {
    const auto &tet_vertices = mesh.tet_vertices;
    const auto &tets = mesh.tets;

    cout << "plot_cover_for_trif " << f_id << endl;
    std::vector<std::array<int, 2>> cut_t_ids;
    for (int t_id = 0; t_id < cut_f_ids.size(); t_id++) {
        if (mesh.tets[t_id].is_removed)
            continue;
        for (int j = 0; j < 4; j++) {
            for (int i:cut_f_ids[t_id][j]) {
                if (i == f_id)
                    cut_t_ids.push_back({{t_id, j}});
            }
        }
    }
    cout << "cut_t_ids.size() = " << cut_t_ids.size() << endl;

    cout << std::setprecision(17);
    cout << input_vertices[input_faces[f_id][0]].transpose() << endl;
    cout << input_vertices[input_faces[f_id][1]].transpose() << endl;
    cout << input_vertices[input_faces[f_id][2]].transpose() << endl;

    for (auto &t_info:cut_t_ids) {
        int t_id = t_info[0];
        int j = t_info[1];
        cout << tet_vertices[tets[t_id][mod4(j + 1)]].pos.transpose() << endl;
        cout << tet_vertices[tets[t_id][mod4(j + 2)]].pos.transpose() << endl;
        cout << tet_vertices[tets[t_id][mod4(j + 3)]].pos.transpose() << endl;
    }
    //pausee();
}

void floatTetWild::find_bad_cover_for_tetf(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
        const Mesh &mesh, const std::vector<std::array<std::vector<int>, 4>>& cut_f_ids, const AABBWrapper& tree, int t_id, int j,
        bool is_on_surface){
    const auto &tet_vertices = mesh.tet_vertices;
    const auto &tets = mesh.tets;

    const auto &v0 = tet_vertices[tets[t_id][mod4(j + 1)]].pos;
    const auto &v1 = tet_vertices[tets[t_id][mod4(j + 2)]].pos;
    const auto &v2 = tet_vertices[tets[t_id][mod4(j + 3)]].pos;

    if(!is_on_surface) {
//        if (!tree.is_out_sf_envelope(ps, mesh.params.eps_2)){
//            cout << "envelope_check says this tet face should be on the surface" << endl;
//            plot_cover_for_tetf(COMMON_INPUT_FOR_CHECK, tree, t_id, j);
//        }
        int t = get_t(v0, v1, v2);

        std::vector<Vector2> cs;
        get_samples({{to_2d(v0, t), to_2d(v1, t), to_2d(v2, t)}}, cs);

        bool is_inside = false;
        for(auto& c: cs) {
            is_inside = false;
            for (int f_id:cut_f_ids[t_id][j]) {
                std::array<Vector2, 3> vs_tri = {{to_2d(input_vertices[input_faces[f_id][0]], t),
                                                         to_2d(input_vertices[input_faces[f_id][1]], t),
                                                         to_2d(input_vertices[input_faces[f_id][2]], t)}};

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
                cout<<"c = "<<c.transpose()<<endl;
                is_inside = true;
                break;
            }
            if(is_inside)
                break;
        }
        if(is_inside) {
            cout << "find a tetf PARTIALLY COVERED by certain trif" << endl;
            plot_cover_for_tetf(COMMON_INPUT_FOR_CHECK, t_id, j);
        }
    } else {
        std::vector<GEO::vec3> ps;
        sample_triangle({{v0, v1, v2}}, ps, mesh.params.dd);
        Scalar dist = tree.dist_sf_envelope(ps, mesh.params.eps_2);
        if (dist > mesh.params.eps_2) {
            cout << "envelope_check says this tetf SHOULD NOT be on the surface" << endl;
            cout << "dist = " << dist << endl;
            cout << "cut_f_ids = ";
            for (int f_id:cut_f_ids[t_id][j])
                cout << f_id << " ";
            cout << endl;
            plot_cover_for_tetf(COMMON_INPUT_FOR_CHECK, t_id, j);
        }
    }
}

void floatTetWild::find_bad_cover_for_trif(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces,
        const Mesh &mesh, const std::vector<std::array<std::vector<int>, 4>>& cut_f_ids, int f_id) {
    const auto &tet_vertices = mesh.tet_vertices;
    const auto &tets = mesh.tets;

    cout << "find_bad_cover_for_trif..." << endl;
    std::vector<std::vector<std::array<int, 2>>> cut_t_ids(input_faces.size());
    for (int t_id = 0; t_id < cut_f_ids.size(); t_id++) {
        if (mesh.tets[t_id].is_removed)
            continue;
        for (int j = 0; j < 4; j++) {
            for (int i:cut_f_ids[t_id][j]) {
                if(f_id < 0)
                    cut_t_ids[i].push_back({{t_id, j}});
                else{
                    if(i == f_id)
                        cut_t_ids[i].push_back({{t_id, j}});
                }
            }
        }
    }

    for (int ff_id = 0; ff_id < cut_t_ids.size(); ff_id++) {
        auto &infos = cut_t_ids[ff_id];
        if (infos.empty())
            continue;
        if (infos.size() == 1) {
            int t_id = infos[0][0];
            int j = infos[0][1];
            std::array<int, 3> tf = {{tets[t_id][mod4(j + 1)], tets[t_id][mod4(j + 2)], tets[t_id][mod4(j + 3)]}};
            std::array<int, 3> rf = {{input_faces[ff_id][0], input_faces[ff_id][1], input_faces[ff_id][2]}};
            std::sort(tf.begin(), tf.end());
            std::sort(rf.begin(), rf.end());
            if (tf == rf)
                continue;
        }

        int t = get_t(input_vertices[input_faces[ff_id][0]], input_vertices[input_faces[ff_id][1]],
                      input_vertices[input_faces[ff_id][2]]);
        std::array<Vector2, 3> rf = {{to_2d(input_vertices[input_faces[ff_id][0]], t),
                                             to_2d(input_vertices[input_faces[ff_id][1]], t),
                                             to_2d(input_vertices[input_faces[ff_id][2]], t),}};

        std::vector<Vector2> cs;
        get_samples(rf, cs);
        for (auto &c:cs) {
            bool is_found = false;
            for (auto &t_info:infos) {
                int t_id = t_info[0];
                int j = t_info[1];
                std::array<Vector2, 3> tf = {{to_2d(tet_vertices[tets[t_id][mod4(j + 1)]].pos, t),
                                                     to_2d(tet_vertices[tets[t_id][mod4(j + 2)]].pos, t),
                                                     to_2d(tet_vertices[tets[t_id][mod4(j + 3)]].pos, t)}};
                int cnt_pos = 0;
                int cnt_neg = 0;
                for (int i = 0; i < 3; i++) {
                    int ori = Predicates::orient_2d(tf[i], tf[mod3(i + 1)], c);
                    if (ori == Predicates::ORI_POSITIVE)
                        cnt_pos++;
                    else if (ori == Predicates::ORI_NEGATIVE)
                        cnt_neg++;
                }
                if (!(cnt_pos == 0 || cnt_neg == 0))
                    continue;

                is_found = true;
                break;
            }
            if (!is_found) {
                cout << "find ERROR" << endl;
                cout << "ff_id = " << ff_id << endl;
                cout << "cut_t_ids[ff_id].size() = " << cut_t_ids[ff_id].size() << endl;

                cout << std::setprecision(17);
                cout << input_vertices[input_faces[ff_id][0]].transpose() << endl;
                cout << input_vertices[input_faces[ff_id][1]].transpose() << endl;
                cout << input_vertices[input_faces[ff_id][2]].transpose() << endl;
                for (auto &t_info:infos) {
                    int t_id = t_info[0];
                    int j = t_info[1];
                    cout << tet_vertices[tets[t_id][mod4(j + 1)]].pos.transpose() << endl;
                    cout << tet_vertices[tets[t_id][mod4(j + 2)]].pos.transpose() << endl;
                    cout << tet_vertices[tets[t_id][mod4(j + 3)]].pos.transpose() << endl;
                }
                //pausee();

                break;
            }
        }
    }
}

void floatTetWild::get_samples(const std::array<Vector2, 3>& rf, std::vector<Vector2>& cs, Scalar ratio){
    Scalar a = (rf[0] - rf[1]).norm();
    Scalar b = (rf[2] - rf[1]).norm();
    Scalar c = (rf[2] - rf[0]).norm();
    Scalar s = (a+b+c)/2;
    Scalar area = sqrt(s*(s-a)*(s-b)*(s-c));
    area = sqrt(area);
    Scalar increment = ratio / area;
    for (Scalar r1 = increment; r1 < 1; r1 += increment) {
        for (Scalar r2 = increment; r2 < 1; r2 += increment  ) {
            Scalar sqrtR = sqrt(r1);
            Scalar A = (1 - sqrtR);
            Scalar B = (sqrtR * (1 - r2));
            Scalar C = (sqrtR * r2);
            cs.emplace_back(A * rf[0][0] + B * rf[1][0] + C * rf[2][0], A * rf[0][1] + B * rf[1][1] + C * rf[2][1]);
        }
    }
}
