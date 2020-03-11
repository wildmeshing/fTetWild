// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include <floattetwild/intersections.h>

#include <floattetwild/Predicates.hpp>
#include <floattetwild/LocalOperations.h>

bool floatTetWild::seg_line_intersection_2d(const std::array<Vector2, 2> &seg, const std::array<Vector2, 2> &line, Scalar& t_seg){
    //assumptions:
    //segs are not degenerate
    //not coplanar

    const Scalar& x1 = seg[0][0];
    const Scalar& y1 = seg[0][1];
    const Scalar& x2 = seg[1][0];
    const Scalar& y2 = seg[1][1];

    const Scalar& x3 = line[0][0];
    const Scalar& y3 = line[0][1];
    const Scalar& x4 = line[1][0];
    const Scalar& y4 = line[1][1];

    Scalar n1 = (y3 - y4) * (x1 - x3) + (x4 - x3) * (y1 - y3);
    Scalar d1 = (x4 - x3) * (y1 - y2) - (x1 - x2) * (y4 - y3);
    if(d1 == 0)
        return false;
    t_seg = n1 / d1;
//    Scalar n2 = (y1 - y2) * (x1 - x3) + (x2 - x1) * (y1 - y3);
    Scalar d2 = (x4 - x3) * (y1 - y2) - (x1 - x2) * (y4 - y3);
    if(d2 == 0) {
//        cout<<"d2==0"<<endl;
        return false;
    }

    if (t_seg < 0 || t_seg > 1){
//        cout<<"t_seg = "<<t_seg<<endl;
        return false;
    }

    return true;
}

bool floatTetWild::seg_seg_intersection_2d(const std::array<Vector2, 2> &seg1, const std::array<Vector2, 2> &seg2, Scalar& t2){
    //assumptions:
    //segs are not degenerate
    //not coplanar

    const Scalar& x1 = seg1[0][0];
    const Scalar& y1 = seg1[0][1];
    const Scalar& x2 = seg1[1][0];
    const Scalar& y2 = seg1[1][1];

    const Scalar& x3 = seg2[0][0];
    const Scalar& y3 = seg2[0][1];
    const Scalar& x4 = seg2[1][0];
    const Scalar& y4 = seg2[1][1];

    Scalar n1 = (y3 - y4) * (x1 - x3) + (x4 - x3) * (y1 - y3);
    Scalar d1 = (x4 - x3) * (y1 - y2) - (x1 - x2) * (y4 - y3);
    if(d1 == 0)
        return false;
    Scalar t1 = n1 / d1;
    Scalar n2 = (y1 - y2) * (x1 - x3) + (x2 - x1) * (y1 - y3);
    Scalar d2 = (x4 - x3) * (y1 - y2) - (x1 - x2) * (y4 - y3);
    if(d2 == 0)
        return false;
    t2 = n2 / d2;

    if (t1 < 0 || t1 > 1 || t2 < 0 || t2 > 1)
        return false;

    return true;
}

floatTetWild::Scalar floatTetWild::seg_seg_squared_dist_3d(const std::array<Vector3, 2> &s1, const std::array<Vector3, 2> &s2) {
    Vector3 w0 = s1[0] - s2[0];
    Vector3 u = (s1[1] - s1[0]).normalized();
    Vector3 v = (s2[1] - s2[0]).normalized();
    Scalar a = u.dot(u);
    Scalar b = u.dot(v);
    Scalar c = v.dot(v);
    Scalar d = u.dot(w0);
    Scalar e = v.dot(w0);

    Scalar dd = a * c - b * b;
    Scalar t1, t2;
    if (dd == 0) {
        t1 = 0;
        t2 = d / b;
    } else {
        t1 = (b * e - c * d) / dd;
        t2 = (a * e - b * d) / dd;
    }

    return (((1 - t1) * s1[0] + t1 * s1[1]) - ((1 - t2) * s2[0] + t2 * s2[1])).squaredNorm();
}

floatTetWild::Scalar floatTetWild::p_line_squared_dist_3d(const Vector3 &v, const Vector3 &a, const Vector3 &b) {
    return ((b - a).cross(a - v)).squaredNorm() / (b - a).squaredNorm();
}

floatTetWild::Scalar floatTetWild::p_seg_squared_dist_3d(const Vector3 &v, const Vector3 &a, const Vector3 &b){
    Vector3 av = v-a;
    Vector3 ab = b-a;
    if(av.dot(ab)<0)
        return av.squaredNorm();
    Vector3 bv = v-b;
    if(bv.dot(-ab)<0)
        return bv.squaredNorm();

    return (ab.cross(-av)).squaredNorm()/ab.squaredNorm();
}

bool floatTetWild::seg_plane_intersection(const Vector3& p1, const Vector3& p2, const Vector3& a, const Vector3& n,
                                          Vector3& p, Scalar& d1) {
    Vector3 u = p2 - p1;
    Vector3 w = p1 - a;

    Scalar D = n.dot(u);
    d1 = -n.dot(w);

//    if (fabs(D) <= SCALAR_ZERO) {// segment is parallel to plane
    Scalar t;
    if (fabs(D) == 0) {// segment is parallel to plane
        if (d1 == 0)// segment lies in plane
            t = INT_MAX;
        else
            t = INT_MIN;
//            return 2;
//        else // no intersection
//            return 0;
        return false;
    }

    t = d1 / D;
    if (t <= 0 || t >= 1) {
        return false;
    }

    p = p1 + t * u;
    return true;
}

int floatTetWild::is_tri_tri_cutted(const std::array<Vector3, 3> &f_tri, const std::array<Vector3, 3> &f_tet,
                                    const std::array<int, 3>& oris_tri) {
    int cnt_pos = 0;
    int cnt_neg = 0;
    int cnt_zero = 0;
    for (int j = 0; j < 3; j++) {
        if (oris_tri[j] == Predicates::ORI_ZERO)
            cnt_zero++;
        else if (oris_tri[j] == Predicates::ORI_POSITIVE)
            cnt_pos++;
        else
            cnt_neg++;
    }
    if(cnt_pos == 3 || cnt_neg == 3)
        return CUT_EMPTY;

    if (cnt_zero == 3) {
        int t = get_t(f_tri[0], f_tri[1], f_tri[2]);
        if (is_tri_tri_cutted_2d({{to_2d(f_tri[0], t), to_2d(f_tri[1], t), to_2d(f_tri[2], t)}},
                                 {{to_2d(f_tet[0], t), to_2d(f_tet[1], t), to_2d(f_tet[2], t)}}))
            return CUT_COPLANAR;
        return CUT_EMPTY;
    }

    if(cnt_zero == 2) {
        int t = get_t(f_tri[0], f_tri[1], f_tri[2]);
        for (int j = 0; j < 3; j++) {
            if (oris_tri[j] == oris_tri[(j + 1) % 3]
                && is_seg_tri_cutted_2d({{to_2d(f_tet[j], t), to_2d(f_tet[(j + 1) % 3], t)}},
                                        {{to_2d(f_tri[0], t), to_2d(f_tri[1], t), to_2d(f_tri[2], t)}})) {
                if(j==0)
                    return CUT_EDGE_0;
                if(j==1)
                    return CUT_EDGE_1;
                if(j==2)
                    return CUT_EDGE_2;
            }
        }
        return CUT_EMPTY;
    }

    if(cnt_pos == 0 || cnt_neg == 0)
        return CUT_EMPTY;

    std::array<int, 3> oris_tet;
    int cnt_pos1 = 0;
    int cnt_neg1 = 0;
    int cnt_zero1 = 0;
    for (int j = 0; j < 3; j++) {
        oris_tet[j] = Predicates::orient_3d(f_tet[0], f_tet[1], f_tet[2], f_tri[j]);
        if (oris_tet[j] == Predicates::ORI_ZERO)
            cnt_zero1++;
        else if (oris_tet[j] == Predicates::ORI_POSITIVE)
            cnt_pos1++;
        else
            cnt_neg1++;
    }
//    if(cnt_pos1 == 3 || cnt_neg1 == 3)
    if(cnt_pos1 == 0 || cnt_neg1 == 0)
        return CUT_EMPTY;

    auto is_f1_cut_f2 = [](const std::array<Vector3, 3> &f1, const std::array<Vector3, 3> &f2,
                           const std::array<int, 3>& oris){//check f2 plane
        Vector3 n = ((f2[1] - f2[2]).cross(f2[0] - f2[2])).normalized();
        int t = get_t(f2[0], f2[1], f2[2]);
        for(int j=0;j<3;j++) {
            if (oris[j] == Predicates::ORI_NEGATIVE && oris[(j + 1) % 3] == Predicates::ORI_POSITIVE
                || oris[(j + 1) % 3] == Predicates::ORI_NEGATIVE && oris[j] == Predicates::ORI_POSITIVE) {
                Vector3 p;
                Scalar _;
                if (seg_plane_intersection(f1[j], f1[(j + 1) % 3], n, f2[0], p, _)
                    && is_p_inside_tri_2d(to_2d(p, t), {{to_2d(f2[0], t), to_2d(f2[1], t), to_2d(f2[2], t)}}))
                    return true;
            }
        }
        return false;
    };
    //check tri plane
    if(is_f1_cut_f2(f_tet, f_tri, oris_tri))
        return CUT_FACE;
    //check tet plane
    if(is_f1_cut_f2(f_tri, f_tet, oris_tet))
        return CUT_FACE;

    return CUT_EMPTY;
}

bool floatTetWild::is_tri_tri_cutted_2d(const std::array<Vector2, 3>& vs_tet, const std::array<Vector2, 3>& vs_tri) {
    std::array<int, 9> tri_tet;
    int cnt_pos0 = 0;
    int cnt_neg0 = 0;
    int cnt_pos1 = 0;
    int cnt_neg1 = 0;
    int cnt_pos2 = 0;
    int cnt_neg2 = 0;
    for (int i = 0; i < 3; i++) {
        tri_tet[i * 3] = Predicates::orient_2d(vs_tri[i], vs_tri[(i + 1) % 3], vs_tet[0]);
        if(tri_tet[i * 3] == Predicates::ORI_POSITIVE)
            cnt_pos0++;
        else if(tri_tet[i * 3] == Predicates::ORI_NEGATIVE)
            cnt_neg0++;

        tri_tet[i * 3 + 1] = Predicates::orient_2d(vs_tri[i], vs_tri[(i + 1) % 3], vs_tet[1]);
        if(tri_tet[i * 3 + 1] == Predicates::ORI_POSITIVE)
            cnt_pos1++;
        else if(tri_tet[i * 3 + 1] == Predicates::ORI_NEGATIVE)
            cnt_neg1++;

        tri_tet[i * 3 + 2] = Predicates::orient_2d(vs_tri[i], vs_tri[(i + 1) % 3], vs_tet[2]);
        if(tri_tet[i * 3 + 2] == Predicates::ORI_POSITIVE)
            cnt_pos2++;
        else if(tri_tet[i * 3 + 2] == Predicates::ORI_NEGATIVE)
            cnt_neg2++;
    }
    if(cnt_neg0 == 3 || cnt_pos0 == 3 || cnt_neg1 == 3 || cnt_pos1 == 3 || cnt_neg2 == 3 || cnt_pos2 == 3)
        return true;//one of the vertices is strictly contained inside the other tri
    if (std::find(tri_tet.begin(), tri_tet.end(), Predicates::ORI_NEGATIVE) == tri_tet.end()
        || std::find(tri_tet.begin(), tri_tet.end(), Predicates::ORI_POSITIVE) == tri_tet.end())
        return true;//tet face is strictly contained by tri face


    std::array<int, 9> tet_tri;
    cnt_neg0 = cnt_pos0 = cnt_neg1 = cnt_pos1 = cnt_neg2 = cnt_pos2 = 0;
    for (int i = 0; i < 3; i++) {
        tet_tri[i * 3] = Predicates::orient_2d(vs_tet[i], vs_tet[(i + 1) % 3], vs_tri[0]);
        if(tet_tri[i * 3] == Predicates::ORI_POSITIVE)
            cnt_pos0++;
        else if(tet_tri[i * 3] == Predicates::ORI_NEGATIVE)
            cnt_neg0++;

        tet_tri[i * 3 + 1] = Predicates::orient_2d(vs_tet[i], vs_tet[(i + 1) % 3], vs_tri[1]);
        if(tet_tri[i * 3 + 1] == Predicates::ORI_POSITIVE)
            cnt_pos1++;
        else if(tet_tri[i * 3 + 1] == Predicates::ORI_NEGATIVE)
            cnt_neg1++;

        tet_tri[i * 3 + 2] = Predicates::orient_2d(vs_tet[i], vs_tet[(i + 1) % 3], vs_tri[2]);
        if(tet_tri[i * 3 + 2] == Predicates::ORI_POSITIVE)
            cnt_pos2++;
        else if(tet_tri[i * 3 + 2] == Predicates::ORI_NEGATIVE)
            cnt_neg2++;
    }
    if(cnt_neg0 == 3 || cnt_pos0 == 3 || cnt_neg1 == 3 || cnt_pos1 == 3 || cnt_neg2 == 3 || cnt_pos2 == 3)
        return true;
    if (std::find(tet_tri.begin(), tet_tri.end(), Predicates::ORI_NEGATIVE) == tet_tri.end()
        || std::find(tet_tri.begin(), tet_tri.end(), Predicates::ORI_POSITIVE) == tet_tri.end())
        return true;//tri face is contained by tet face

    for (int tri_e_id = 0; tri_e_id < 3; tri_e_id++) {
        for (int tet_e_id = 0; tet_e_id < 3; tet_e_id++) {
            if (is_crossing(tri_tet[tri_e_id * 3 + tet_e_id], tri_tet[tri_e_id * 3 + (tet_e_id + 1) % 3])
                && is_crossing(tet_tri[tet_e_id * 3 + tri_e_id], tet_tri[tet_e_id * 3 + (tri_e_id + 1) % 3]))
                return true;
        }
    }

    return false;
}

bool floatTetWild::is_seg_tri_cutted_2d(const std::array<Vector2, 2> &seg, const std::array<Vector2, 3> &tri) {
    std::vector<Vector2> ps;
    for(int i=0;i<3;i++){
        Scalar t2;
        if(seg_seg_intersection_2d(seg, {{tri[i], tri[(i+1)%3]}}, t2)){
            ps.push_back(t2 * tri[i] + (1-t2) * tri[(i+1)%3]);//todo double check
        }
    }
    if(ps.size()<2)
        return false;
    for(int i=0;i<ps.size();i++) {
        Vector2 v = ps[i] - ps[(i + 1) % ps.size()];
        if (v[0] < SCALAR_ZERO && v[1] < SCALAR_ZERO) {
            ps.erase(ps.begin() + i);
            i--;
        }
    }
    assert(ps.size()<3);
    if(ps.size()<2)
        return false;

    return true;

    //////////

    std::array<int, 3> tri_seg;
    for (int i = 0; i < 3; i++) {
        tri_seg[i] = Predicates::orient_2d(seg[0], seg[1], tri[i]);
    }
    if (std::find(tri_seg.begin(), tri_seg.end(), Predicates::ORI_NEGATIVE) == tri_seg.end()
        || std::find(tri_seg.begin(), tri_seg.end(), Predicates::ORI_POSITIVE) == tri_seg.end())
        return false; //tri on the same side of seg

    std::array<int, 6> seg_tri;
    int cnt_pos0 = 0;
    int cnt_neg0 = 0;
    int cnt_pos1 = 0;
    int cnt_neg1 = 0;
    for (int i = 0; i < 3; i++) {
        seg_tri[i * 2] = Predicates::orient_2d(tri[i], tri[(i + 1) % 3], seg[0]);
        if(seg_tri[i * 2] == Predicates::ORI_POSITIVE)
            cnt_pos0++;
        else if(seg_tri[i * 2] == Predicates::ORI_NEGATIVE)
            cnt_neg0++;
        seg_tri[i * 2 + 1] = Predicates::orient_2d(tri[i], tri[(i + 1) % 3], seg[1]);
        if(seg_tri[i * 2 + 1] == Predicates::ORI_POSITIVE)
            cnt_pos1++;
        else if(seg_tri[i * 2 + 1] == Predicates::ORI_NEGATIVE)
            cnt_neg1++;
    }
    if (cnt_neg0 == 3 || cnt_pos0 == 3 || cnt_neg1 == 3 || cnt_pos1 == 3) // one of the endpoints is contained inside tri
        return true;

    for (int tri_e_id = 0; tri_e_id < 3; tri_e_id++) {
        if (is_crossing(seg_tri[tri_e_id], seg_tri[(tri_e_id + 1) % 3])
            && is_crossing(tri_seg[tri_e_id], tri_seg[(tri_e_id + 1) % 3]))
            return true;
    }

    return false;
}

bool floatTetWild::is_p_inside_tri_2d(const Vector2& p, const std::array<Vector2, 3> &tri) {
    int cnt_pos = 0;
    int cnt_neg = 0;

    for (int i = 0; i < 3; i++) {
        int ori = Predicates::orient_2d(tri[i], tri[(i + 1) % 3], p);
        if (ori == Predicates::ORI_POSITIVE)
            cnt_pos++;
        else if (ori == Predicates::ORI_NEGATIVE)
            cnt_neg++;
    }
//    if(cnt_neg==0 || cnt_pos==0)
    if (cnt_neg == 3 || cnt_pos == 3) //strict inside
        return true;
    return false;
}

int floatTetWild::get_t(const Vector3 &p0, const Vector3 &p1, const Vector3 &p2) {
    static const std::array<Vector3, 3> ns = {{Vector3(1, 0, 0), Vector3(0, 1, 0), Vector3(0, 0, 1)}};

    Vector3 n = (p1 - p2).cross(p0 - p2);
    Scalar max = 0;
    int t = 0;
    for (int i = 0; i < 3; i++) {
//        Scalar cos_a = abs(n.dot(ns[i]));
        Scalar cos_a = abs(n[i]);
        if (cos_a > max) {
            max = cos_a;
            t = i;
        }
    }
    return t;

//    int t = 2;
//    for (int k = 0; k < 3; k++) {
//        if (p2[k] == p0[k] && p0[k] == p1[k]) {
//            t = k;
//            break;
//        }
//    }
//    return t;
}

floatTetWild::Vector2 floatTetWild::to_2d(const Vector3 &p, int t) {
    return Vector2(p[(t + 1) % 3], p[(t + 2) % 3]);
}

floatTetWild::Vector2 floatTetWild::to_2d(const Vector3 &p, const Vector3& n, const Vector3& pp, int t) {
    Scalar dist = n.dot(p - pp);
    Vector3 proj_p = p - dist * n;
    return Vector2(proj_p[(t + 1) % 3], proj_p[(t + 2) % 3]);
}

bool floatTetWild::is_crossing(int s1, int s2) {
    if (s1 == Predicates::ORI_POSITIVE && s2 == Predicates::ORI_NEGATIVE
        || s2 == Predicates::ORI_POSITIVE && s1 == Predicates::ORI_NEGATIVE)
        return true;
    return false;
}



extern "C++" int tri_tri_intersection_test_3d(floatTetWild::Scalar p1[3], floatTetWild::Scalar q1[3], floatTetWild::Scalar r1[3],
                                             floatTetWild::Scalar p2[3], floatTetWild::Scalar q2[3], floatTetWild::Scalar r2[3],
                                             int * coplanar,
                                             floatTetWild::Scalar source[3], floatTetWild::Scalar target[3]);

int floatTetWild::is_tri_tri_cutted(const Vector3& p1, const Vector3& p2, const Vector3& p3,
                                    const Vector3& q1, const Vector3& q2, const Vector3& q3) {
    std::array<Scalar, 3> p_1 = {{0, 0, 0}}, q_1 = {{0, 0, 0}}, r_1 = {{0, 0, 0}};
    std::array<Scalar, 3> p_2 = {{0, 0, 0}}, q_2 = {{0, 0, 0}}, r_2 = {{0, 0, 0}};
    int coplanar = 0;
    std::array<Scalar, 3> s = {{0,0,0}}, t = {{0,0,0}};
    for (int j = 0; j < 3; j++) {
        p_1[j] = p1[j];
        q_1[j] = p2[j];
        r_1[j] = p3[j];
        p_2[j] = q1[j];
        q_2[j] = q2[j];
        r_2[j] = q3[j];
    }

    if (!tri_tri_intersection_test_3d(&p_1[0], &q_1[0], &r_1[0], &p_2[0], &q_2[0], &r_2[0], &coplanar, &s[0], &t[0]))
        return CUT_EMPTY;

    if(coplanar == 1) {
        // if(tri_tri_overlap_test_3d(&p_1[0], &q_1[0], &r_1[0], &p_2[0], &q_2[0], &r_2[0])) //TODO?
            return CUT_COPLANAR;
    }

    if (s[0] == t[0] && s[1] == t[1] && s[2] == t[2])
        return CUT_EMPTY;

    auto is_collinear = [&](const Vector3 &st, const Vector3 &q12) {
        Vector3 cross = st.cross(q12);
        if (cross[0] <= SCALAR_ZERO && cross[1] <= SCALAR_ZERO && cross[2] <= SCALAR_ZERO)
            return true;
        return false;
    };

    Vector3 st(t[0] - s[0], t[1] - s[1], t[2] - s[2]);
    if (is_collinear(st, q1 - q2))
        return CUT_EDGE_0;
    if (is_collinear(st, q2 - q3))
        return CUT_EDGE_1;
    if (is_collinear(st, q3 - q1))
        return CUT_EDGE_2;

    return CUT_FACE;
}

int floatTetWild::is_tri_tri_cutted_hint(const Vector3& p1, const Vector3& p2, const Vector3& p3,
                                         const Vector3& q1, const Vector3& q2, const Vector3& q3, int hint, bool is_debug) {
    std::array<Scalar, 3> p_1 = {{0, 0, 0}}, q_1 = {{1, 0, 0}}, r_1 = {{0, 1, 0}};
    std::array<Scalar, 3> p_2 = {{0, 0, 0}}, q_2 = {{-1, -1, 0}}, r_2 = {{1, 1, 0}};
    int coplanar = 0;
    std::array<Scalar, 3> s = {{0,0,0}}, t = {{0,0,0}};

    for (int j = 0; j < 3; j++) {
        p_1[j] = p1[j];
        q_1[j] = p2[j];
        r_1[j] = p3[j];
        p_2[j] = q1[j];
        q_2[j] = q2[j];
        r_2[j] = q3[j];
    }

    if(hint == CUT_COPLANAR){
//        if(tri_tri_overlap_test_3d(&p_1[0], &q_1[0], &r_1[0], &p_2[0], &q_2[0], &r_2[0]))
//            return CUT_COPLANAR;

        //todo: to2d
//        auto to_2d = [](const Vector3& v, int t) {
//            return Vector2(v[(t + 1) % 3], v[(t + 2) % 3]);
//        };
//        int t = 2;
//        for (int k = 0; k < 3; k++) {
//            if (p3[k] == p1[k] && p1[k] == p2[k]) {
//                t = k;
//                break;
//            }
//        }

        int t = get_t(p1, p2, p3);

        if(is_tri_tri_cutted_2d({{to_2d(p1, t), to_2d(p2, t), to_2d(p3, t)}}, {{to_2d(q1, t), to_2d(q2, t), to_2d(q3, t)}}))
            return CUT_COPLANAR;
        return CUT_EMPTY;
    }

    if(is_debug){
//        p_1 = {{22, -27.478179999999998, 3.5}};
//        q_1 = {{25.5, -27.478179999999998, 0}};
//        r_1 = {{-25.5, -27.478179999999998, 0}};
//        p_2 = {{28.050000000000001, -30.028179999999999, -2.5500000000000003}};
//        q_2 = {{25.5, -13.183373550207467, 0}};
//        r_2 = {{25.5, -19.059304172821577, 0}};
//
//        cout<<"res = "<<tri_tri_intersection_test_3d(&p_1[0], &q_1[0], &r_1[0], &p_2[0], &q_2[0], &r_2[0], &coplanar, &s[0], &t[0])<<endl;
//        cout<<std::setprecision(17);
//        cout<<"s"<<endl;
//        cout<<s[0]<<" "<<s[1]<<" "<<s[2]<<endl;
//        cout<<"t"<<endl;
//        cout<<t[0]<<" "<<t[1]<<" "<<t[2]<<endl;
//        cout<<"coplanar = "<<coplanar<<endl;
//        cout<<"pqr 1"<<endl;
//        for(int j=0;j<3;j++)
//            cout<<p_1[j]<<" ";
//        cout<<endl;
//        for(int j=0;j<3;j++)
//            cout<<q_1[j]<<" ";
//        cout<<endl;
//        for(int j=0;j<3;j++)
//            cout<<r_1[j]<<" ";
//        cout<<endl;
//        cout<<"pqr 2"<<endl;
//        for(int j=0;j<3;j++)
//            cout<<p_2[j]<<" ";
//        cout<<endl;
//        for(int j=0;j<3;j++)
//            cout<<q_2[j]<<" ";
//        cout<<endl;
//        for(int j=0;j<3;j++)
//            cout<<r_2[j]<<" ";
//        cout<<endl;
        //pausee();
    }

    int result = tri_tri_intersection_test_3d(&p_1[0], &q_1[0], &r_1[0], &p_2[0], &q_2[0], &r_2[0], &coplanar, &s[0], &t[0]);
    if(is_debug){
        cout<<">>result = "<<result<<endl;
    }
    if (result != 1) {
//        if (coplanar != 0) {
//            cout<<"CUT_EMPTY "<<result<<endl;
//            //pausee();
//        }
        return CUT_EMPTY;
    }
//    if (coplanar != 0) {
//        cout<<"NOT CUT_EMPTY "<<coplanar<<" "<<result<<endl;
//        //pausee();
//    }

    if (std::abs(s[0] - t[0]) <= SCALAR_ZERO && std::abs(s[1] - t[1]) <= SCALAR_ZERO && std::abs(s[2] - t[2]) <= SCALAR_ZERO)
        return CUT_EMPTY;

//    cout<<"s"<<s[0]<<" "<<s[1]<<" "<<s[2]<<endl;
//    cout<<"t"<<t[0]<<" "<<t[1]<<" "<<t[2]<<endl;
//    cout<<(s[0]-t[0])<<" "<<(s[1]-t[1])<<" "<<(s[2]-t[2])<<endl;

    if (hint == CUT_EDGE_0)
        return CUT_EDGE_0;
    if (hint == CUT_EDGE_1)
        return CUT_EDGE_1;
    if (hint == CUT_EDGE_2)
        return CUT_EDGE_2;

    return CUT_FACE;
}

void floatTetWild::get_bbox_face(const Vector3& p0, const Vector3& p1, const Vector3& p2,
        Vector3& min, Vector3& max, Scalar eps) {
    min = p0;
    max = p0;
    for (int j = 0; j < 3; j++) {
        Scalar tmp_min = std::min(p1[j], p2[j]);
        if (tmp_min < min[j])
            min[j] = tmp_min;
        Scalar tmp_max = std::max(p1[j], p2[j]);
        if (tmp_max > max[j])
            max[j] = tmp_max;
    }

    if (eps != 0) {
        for (int j = 0; j < 3; j++) {
            min[j] -= eps;
            max[j] += eps;
        }
    }
}

void floatTetWild::get_bbox_tet(const Vector3& p0, const Vector3& p1, const Vector3& p2, const Vector3& p3,
        Vector3& min, Vector3& max, Scalar eps) {
    min = p0;
    max = p0;
    for (int j = 0; j < 3; j++) {
        if (p1[j] < min[j])
            min[j] = p1[j];
        if (p2[j] < min[j])
            min[j] = p2[j];
        if (p3[j] < min[j])
            min[j] = p3[j];

        if (p1[j] > max[j])
            max[j] = p1[j];
        if (p2[j] > max[j])
            max[j] = p2[j];
        if (p3[j] > max[j])
            max[j] = p3[j];
    }

    if (eps != 0) {
        for (int j = 0; j < 3; j++) {
            min[j] -= eps;
            max[j] += eps;
        }
    }
}

bool floatTetWild::is_bbox_intersected(const Vector3& min1, const Vector3& max1, const Vector3& min2, const Vector3& max2) {
    for (int j = 0; j < 3; j++) {
        if (min1[j] > max2[j])
            return false;
        if (max1[j] < min2[j])
            return false;
    }
    return true;
}

bool floatTetWild::is_tri_inside_tet(const std::array<Vector3, 3>& ps,
        const Vector3& p0t, const Vector3& p1t, const Vector3& p2t, const Vector3& p3t) {
    int cnt_pos = 0;
    int cnt_neg = 0;

    for (int i = 0; i < 3; i++) {
        int ori = Predicates::orient_3d(ps[i], p1t, p2t, p3t);
        if (ori == Predicates::ORI_POSITIVE)
            cnt_pos++;
        else if (ori == Predicates::ORI_NEGATIVE)
            cnt_neg++;

        ori = Predicates::orient_3d(p0t, ps[i], p2t, p3t);
        if (ori == Predicates::ORI_POSITIVE)
            cnt_pos++;
        else if (ori == Predicates::ORI_NEGATIVE)
            cnt_neg++;

        ori = Predicates::orient_3d(p0t, p1t, ps[i], p3t);
        if (ori == Predicates::ORI_POSITIVE)
            cnt_pos++;
        else if (ori == Predicates::ORI_NEGATIVE)
            cnt_neg++;

        ori = Predicates::orient_3d(p0t, p1t, p2t, ps[i]);
        if (ori == Predicates::ORI_POSITIVE)
            cnt_pos++;
        else if (ori == Predicates::ORI_NEGATIVE)
            cnt_neg++;
    }

    if(cnt_pos == 0 || cnt_neg == 0)
        return true;

    return false;
}

bool floatTetWild::is_point_inside_tet(const Vector3& p, const Vector3& p0t, const Vector3& p1t, const Vector3& p2t, const Vector3& p3t) {///inside or on
    int cnt_pos = 0;
    int cnt_neg = 0;

    int ori = Predicates::orient_3d(p, p1t, p2t, p3t);
    if (ori == Predicates::ORI_POSITIVE)
        cnt_pos++;
    else if (ori == Predicates::ORI_NEGATIVE)
        cnt_neg++;

    ori = Predicates::orient_3d(p0t, p, p2t, p3t);
    if (ori == Predicates::ORI_POSITIVE)
        cnt_pos++;
    else if (ori == Predicates::ORI_NEGATIVE)
        cnt_neg++;

    ori = Predicates::orient_3d(p0t, p1t, p, p3t);
    if (ori == Predicates::ORI_POSITIVE)
        cnt_pos++;
    else if (ori == Predicates::ORI_NEGATIVE)
        cnt_neg++;

    ori = Predicates::orient_3d(p0t, p1t, p2t, p);
    if (ori == Predicates::ORI_POSITIVE)
        cnt_pos++;
    else if (ori == Predicates::ORI_NEGATIVE)
        cnt_neg++;

    if (cnt_pos == 0 || cnt_neg == 0)
        return true;

    return false;
}