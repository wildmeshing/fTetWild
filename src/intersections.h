// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#ifndef FLOATTETWILD_INTERSECTIONS_H
#define FLOATTETWILD_INTERSECTIONS_H

#include <floattetwild/Mesh.hpp>

namespace floatTetWild {
#define CUT_EDGE_0 0
#define CUT_EDGE_1 1
#define CUT_EDGE_2 2
#define CUT_FACE 3
#define CUT_COPLANAR 4
#define CUT_EMPTY -1

    int is_tri_tri_cutted(const std::array<Vector3, 3> &f_tri, const std::array<Vector3, 3> &f_tet,
                          const std::array<int, 3>& oris_tri);

    Scalar seg_seg_squared_dist_3d(const std::array<Vector3, 2> &s1, const std::array<Vector3, 2> &s2);

    Scalar p_seg_squared_dist_3d(const Vector3 &p, const Vector3 &a, const Vector3 &b);
    Scalar p_line_squared_dist_3d(const Vector3 &p, const Vector3 &a, const Vector3 &b);

    bool is_p_inside_tri_2d(const Vector2& p, const std::array<Vector2, 3> &tri);
    bool is_seg_tri_cutted_2d(const std::array<Vector2, 2> &seg, const std::array<Vector2, 3> &tri);
    bool is_tri_tri_cutted_2d(const std::array<Vector2, 3> &p_tet, const std::array<Vector2, 3> &p_tri);

    bool seg_seg_intersection_2d(const std::array<Vector2, 2> &seg1, const std::array<Vector2, 2> &seg2, Scalar& t2);
    bool seg_line_intersection_2d(const std::array<Vector2, 2> &seg, const std::array<Vector2, 2> &line, Scalar& t_seg);
    bool seg_plane_intersection(const Vector3 &p1, const Vector3 &p2, const Vector3 &a, const Vector3 &n,
                                Vector3 &p, Scalar &d1);

    int get_t(const Vector3 &p0, const Vector3 &p1, const Vector3 &p2);
    Vector2 to_2d(const Vector3 &p, int t);
    Vector2 to_2d(const Vector3 &p, const Vector3& n, const Vector3& pp, int t);

    bool is_crossing(int s1, int s2);

    int is_tri_tri_cutted(const Vector3 &p1, const Vector3 &p2, const Vector3 &p3,//cutting tri
                          const Vector3 &q1, const Vector3 &q2, const Vector3 &q3);//face of tet
    int is_tri_tri_cutted_hint(const Vector3 &p1, const Vector3 &p2, const Vector3 &p3,//cutting tri
                               const Vector3 &q1, const Vector3 &q2, const Vector3 &q3, int hint,
                               bool is_debug = false);//face of tet

    void get_bbox_face(const Vector3& p0, const Vector3& p1, const Vector3& p2, Vector3& min, Vector3& max, Scalar eps = 0);
    void get_bbox_tet(const Vector3& p0, const Vector3& p1, const Vector3& p2, const Vector3& p3, Vector3& min, Vector3& max, Scalar eps = 0);

    bool is_bbox_intersected(const Vector3& min1, const Vector3& max1, const Vector3& min2, const Vector3& max2);

    bool is_tri_inside_tet(const std::array<Vector3, 3>& ps,
                           const Vector3& p0t, const Vector3& p1t, const Vector3& p2t, const Vector3& p3t);
    bool is_point_inside_tet(const Vector3& p, const Vector3& p0t, const Vector3& p1t, const Vector3& p2t, const Vector3& p3t);
}

#endif //FLOATTETWILD_INTERSECTIONS_H
