/*
*  Triangle-Triangle Overlap Test Routines
*  July, 2002
*  Updated December 2003
*
*  This file contains C implementation of algorithms for
*  performing two and three-dimensional triangle-triangle intersection test
*  The algorithms and underlying theory are described in
*
* "Fast and Robust Triangle-Triangle Overlap Test
*  Using Orientation Predicates"  P. Guigue - O. Devillers
*
*  Journal of Graphics Tools, 8(1), 2003
*
*  Several geometric predicates are defined.  Their parameters are all
*  points.  Each point is an array of two or three real precision
*  floating point numbers. The geometric predicates implemented in
*  this file are:
*
*
*    int tri_tri_intersection_test_3d(p1,q1,r1,p2,q2,r2,
*                                     coplanar,source,target)
*
*       is a version that computes the segment of intersection when
*       the triangles overlap (and are not coplanar)
*
*    each function returns 1 if the triangles (including their
*    boundary) intersect, otherwise 0
*
*
*  Other information are available from the Web page
*  http:<i>//www.acm.org/jgt/papers/GuigueDevillers03/
*
*/

// modified by Aaron to better detect coplanarity


#ifdef FLOAT_TETWILD_USE_FLOAT
    typedef float real;                      // float
#else
    typedef double real;                      // double
#endif

/* function prototype */



int tri_tri_intersection_test_3d(real p1[3], real q1[3], real r1[3],
                                 real p2[3], real q2[3], real r2[3],
                                 int * coplanar,
                                 real source[3],real target[3]);

int sub_sub_cross_sub_dot(real a[3], real b[3], real c[3], real d[3]);



/* coplanar returns whether the triangles are coplanar
*  source and target are the endpoints of the segment of
*  intersection if it exists)
*/


/* some 3D macros */

#define CROSS(dest,v1,v2)                       \
               dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
               dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
               dest[2]=v1[0]*v2[1]-v1[1]*v2[0];


#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])



#define SUB(dest,v1,v2) dest[0]=v1[0]-v2[0]; \
                        dest[1]=v1[1]-v2[1]; \
                        dest[2]=v1[2]-v2[2];


#define SCALAR(dest,alpha,v) dest[0] = alpha * v[0]; \
                             dest[1] = alpha * v[1]; \
                             dest[2] = alpha * v[2];


/*
*
*  Three-dimensional Triangle-Triangle Intersection
*
*/

/*
   This macro is called when the triangles surely intersect
   It constructs the segment of intersection of the two triangles
   if they are not coplanar.
*/

#define CONSTRUCT_INTERSECTION(p1,q1,r1,p2,q2,r2) { \
  if (sub_sub_cross_sub_dot(q1, r2, p1, p2) > 0) {\
    if (sub_sub_cross_sub_dot(r1, r2, p1, p2) <=  0) {\
      if (sub_sub_cross_sub_dot(r1, q2, p1, p2) >  0) {\
         SUB(v1,p1,p2) \
         SUB(v2,p1,r1) \
         alpha = DOT(v1,N2) / DOT(v2,N2); \
         SCALAR(v1,alpha,v2) \
         SUB(source,p1,v1) \
         SUB(v1,p2,p1) \
         SUB(v2,p2,r2) \
         alpha = DOT(v1,N1) / DOT(v2,N1); \
         SCALAR(v1,alpha,v2) \
         SUB(target,p2,v1) \
         return 1; \
      }\
      else { \
         SUB(v1,p2,p1) \
         SUB(v2,p2,q2) \
         alpha = DOT(v1,N1) / DOT(v2,N1); \
         SCALAR(v1,alpha,v2) \
         SUB(source,p2,v1) \
         SUB(v1,p2,p1) \
         SUB(v2,p2,r2) \
         alpha = DOT(v1,N1) / DOT(v2,N1); \
         SCALAR(v1,alpha,v2) \
         SUB(target,p2,v1) \
         return 1; \
      } \
    } \
    else { \
      return 0; \
    } \
  } \
  else { \
    if (sub_sub_cross_sub_dot(q1, q2, p1, p2) <  0) {\
      return 0; \
    } \
    else { \
      if (sub_sub_cross_sub_dot(r1, q2, p1, p2) >= 0) {\
        SUB(v1,p1,p2) \
        SUB(v2,p1,r1) \
        alpha = DOT(v1,N2) / DOT(v2,N2); \
        SCALAR(v1,alpha,v2) \
        SUB(source,p1,v1) \
        SUB(v1,p1,p2) \
        SUB(v2,p1,q1) \
        alpha = DOT(v1,N2) / DOT(v2,N2); \
        SCALAR(v1,alpha,v2) \
        SUB(target,p1,v1) \
        return 1; \
      } \
      else { \
        SUB(v1,p2,p1) \
        SUB(v2,p2,q2) \
        alpha = DOT(v1,N1) / DOT(v2,N1); \
        SCALAR(v1,alpha,v2) \
        SUB(source,p2,v1) \
        SUB(v1,p1,p2) \
        SUB(v2,p1,q1) \
        alpha = DOT(v1,N2) / DOT(v2,N2); \
        SCALAR(v1,alpha,v2) \
        SUB(target,p1,v1) \
        return 1; \
      } \
    } \
  } \
}



#define TRI_TRI_INTER_3D(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2) { \
  if (dp2 > 0) { \
     if (dq2 > 0) CONSTRUCT_INTERSECTION(p1,r1,q1,r2,p2,q2) \
     else if (dr2 > 0) CONSTRUCT_INTERSECTION(p1,r1,q1,q2,r2,p2)\
     else CONSTRUCT_INTERSECTION(p1,q1,r1,p2,q2,r2) }\
  else if (dp2 < 0) { \
    if (dq2 < 0) CONSTRUCT_INTERSECTION(p1,q1,r1,r2,p2,q2)\
    else if (dr2 < 0) CONSTRUCT_INTERSECTION(p1,q1,r1,q2,r2,p2)\
    else CONSTRUCT_INTERSECTION(p1,r1,q1,p2,q2,r2)\
  } else { \
    if (dq2 < 0) { \
      if (dr2 >= 0)  CONSTRUCT_INTERSECTION(p1,r1,q1,q2,r2,p2)\
      else CONSTRUCT_INTERSECTION(p1,q1,r1,p2,q2,r2)\
    } \
    else if (dq2 > 0) { \
      if (dr2 > 0) CONSTRUCT_INTERSECTION(p1,r1,q1,p2,q2,r2)\
      else  CONSTRUCT_INTERSECTION(p1,q1,r1,q2,r2,p2)\
    } \
    else  { \
      if (dr2 > 0) CONSTRUCT_INTERSECTION(p1,q1,r1,r2,p2,q2)\
      else if (dr2 < 0) CONSTRUCT_INTERSECTION(p1,r1,q1,r2,p2,q2)\
      else { \
        *coplanar = 1; \
       return -1;\
     } \
  }} }


/*
   The following version computes the segment of intersection of the
   two triangles if it exists.
   coplanar returns whether the triangles are coplanar
   source and target are the endpoints of the line segment of intersection
*/

//extern "C" real orient3d(const real *pa, const real *pb, const real *pc, const real *pd);

#include <geogram/delaunay/delaunay_3d.h>
inline int sub_sub_cross_sub_dot(real pa[3], real pb[3], real pc[3], real pd[3]) {
//    const real result = orient3d(pa, pb, pc, pd);
    auto result = -GEO::PCK::orient_3d(pa, pb, pc, pd);
    if (result > 0)
        return 1;
    else if (result < 0)
        return -1;
    return 0;
}

int tri_tri_intersection_test_3d(real p1[3], real q1[3], real r1[3],
                                 real p2[3], real q2[3], real r2[3],
                                 int* coplanar,
                                 real source[3], real target[3] )

{
    int dp1, dq1, dr1, dp2, dq2, dr2;
    real v1[3], v2[3], v[3];
    real N1[3], N2[3], N[3];
    real alpha;

	SUB(v1,q1,p1)
    SUB(v2,r1,p1)
    CROSS(N1,v1,v2)

	SUB(v1,p2,r2)
    SUB(v2,q2,r2)
    CROSS(N2,v1,v2)


    *coplanar = 0;

    // Compute distance signs  of p1, q1 and r1
    // to the plane of triangle(p2,q2,r2)


    dp1 = sub_sub_cross_sub_dot(p2, q2, r2, p1);
    dq1 = sub_sub_cross_sub_dot(p2, q2, r2, q1);
    dr1 = sub_sub_cross_sub_dot(p2, q2, r2, r1);

    if (((dp1 * dq1) > 0) && ((dp1 * dr1) > 0))  return 666;

    // Compute distance signs  of p2, q2 and r2
    // to the plane of triangle(p1,q1,r1)

    dp2 = sub_sub_cross_sub_dot(p1, q1, r1, p2);
    dq2 = sub_sub_cross_sub_dot(p1, q1, r1, q2);
    dr2 = sub_sub_cross_sub_dot(p1, q1, r1, r2);

    if (((dp2 * dq2) > 0) && ((dp2 * dr2) > 0)) return 666;

    // Permutation in a canonical form of T1's vertices

    if (dp1 > 0) {
        if (dq1 > 0) TRI_TRI_INTER_3D(r1,p1,q1,p2,r2,q2,dp2,dr2,dq2)
        else if (dr1 > 0) TRI_TRI_INTER_3D(q1,r1,p1,p2,r2,q2,dp2,dr2,dq2)

        else TRI_TRI_INTER_3D(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2)
    } else if (dp1 < 0) {
        if (dq1 < 0) TRI_TRI_INTER_3D(r1,p1,q1,p2,q2,r2,dp2,dq2,dr2)
        else if (dr1 < 0) TRI_TRI_INTER_3D(q1,r1,p1,p2,q2,r2,dp2,dq2,dr2)
        else TRI_TRI_INTER_3D(p1,q1,r1,p2,r2,q2,dp2,dr2,dq2)
    } else {
        if (dq1 < 0) {
            if (dr1 >= 0) TRI_TRI_INTER_3D(q1,r1,p1,p2,r2,q2,dp2,dr2,dq2)
            else TRI_TRI_INTER_3D(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2)
        }
        else if (dq1 > 0) {
            if (dr1 > 0) TRI_TRI_INTER_3D(p1,q1,r1,p2,r2,q2,dp2,dr2,dq2)
            else TRI_TRI_INTER_3D(q1,r1,p1,p2,q2,r2,dp2,dq2,dr2)
        }
        else  {
            if (dr1 > 0) TRI_TRI_INTER_3D(r1,p1,q1,p2,q2,r2,dp2,dq2,dr2)
            else if (dr1 < 0) TRI_TRI_INTER_3D(r1,p1,q1,p2,r2,q2,dp2,dr2,dq2)
            else {
                // triangles are co-planar
				
                *coplanar = 1;
			
                return -1;
            }
        }
    }
};
