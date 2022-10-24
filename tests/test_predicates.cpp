////////////////////////////////////////////////////////////////////////////////
#include <floattetwild/Predicates.hpp>
#include <catch2/catch_all.hpp>
#include <iostream>
////////////////////////////////////////////////////////////////////////////////

using namespace floatTetWild;

extern "C++" int tri_tri_intersection_test_3d(floatTetWild::Scalar p1[3], floatTetWild::Scalar q1[3], floatTetWild::Scalar r1[3],
                                             floatTetWild::Scalar p2[3], floatTetWild::Scalar q2[3], floatTetWild::Scalar r2[3],
                                             int * coplanar,
                                             floatTetWild::Scalar source[3], floatTetWild::Scalar target[3]);

TEST_CASE("orient_3d", "[predicates]") {
	Vector3 p1; p1 << 0, 0, 0;
	Vector3 p2; p2 << 1, 0, 0;
	Vector3 p3; p3 << 0, 1, 0;
	Vector3 p4; p4 << 0, 0, 1;

	int orientation;

	orientation = Predicates::orient_3d(p1, p2, p3, p4);
	REQUIRE(orientation == Predicates::ORI_NEGATIVE);


	orientation = Predicates::orient_3d(p1,p4, p3, p2);
	REQUIRE(orientation == Predicates::ORI_POSITIVE);
}


TEST_CASE("tri_tri_intersection_test_3d_rounded", "[predicates]") {
    Scalar p1[3] = {22, -27.47818, 3.5};
    Scalar q1[3] = {25.5, -27.47818, 0};
    Scalar r1[3] = {-25.5, -27.47818, 0};

    Scalar p2[3]= {28.05000000000000, -30.02818, -2.55};
    Scalar q2[3]= {25.5, -13.183373550207467, 0};
    Scalar r2[3]= {25.5, -19.059304172821577, 0};

    int coplanar = 0;
    Scalar s[3] = {0, 0, 0};
    Scalar t[3] = {0, 0, 0};
    int res = tri_tri_intersection_test_3d(
      p1, q1, r1,
      p2, q2, r2,
      &coplanar,
      s, t
      );

    REQUIRE(coplanar == 0);
    REQUIRE(res == 0);
}


TEST_CASE("tri_tri_intersection_test_3d_floating", "[predicates]") {
    Scalar p1[3] = {22, -27.478179999999998, 3.5};
    Scalar q1[3] = {25.5, -27.478179999999998, 0};
    Scalar r1[3] = {-25.5, -27.478179999999998, 0};

    Scalar p2[3]= {28.050000000000001, -30.028179999999999, -2.5500000000000003};
    Scalar q2[3]= {25.5, -13.183373550207467, 0};
    Scalar r2[3]= {25.5, -19.059304172821577, 0};

    int coplanar = 0;
    Scalar s[3] = {0, 0, 0};
    Scalar t[3] = {0, 0, 0};
    int res = tri_tri_intersection_test_3d(
      p1, q1, r1,
      p2, q2, r2,
      &coplanar,
      s, t
      );



    REQUIRE(coplanar == 0);
    REQUIRE(res == 666);
}
