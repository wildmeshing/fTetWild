////////////////////////////////////////////////////////////////////////////////
#include <catch2/catch.hpp>
#include <iostream>
#include <floattetwild/ElementSubdivision.h>
////////////////////////////////////////////////////////////////////////////////

using namespace floatTetWild;
void SLtest() {
	const std::array<Vector3, 3> cutface = { Vector3(1,1,1), Vector3(1, 1,0.5), Vector3(0.11, 1, 0.9) };
	//cutface[0] << 0,0,0;



	const Vector3 linepoints0(0, 0, 0);
	const Vector3 linepoints1(1, 1, 1);
	int cutOrnot;
	Vector3 interp;
	Scalar t;
	int ion;
	Parameters params;
	Subdivision::SLIntersection(params, cutface, linepoints0, linepoints1, cutOrnot, interp, t, ion);
	std::cout << cutOrnot << "\n" << ion << "\n" << t << std::endl;
}



TEST_CASE("surface cut line", "[subdivisio]]]]") {
	SLtest();
}


