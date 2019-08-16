#include <floattetwild/Predicates.hpp>

extern "C" floatTetWild::Scalar orient3d(const floatTetWild::Scalar *pa, const floatTetWild::Scalar *pb, const floatTetWild::Scalar *pc, const floatTetWild::Scalar *pd);
extern "C" floatTetWild::Scalar orient2d(const floatTetWild::Scalar *pa, const floatTetWild::Scalar *pb, const floatTetWild::Scalar *pc);

#include <geogram/delaunay/delaunay_3d.h>
namespace floatTetWild {
#define GEO_PREDICATES true
	const int Predicates::ORI_POSITIVE;
	const int Predicates::ORI_ZERO;
	const int Predicates::ORI_NEGATIVE;
	const int Predicates::ORI_UNKNOWN;

	int Predicates::orient_3d(const Vector3& p1, const Vector3& p2, const Vector3& p3, const Vector3& p4) {
#if GEO_PREDICATES
		const int result = -GEO::PCK::orient_3d(p1.data(), p2.data(), p3.data(), p4.data());
#else
		const Scalar result = orient3d(p1.data(), p2.data(), p3.data(), p4.data());
#endif

//		if (result > SCALAR_ZERO)
//			return ORI_POSITIVE;
//		else if (result < -SCALAR_ZERO)
//			return ORI_NEGATIVE;
//		else
//			return ORI_ZERO;

		if (result > 0)
			return ORI_POSITIVE;
		else if (result < 0)
			return ORI_NEGATIVE;
		else
			return ORI_ZERO;
	}

	int Predicates::orient_3d_tolerance(const Vector3& p1, const Vector3& p2, const Vector3& p3, const Vector3& p) {
		const Scalar result = orient3d(p1.data(), p2.data(), p3.data(), p.data());
		if (result == 0)
			return ORI_ZERO;

		Vector3 n = ((p2 - p3).cross(p1 - p3)).normalized();
		Scalar d = std::abs(n.dot(p - p1));
		if (d <= SCALAR_ZERO)
			return Predicates::ORI_ZERO;

		if (result > 0)
			return ORI_POSITIVE;
		else
			return ORI_NEGATIVE;
	}

	Scalar Predicates::orient_3d_volume(const Vector3& p1, const Vector3& p2, const Vector3& p3, const Vector3& p4) {
#if GEO_PREDICATES
		const int ori = -GEO::PCK::orient_3d(p1.data(), p2.data(), p3.data(), p4.data());
#else
		const Scalar ori = orient3d(p1.data(), p2.data(), p3.data(), p4.data());
#endif
		if (ori <= 0)
			return ori;
		else
			return (p1 - p4).dot((p2 - p4).cross(p3 - p4)) / 6;
	}

	int Predicates::orient_2d(const Vector2& p1, const Vector2& p2, const Vector2& p3) {
#if GEO_PREDICATES
		const int result = -GEO::PCK::orient_2d(p1.data(), p2.data(), p3.data());
#else
		const Scalar result = orient2d(p1.data(), p2.data(), p3.data());
#endif
		if (result > 0)
			return ORI_POSITIVE;
		else if (result < 0)
			return ORI_NEGATIVE;
		else
			return ORI_ZERO;
	}

} // namespace floatTetWild
