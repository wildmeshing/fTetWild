#pragma once

#include <floattetwild/Types.hpp>

namespace floatTetWild {

	class Predicates
	{
	public:
		static const int ORI_POSITIVE = 1;
        static const int ORI_ZERO = 0;
        static const int ORI_NEGATIVE = -1;
		static const int ORI_UNKNOWN = INT_MAX;

		static int orient_3d(const Vector3& p1, const Vector3& p2, const Vector3& p3, const Vector3& p4);
		static int orient_3d_tolerance(const Vector3& p1, const Vector3& p2, const Vector3& p3, const Vector3& p4);
		static Scalar orient_3d_volume(const Vector3& p1, const Vector3& p2, const Vector3& p3, const Vector3& p4);

        static int orient_2d(const Vector2& p1, const Vector2& p2, const Vector2& p3);
	};

} // namespace floattetwild
