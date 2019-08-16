#pragma once

// to generate delaunay mesh using the 3d model along with points of the box

#include <floattetwild/Mesh.hpp>
#include <floattetwild/Types.hpp>

#include <geogram/delaunay/delaunay_3d.h>
#include <floattetwild/AABBWrapper.h>

namespace floatTetWild {

	class FloatTetDelaunay {
	public:
		// to generate delaunay meshes with the model vertices and bounding box points.
		static void tetrahedralize(const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces, const AABBWrapper &tree,
								   Mesh &mesh, std::vector<bool> &is_face_inserted);

	};
}
