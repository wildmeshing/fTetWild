// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Teseo Schneider <teseo.schneider@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#pragma once

#include <floattetwild/Types.hpp>

#include <array>
#include <vector>

namespace floatTetWild {
	class CutTable {
public:
		static const std::vector<std::vector<Vector4i>>& get_tet_confs(const int idx);
		static const std::vector<std::vector<Vector2i>>& get_diag_confs(const int idx);
		static const std::vector<std::vector<std::array<bool, 4>>>& get_surface_conf(const int idx);
		static const std::vector<std::vector<Vector4i>>& get_face_id_conf(const int idx);
		static inline const std::vector<Vector4i>& get_tet_conf(const int idx, const int cfg){ return get_tet_confs(idx)[cfg]; }
		static inline const std::vector<std::array<bool, 4>>& get_surface_conf(const int idx, const int cfg){ return get_surface_conf(idx)[cfg]; }
		static inline const std::vector<Vector4i>& get_face_id_conf(const int idx, const int cfg){ return get_face_id_conf(idx)[cfg]; }
	};
}
