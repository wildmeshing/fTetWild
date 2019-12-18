// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Teseo Schneider <teseo.schneider@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include "auto_table.hpp"

#include <cassert>

namespace floatTetWild {
	const std::vector<std::vector<Vector4i>>& CutTable::get_tet_confs(const int idx) {
		static const std::array<std::vector<std::vector<Vector4i>>, 64> table= {{

			{

				{
					Vector4i(0, 1, 2, 3)
				}
			},
			{

				{
					Vector4i(4, 3, 2, 0),Vector4i(4, 1, 2, 3)
				}
			},
			{

				{
					Vector4i(4, 3, 0, 1),Vector4i(4, 2, 0, 3)
				}
			},
			{

				{
					Vector4i(1, 3, 5, 4),Vector4i(4, 0, 3, 5),Vector4i(5, 0, 3, 2)
				},
				{
					Vector4i(1, 3, 5, 4),Vector4i(4, 2, 3, 5),Vector4i(4, 0, 3, 2)
				}
			},
			{

				{
					Vector4i(4, 3, 1, 2),Vector4i(4, 0, 1, 3)
				}
			},
			{

				{
					Vector4i(0, 4, 5, 3),Vector4i(4, 5, 3, 1),Vector4i(5, 2, 3, 1)
				},
				{
					Vector4i(0, 4, 5, 3),Vector4i(4, 5, 3, 2),Vector4i(4, 2, 3, 1)
				}
			},
			{

				{
					Vector4i(2, 3, 5, 4),Vector4i(4, 1, 3, 5),Vector4i(5, 1, 3, 0)
				},
				{
					Vector4i(2, 3, 5, 4),Vector4i(4, 0, 3, 5),Vector4i(4, 1, 3, 0)
				}
			},
			{},
			{

				{
					Vector4i(4, 2, 1, 0),Vector4i(4, 3, 1, 2)
				}
			},
			{

				{
					Vector4i(0, 2, 5, 4),Vector4i(4, 1, 2, 5),Vector4i(5, 1, 2, 3)
				},
				{
					Vector4i(0, 2, 5, 4),Vector4i(4, 3, 2, 5),Vector4i(4, 1, 2, 3)
				}
			},
			{

				{
					Vector4i(4, 5, 0, 1),Vector4i(4, 2, 0, 5),Vector4i(4, 1, 3, 5),Vector4i(4, 5, 3, 2)
				}
			},
			{

				{
					Vector4i(1, 7, 5, 4),Vector4i(0, 7, 6, 4),Vector4i(2, 3, 6, 7),Vector4i(2, 7, 6, 0),Vector4i(2, 7, 5, 3),Vector4i(3, 7, 5, 1),Vector4i(0, 4, 5, 7),Vector4i(0, 5, 2, 7),Vector4i(1, 4, 6, 7),Vector4i(1, 6, 3, 7)
				},
				{
					Vector4i(1, 7, 5, 4),Vector4i(0, 7, 6, 4),Vector4i(2, 3, 6, 7),Vector4i(2, 7, 6, 0),Vector4i(2, 7, 5, 3),Vector4i(3, 7, 5, 1),Vector4i(0, 4, 2, 7),Vector4i(4, 5, 2, 7),Vector4i(1, 4, 6, 7),Vector4i(1, 6, 3, 7)
				},
				{
					Vector4i(1, 7, 5, 4),Vector4i(0, 7, 6, 4),Vector4i(2, 3, 6, 7),Vector4i(2, 7, 6, 0),Vector4i(2, 7, 5, 3),Vector4i(3, 7, 5, 1),Vector4i(0, 4, 5, 7),Vector4i(0, 5, 2, 7),Vector4i(1, 4, 3, 7),Vector4i(4, 6, 3, 7)
				},
				{
					Vector4i(1, 7, 5, 4),Vector4i(0, 7, 6, 4),Vector4i(2, 3, 6, 7),Vector4i(2, 7, 6, 0),Vector4i(2, 7, 5, 3),Vector4i(3, 7, 5, 1),Vector4i(0, 4, 2, 7),Vector4i(4, 5, 2, 7),Vector4i(1, 4, 3, 7),Vector4i(4, 6, 3, 7)
				}
			},
			{

				{
					Vector4i(0, 4, 5, 1),Vector4i(4, 5, 1, 2),Vector4i(5, 3, 1, 2)
				},
				{
					Vector4i(0, 4, 5, 1),Vector4i(4, 5, 1, 3),Vector4i(4, 3, 1, 2)
				}
			},
			{

				{
					Vector4i(4, 0, 6, 5),Vector4i(4, 5, 6, 2),Vector4i(6, 2, 3, 1),Vector4i(4, 2, 6, 1)
				},
				{
					Vector4i(4, 0, 6, 5),Vector4i(4, 5, 6, 2),Vector4i(6, 2, 3, 4),Vector4i(4, 2, 3, 1)
				},
				{
					Vector4i(4, 0, 6, 5),Vector4i(1, 5, 6, 2),Vector4i(6, 2, 3, 1),Vector4i(4, 5, 6, 1)
				},
				{
					Vector4i(4, 0, 6, 5),Vector4i(4, 5, 3, 2),Vector4i(6, 5, 3, 4),Vector4i(4, 2, 3, 1)
				},
				{
					Vector4i(4, 0, 6, 5),Vector4i(1, 5, 3, 2),Vector4i(1, 5, 6, 3),Vector4i(4, 5, 6, 1)
				},
				{
					Vector4i(4, 0, 6, 5),Vector4i(1, 5, 3, 2),Vector4i(4, 5, 6, 3),Vector4i(4, 5, 3, 1)
				},
				{
					Vector4i(4, 0, 6, 5),Vector4i(1, 4, 3, 7),Vector4i(3, 7, 4, 6),Vector4i(5, 7, 2, 6),Vector4i(3, 6, 2, 7),Vector4i(1, 7, 2, 5),Vector4i(4, 7, 1, 5),Vector4i(1, 7, 3, 2),Vector4i(4, 5, 6, 7)
				},
				{
					Vector4i(4, 0, 6, 5),Vector4i(1, 4, 6, 7),Vector4i(1, 6, 3, 7),Vector4i(5, 7, 3, 6),Vector4i(3, 5, 2, 7),Vector4i(4, 7, 1, 2),Vector4i(4, 7, 2, 5),Vector4i(1, 7, 3, 2),Vector4i(4, 5, 6, 7)
				}
			},
			{

				{
					Vector4i(2, 7, 5, 4),Vector4i(0, 5, 6, 7),Vector4i(1, 7, 6, 3),Vector4i(1, 0, 6, 7),Vector4i(1, 3, 4, 7),Vector4i(3, 2, 4, 7),Vector4i(1, 4, 5, 7),Vector4i(1, 5, 0, 7),Vector4i(2, 7, 6, 5),Vector4i(2, 7, 3, 6)
				},
				{
					Vector4i(2, 7, 5, 4),Vector4i(0, 5, 6, 7),Vector4i(1, 7, 6, 3),Vector4i(1, 0, 6, 7),Vector4i(1, 3, 4, 7),Vector4i(3, 2, 4, 7),Vector4i(1, 4, 0, 7),Vector4i(4, 5, 0, 7),Vector4i(2, 7, 6, 5),Vector4i(2, 7, 3, 6)
				},
				{
					Vector4i(2, 7, 5, 4),Vector4i(0, 5, 6, 7),Vector4i(1, 7, 6, 3),Vector4i(1, 0, 6, 7),Vector4i(1, 3, 4, 7),Vector4i(3, 2, 4, 7),Vector4i(1, 4, 5, 7),Vector4i(1, 5, 0, 7),Vector4i(2, 7, 3, 5),Vector4i(5, 7, 3, 6)
				},
				{
					Vector4i(2, 7, 5, 4),Vector4i(0, 5, 6, 7),Vector4i(1, 7, 6, 3),Vector4i(1, 0, 6, 7),Vector4i(1, 3, 4, 7),Vector4i(3, 2, 4, 7),Vector4i(1, 4, 0, 7),Vector4i(4, 5, 0, 7),Vector4i(2, 7, 3, 5),Vector4i(5, 7, 3, 6)
				}
			},
			{},
			{

				{
					Vector4i(4, 0, 2, 1),Vector4i(4, 3, 2, 0)
				}
			},
			{

				{
					Vector4i(1, 4, 5, 2),Vector4i(4, 5, 2, 0),Vector4i(5, 3, 2, 0)
				},
				{
					Vector4i(1, 4, 5, 2),Vector4i(4, 5, 2, 3),Vector4i(4, 3, 2, 0)
				}
			},
			{

				{
					Vector4i(1, 0, 5, 4),Vector4i(4, 2, 0, 5),Vector4i(5, 2, 0, 3)
				},
				{
					Vector4i(1, 0, 5, 4),Vector4i(4, 3, 0, 5),Vector4i(4, 2, 0, 3)
				}
			},
			{

				{
					Vector4i(4, 5, 6, 1),Vector4i(4, 2, 6, 5),Vector4i(6, 0, 3, 2),Vector4i(4, 0, 6, 2)
				},
				{
					Vector4i(4, 5, 6, 1),Vector4i(4, 2, 6, 5),Vector4i(6, 4, 3, 2),Vector4i(4, 0, 3, 2)
				},
				{
					Vector4i(4, 5, 6, 1),Vector4i(0, 2, 6, 5),Vector4i(6, 0, 3, 2),Vector4i(4, 0, 6, 5)
				},
				{
					Vector4i(4, 5, 6, 1),Vector4i(4, 2, 3, 5),Vector4i(6, 4, 3, 5),Vector4i(4, 0, 3, 2)
				},
				{
					Vector4i(4, 5, 6, 1),Vector4i(0, 2, 3, 5),Vector4i(0, 3, 6, 5),Vector4i(4, 0, 6, 5)
				},
				{
					Vector4i(4, 5, 6, 1),Vector4i(0, 2, 3, 5),Vector4i(4, 3, 6, 5),Vector4i(4, 0, 3, 5)
				},
				{
					Vector4i(4, 5, 6, 1),Vector4i(0, 7, 3, 4),Vector4i(3, 6, 4, 7),Vector4i(5, 6, 2, 7),Vector4i(3, 7, 2, 6),Vector4i(0, 5, 2, 7),Vector4i(4, 5, 0, 7),Vector4i(0, 2, 3, 7),Vector4i(4, 7, 6, 5)
				},
				{
					Vector4i(4, 5, 6, 1),Vector4i(0, 7, 6, 4),Vector4i(0, 7, 3, 6),Vector4i(5, 6, 3, 7),Vector4i(3, 7, 2, 5),Vector4i(4, 2, 0, 7),Vector4i(4, 5, 2, 7),Vector4i(0, 2, 3, 7),Vector4i(4, 7, 6, 5)
				}
			},
			{

				{
					Vector4i(4, 5, 1, 2),Vector4i(4, 0, 1, 5),Vector4i(4, 2, 3, 5),Vector4i(4, 5, 3, 0)
				}
			},
			{

				{
					Vector4i(0, 4, 5, 7),Vector4i(1, 4, 6, 7),Vector4i(2, 3, 5, 7),Vector4i(3, 0, 5, 7),Vector4i(2, 7, 6, 3),Vector4i(2, 1, 6, 7),Vector4i(1, 7, 5, 4),Vector4i(1, 7, 2, 5),Vector4i(0, 7, 6, 4),Vector4i(0, 7, 3, 6)
				},
				{
					Vector4i(0, 4, 5, 7),Vector4i(1, 4, 6, 7),Vector4i(2, 3, 5, 7),Vector4i(3, 0, 5, 7),Vector4i(2, 7, 6, 3),Vector4i(2, 1, 6, 7),Vector4i(1, 7, 2, 4),Vector4i(4, 7, 2, 5),Vector4i(0, 7, 6, 4),Vector4i(0, 7, 3, 6)
				},
				{
					Vector4i(0, 4, 5, 7),Vector4i(1, 4, 6, 7),Vector4i(2, 3, 5, 7),Vector4i(3, 0, 5, 7),Vector4i(2, 7, 6, 3),Vector4i(2, 1, 6, 7),Vector4i(1, 7, 5, 4),Vector4i(1, 7, 2, 5),Vector4i(0, 7, 3, 4),Vector4i(4, 7, 3, 6)
				},
				{
					Vector4i(0, 4, 5, 7),Vector4i(1, 4, 6, 7),Vector4i(2, 3, 5, 7),Vector4i(3, 0, 5, 7),Vector4i(2, 7, 6, 3),Vector4i(2, 1, 6, 7),Vector4i(1, 7, 2, 4),Vector4i(4, 7, 2, 5),Vector4i(0, 7, 3, 4),Vector4i(4, 7, 3, 6)
				}
			},
			{

				{
					Vector4i(2, 7, 5, 4),Vector4i(1, 7, 6, 4),Vector4i(0, 3, 6, 7),Vector4i(0, 7, 6, 1),Vector4i(0, 7, 5, 3),Vector4i(3, 7, 5, 2),Vector4i(1, 4, 5, 7),Vector4i(1, 5, 0, 7),Vector4i(2, 4, 6, 7),Vector4i(2, 6, 3, 7)
				},
				{
					Vector4i(2, 7, 5, 4),Vector4i(1, 7, 6, 4),Vector4i(0, 3, 6, 7),Vector4i(0, 7, 6, 1),Vector4i(0, 7, 5, 3),Vector4i(3, 7, 5, 2),Vector4i(1, 4, 0, 7),Vector4i(4, 5, 0, 7),Vector4i(2, 4, 6, 7),Vector4i(2, 6, 3, 7)
				},
				{
					Vector4i(2, 7, 5, 4),Vector4i(1, 7, 6, 4),Vector4i(0, 3, 6, 7),Vector4i(0, 7, 6, 1),Vector4i(0, 7, 5, 3),Vector4i(3, 7, 5, 2),Vector4i(1, 4, 5, 7),Vector4i(1, 5, 0, 7),Vector4i(2, 4, 3, 7),Vector4i(4, 6, 3, 7)
				},
				{
					Vector4i(2, 7, 5, 4),Vector4i(1, 7, 6, 4),Vector4i(0, 3, 6, 7),Vector4i(0, 7, 6, 1),Vector4i(0, 7, 5, 3),Vector4i(3, 7, 5, 2),Vector4i(1, 4, 0, 7),Vector4i(4, 5, 0, 7),Vector4i(2, 4, 3, 7),Vector4i(4, 6, 3, 7)
				}
			},
			{},
			{

				{
					Vector4i(3, 2, 5, 4),Vector4i(4, 0, 2, 5),Vector4i(5, 0, 2, 1)
				},
				{
					Vector4i(3, 2, 5, 4),Vector4i(4, 1, 2, 5),Vector4i(4, 0, 2, 1)
				}
			},
			{},
			{

				{
					Vector4i(1, 7, 6, 4),Vector4i(3, 7, 6, 5),Vector4i(0, 7, 4, 2),Vector4i(0, 1, 4, 7),Vector4i(0, 2, 5, 7),Vector4i(2, 3, 5, 7),Vector4i(2, 4, 6, 7),Vector4i(2, 6, 3, 7),Vector4i(0, 5, 6, 7),Vector4i(0, 6, 1, 7)
				},
				{
					Vector4i(1, 7, 6, 4),Vector4i(3, 7, 6, 5),Vector4i(0, 7, 4, 2),Vector4i(0, 1, 4, 7),Vector4i(0, 2, 5, 7),Vector4i(2, 3, 5, 7),Vector4i(2, 4, 3, 7),Vector4i(4, 6, 3, 7),Vector4i(0, 5, 6, 7),Vector4i(0, 6, 1, 7)
				},
				{
					Vector4i(1, 7, 6, 4),Vector4i(3, 7, 6, 5),Vector4i(0, 7, 4, 2),Vector4i(0, 1, 4, 7),Vector4i(0, 2, 5, 7),Vector4i(2, 3, 5, 7),Vector4i(2, 4, 6, 7),Vector4i(2, 6, 3, 7),Vector4i(0, 5, 1, 7),Vector4i(5, 6, 1, 7)
				},
				{
					Vector4i(1, 7, 6, 4),Vector4i(3, 7, 6, 5),Vector4i(0, 7, 4, 2),Vector4i(0, 1, 4, 7),Vector4i(0, 2, 5, 7),Vector4i(2, 3, 5, 7),Vector4i(2, 4, 3, 7),Vector4i(4, 6, 3, 7),Vector4i(0, 5, 1, 7),Vector4i(5, 6, 1, 7)
				}
			},
			{},
			{

				{
					Vector4i(0, 4, 5, 7),Vector4i(3, 7, 6, 5),Vector4i(1, 2, 4, 7),Vector4i(1, 7, 4, 0),Vector4i(1, 7, 6, 2),Vector4i(2, 7, 6, 3),Vector4i(2, 7, 5, 4),Vector4i(2, 7, 3, 5),Vector4i(0, 5, 6, 7),Vector4i(0, 6, 1, 7)
				},
				{
					Vector4i(0, 4, 5, 7),Vector4i(3, 7, 6, 5),Vector4i(1, 2, 4, 7),Vector4i(1, 7, 4, 0),Vector4i(1, 7, 6, 2),Vector4i(2, 7, 6, 3),Vector4i(2, 7, 3, 4),Vector4i(4, 7, 3, 5),Vector4i(0, 5, 6, 7),Vector4i(0, 6, 1, 7)
				},
				{
					Vector4i(0, 4, 5, 7),Vector4i(3, 7, 6, 5),Vector4i(1, 2, 4, 7),Vector4i(1, 7, 4, 0),Vector4i(1, 7, 6, 2),Vector4i(2, 7, 6, 3),Vector4i(2, 7, 5, 4),Vector4i(2, 7, 3, 5),Vector4i(0, 5, 1, 7),Vector4i(5, 6, 1, 7)
				},
				{
					Vector4i(0, 4, 5, 7),Vector4i(3, 7, 6, 5),Vector4i(1, 2, 4, 7),Vector4i(1, 7, 4, 0),Vector4i(1, 7, 6, 2),Vector4i(2, 7, 6, 3),Vector4i(2, 7, 3, 4),Vector4i(4, 7, 3, 5),Vector4i(0, 5, 1, 7),Vector4i(5, 6, 1, 7)
				}
			},
			{},
			{

				{
					Vector4i(1, 5, 7, 4),Vector4i(7, 0, 6, 5),Vector4i(1, 0, 7, 5),Vector4i(3, 5, 7, 6),Vector4i(7, 2, 4, 5),Vector4i(3, 2, 7, 5)
				},
				{
					Vector4i(1, 5, 7, 4),Vector4i(7, 0, 6, 5),Vector4i(1, 0, 7, 5),Vector4i(3, 5, 7, 6),Vector4i(7, 3, 4, 5),Vector4i(3, 2, 4, 5)
				},
				{
					Vector4i(1, 5, 7, 4),Vector4i(7, 0, 6, 5),Vector4i(1, 0, 7, 5),Vector4i(2, 5, 7, 6),Vector4i(7, 2, 4, 5),Vector4i(3, 2, 7, 6)
				},
				{
					Vector4i(1, 5, 7, 4),Vector4i(7, 0, 6, 5),Vector4i(1, 0, 7, 5),Vector4i(2, 8, 4, 3),Vector4i(4, 7, 3, 8),Vector4i(6, 7, 5, 8),Vector4i(4, 8, 5, 7),Vector4i(2, 6, 5, 8),Vector4i(3, 6, 2, 8),Vector4i(2, 5, 4, 8),Vector4i(3, 8, 7, 6)
				},
				{
					Vector4i(1, 5, 7, 4),Vector4i(7, 1, 6, 5),Vector4i(1, 0, 6, 5),Vector4i(3, 5, 7, 6),Vector4i(7, 2, 4, 5),Vector4i(3, 2, 7, 5)
				},
				{
					Vector4i(1, 5, 7, 4),Vector4i(7, 1, 6, 5),Vector4i(1, 0, 6, 5),Vector4i(3, 5, 7, 6),Vector4i(7, 3, 4, 5),Vector4i(3, 2, 4, 5)
				},
				{
					Vector4i(1, 5, 7, 4),Vector4i(7, 1, 6, 5),Vector4i(1, 0, 6, 5),Vector4i(2, 5, 7, 6),Vector4i(7, 2, 4, 5),Vector4i(3, 2, 7, 6)
				},
				{
					Vector4i(1, 5, 7, 4),Vector4i(7, 1, 6, 5),Vector4i(1, 0, 6, 5),Vector4i(2, 8, 4, 3),Vector4i(4, 7, 3, 8),Vector4i(6, 7, 5, 8),Vector4i(4, 8, 5, 7),Vector4i(2, 6, 5, 8),Vector4i(3, 6, 2, 8),Vector4i(2, 5, 4, 8),Vector4i(3, 8, 7, 6)
				},
				{
					Vector4i(0, 5, 7, 4),Vector4i(7, 0, 6, 5),Vector4i(1, 0, 7, 4),Vector4i(3, 5, 7, 6),Vector4i(7, 2, 4, 5),Vector4i(3, 2, 7, 5)
				},
				{
					Vector4i(0, 5, 7, 4),Vector4i(7, 0, 6, 5),Vector4i(1, 0, 7, 4),Vector4i(3, 5, 7, 6),Vector4i(7, 3, 4, 5),Vector4i(3, 2, 4, 5)
				},
				{
					Vector4i(0, 5, 7, 4),Vector4i(7, 0, 6, 5),Vector4i(1, 0, 7, 4),Vector4i(2, 5, 7, 6),Vector4i(7, 2, 4, 5),Vector4i(3, 2, 7, 6)
				},
				{
					Vector4i(0, 5, 7, 4),Vector4i(7, 0, 6, 5),Vector4i(1, 0, 7, 4),Vector4i(2, 8, 4, 3),Vector4i(4, 7, 3, 8),Vector4i(6, 7, 5, 8),Vector4i(4, 8, 5, 7),Vector4i(2, 6, 5, 8),Vector4i(3, 6, 2, 8),Vector4i(2, 5, 4, 8),Vector4i(3, 8, 7, 6)
				},
				{
					Vector4i(1, 5, 6, 4),Vector4i(7, 1, 6, 4),Vector4i(1, 0, 6, 5),Vector4i(3, 5, 4, 6),Vector4i(7, 3, 4, 6),Vector4i(3, 2, 4, 5)
				},
				{
					Vector4i(1, 5, 6, 4),Vector4i(7, 1, 6, 4),Vector4i(1, 0, 6, 5),Vector4i(2, 5, 4, 6),Vector4i(2, 4, 7, 6),Vector4i(3, 2, 7, 6)
				},
				{
					Vector4i(1, 5, 6, 4),Vector4i(7, 1, 6, 4),Vector4i(1, 0, 6, 5),Vector4i(2, 5, 4, 6),Vector4i(3, 4, 7, 6),Vector4i(3, 2, 4, 6)
				},
				{
					Vector4i(1, 5, 6, 4),Vector4i(7, 1, 6, 4),Vector4i(1, 0, 6, 5),Vector4i(2, 8, 7, 3),Vector4i(2, 8, 4, 7),Vector4i(6, 7, 4, 8),Vector4i(4, 8, 5, 6),Vector4i(3, 5, 2, 8),Vector4i(3, 6, 5, 8),Vector4i(2, 5, 4, 8),Vector4i(3, 8, 7, 6)
				},
				{
					Vector4i(0, 5, 6, 4),Vector4i(0, 6, 7, 4),Vector4i(1, 0, 7, 4),Vector4i(3, 5, 4, 6),Vector4i(7, 3, 4, 6),Vector4i(3, 2, 4, 5)
				},
				{
					Vector4i(0, 5, 6, 4),Vector4i(0, 6, 7, 4),Vector4i(1, 0, 7, 4),Vector4i(2, 5, 4, 6),Vector4i(2, 4, 7, 6),Vector4i(3, 2, 7, 6)
				},
				{
					Vector4i(0, 5, 6, 4),Vector4i(0, 6, 7, 4),Vector4i(1, 0, 7, 4),Vector4i(2, 5, 4, 6),Vector4i(3, 4, 7, 6),Vector4i(3, 2, 4, 6)
				},
				{
					Vector4i(0, 5, 6, 4),Vector4i(0, 6, 7, 4),Vector4i(1, 0, 7, 4),Vector4i(2, 8, 7, 3),Vector4i(2, 8, 4, 7),Vector4i(6, 7, 4, 8),Vector4i(4, 8, 5, 6),Vector4i(3, 5, 2, 8),Vector4i(3, 6, 5, 8),Vector4i(2, 5, 4, 8),Vector4i(3, 8, 7, 6)
				},
				{
					Vector4i(0, 5, 6, 4),Vector4i(1, 6, 7, 4),Vector4i(1, 0, 6, 4),Vector4i(3, 5, 4, 6),Vector4i(7, 3, 4, 6),Vector4i(3, 2, 4, 5)
				},
				{
					Vector4i(0, 5, 6, 4),Vector4i(1, 6, 7, 4),Vector4i(1, 0, 6, 4),Vector4i(2, 5, 4, 6),Vector4i(2, 4, 7, 6),Vector4i(3, 2, 7, 6)
				},
				{
					Vector4i(0, 5, 6, 4),Vector4i(1, 6, 7, 4),Vector4i(1, 0, 6, 4),Vector4i(2, 5, 4, 6),Vector4i(3, 4, 7, 6),Vector4i(3, 2, 4, 6)
				},
				{
					Vector4i(0, 5, 6, 4),Vector4i(1, 6, 7, 4),Vector4i(1, 0, 6, 4),Vector4i(2, 8, 7, 3),Vector4i(2, 8, 4, 7),Vector4i(6, 7, 4, 8),Vector4i(4, 8, 5, 6),Vector4i(3, 5, 2, 8),Vector4i(3, 6, 5, 8),Vector4i(2, 5, 4, 8),Vector4i(3, 8, 7, 6)
				},
				{
					Vector4i(0, 8, 6, 1),Vector4i(6, 7, 1, 8),Vector4i(4, 7, 5, 8),Vector4i(6, 8, 5, 7),Vector4i(0, 4, 5, 8),Vector4i(1, 4, 0, 8),Vector4i(0, 5, 6, 8),Vector4i(1, 8, 7, 4),Vector4i(3, 5, 7, 6),Vector4i(7, 2, 4, 5),Vector4i(3, 2, 7, 5)
				},
				{
					Vector4i(0, 8, 6, 1),Vector4i(6, 7, 1, 8),Vector4i(4, 7, 5, 8),Vector4i(6, 8, 5, 7),Vector4i(0, 4, 5, 8),Vector4i(1, 4, 0, 8),Vector4i(0, 5, 6, 8),Vector4i(1, 8, 7, 4),Vector4i(3, 5, 7, 6),Vector4i(7, 3, 4, 5),Vector4i(3, 2, 4, 5)
				},
				{
					Vector4i(0, 8, 6, 1),Vector4i(6, 7, 1, 8),Vector4i(4, 7, 5, 8),Vector4i(6, 8, 5, 7),Vector4i(0, 4, 5, 8),Vector4i(1, 4, 0, 8),Vector4i(0, 5, 6, 8),Vector4i(1, 8, 7, 4),Vector4i(2, 5, 7, 6),Vector4i(7, 2, 4, 5),Vector4i(3, 2, 7, 6)
				},
				{
					Vector4i(0, 8, 6, 1),Vector4i(6, 7, 1, 8),Vector4i(4, 7, 5, 8),Vector4i(6, 8, 5, 7),Vector4i(0, 4, 5, 8),Vector4i(1, 4, 0, 8),Vector4i(0, 5, 6, 8),Vector4i(1, 8, 7, 4),Vector4i(2, 9, 4, 3),Vector4i(4, 7, 3, 9),Vector4i(6, 7, 5, 9),Vector4i(4, 9, 5, 7),Vector4i(2, 6, 5, 9),Vector4i(3, 6, 2, 9),Vector4i(2, 5, 4, 9),Vector4i(3, 9, 7, 6)
				},
				{
					Vector4i(0, 8, 7, 1),Vector4i(0, 8, 6, 7),Vector4i(4, 7, 6, 8),Vector4i(6, 8, 5, 4),Vector4i(1, 5, 0, 8),Vector4i(1, 4, 5, 8),Vector4i(0, 5, 6, 8),Vector4i(1, 8, 7, 4),Vector4i(3, 5, 4, 6),Vector4i(7, 3, 4, 6),Vector4i(3, 2, 4, 5)
				},
				{
					Vector4i(0, 8, 7, 1),Vector4i(0, 8, 6, 7),Vector4i(4, 7, 6, 8),Vector4i(6, 8, 5, 4),Vector4i(1, 5, 0, 8),Vector4i(1, 4, 5, 8),Vector4i(0, 5, 6, 8),Vector4i(1, 8, 7, 4),Vector4i(2, 5, 4, 6),Vector4i(2, 4, 7, 6),Vector4i(3, 2, 7, 6)
				},
				{
					Vector4i(0, 8, 7, 1),Vector4i(0, 8, 6, 7),Vector4i(4, 7, 6, 8),Vector4i(6, 8, 5, 4),Vector4i(1, 5, 0, 8),Vector4i(1, 4, 5, 8),Vector4i(0, 5, 6, 8),Vector4i(1, 8, 7, 4),Vector4i(2, 5, 4, 6),Vector4i(3, 4, 7, 6),Vector4i(3, 2, 4, 6)
				},
				{
					Vector4i(0, 8, 7, 1),Vector4i(0, 8, 6, 7),Vector4i(4, 7, 6, 8),Vector4i(6, 8, 5, 4),Vector4i(1, 5, 0, 8),Vector4i(1, 4, 5, 8),Vector4i(0, 5, 6, 8),Vector4i(1, 8, 7, 4),Vector4i(2, 9, 7, 3),Vector4i(2, 9, 4, 7),Vector4i(6, 7, 4, 9),Vector4i(4, 9, 5, 6),Vector4i(3, 5, 2, 9),Vector4i(3, 6, 5, 9),Vector4i(2, 5, 4, 9),Vector4i(3, 9, 7, 6)
				}
			},
			{},
			{

				{
					Vector4i(4, 1, 0, 2),Vector4i(4, 3, 0, 1)
				}
			},
			{

				{
					Vector4i(4, 5, 2, 0),Vector4i(4, 1, 2, 5),Vector4i(4, 0, 3, 5),Vector4i(4, 5, 3, 1)
				}
			},
			{

				{
					Vector4i(2, 4, 5, 0),Vector4i(4, 5, 0, 1),Vector4i(5, 3, 0, 1)
				},
				{
					Vector4i(2, 4, 5, 0),Vector4i(4, 5, 0, 3),Vector4i(4, 3, 0, 1)
				}
			},
			{

				{
					Vector4i(1, 7, 5, 4),Vector4i(2, 5, 6, 7),Vector4i(0, 3, 4, 7),Vector4i(3, 1, 4, 7),Vector4i(0, 7, 6, 3),Vector4i(0, 2, 6, 7),Vector4i(0, 4, 5, 7),Vector4i(0, 5, 2, 7),Vector4i(1, 7, 6, 5),Vector4i(1, 7, 3, 6)
				},
				{
					Vector4i(1, 7, 5, 4),Vector4i(2, 5, 6, 7),Vector4i(0, 3, 4, 7),Vector4i(3, 1, 4, 7),Vector4i(0, 7, 6, 3),Vector4i(0, 2, 6, 7),Vector4i(0, 4, 2, 7),Vector4i(4, 5, 2, 7),Vector4i(1, 7, 6, 5),Vector4i(1, 7, 3, 6)
				},
				{
					Vector4i(1, 7, 5, 4),Vector4i(2, 5, 6, 7),Vector4i(0, 3, 4, 7),Vector4i(3, 1, 4, 7),Vector4i(0, 7, 6, 3),Vector4i(0, 2, 6, 7),Vector4i(0, 4, 5, 7),Vector4i(0, 5, 2, 7),Vector4i(1, 7, 3, 5),Vector4i(5, 7, 3, 6)
				},
				{
					Vector4i(1, 7, 5, 4),Vector4i(2, 5, 6, 7),Vector4i(0, 3, 4, 7),Vector4i(3, 1, 4, 7),Vector4i(0, 7, 6, 3),Vector4i(0, 2, 6, 7),Vector4i(0, 4, 2, 7),Vector4i(4, 5, 2, 7),Vector4i(1, 7, 3, 5),Vector4i(5, 7, 3, 6)
				}
			},
			{

				{
					Vector4i(2, 1, 5, 4),Vector4i(4, 0, 1, 5),Vector4i(5, 0, 1, 3)
				},
				{
					Vector4i(2, 1, 5, 4),Vector4i(4, 3, 1, 5),Vector4i(4, 0, 1, 3)
				}
			},
			{

				{
					Vector4i(0, 4, 5, 7),Vector4i(2, 7, 6, 5),Vector4i(1, 7, 4, 3),Vector4i(3, 7, 4, 0),Vector4i(1, 3, 6, 7),Vector4i(1, 7, 6, 2),Vector4i(1, 7, 5, 4),Vector4i(1, 7, 2, 5),Vector4i(0, 5, 6, 7),Vector4i(0, 6, 3, 7)
				},
				{
					Vector4i(0, 4, 5, 7),Vector4i(2, 7, 6, 5),Vector4i(1, 7, 4, 3),Vector4i(3, 7, 4, 0),Vector4i(1, 3, 6, 7),Vector4i(1, 7, 6, 2),Vector4i(1, 7, 2, 4),Vector4i(4, 7, 2, 5),Vector4i(0, 5, 6, 7),Vector4i(0, 6, 3, 7)
				},
				{
					Vector4i(0, 4, 5, 7),Vector4i(2, 7, 6, 5),Vector4i(1, 7, 4, 3),Vector4i(3, 7, 4, 0),Vector4i(1, 3, 6, 7),Vector4i(1, 7, 6, 2),Vector4i(1, 7, 5, 4),Vector4i(1, 7, 2, 5),Vector4i(0, 5, 3, 7),Vector4i(5, 6, 3, 7)
				},
				{
					Vector4i(0, 4, 5, 7),Vector4i(2, 7, 6, 5),Vector4i(1, 7, 4, 3),Vector4i(3, 7, 4, 0),Vector4i(1, 3, 6, 7),Vector4i(1, 7, 6, 2),Vector4i(1, 7, 2, 4),Vector4i(4, 7, 2, 5),Vector4i(0, 5, 3, 7),Vector4i(5, 6, 3, 7)
				}
			},
			{

				{
					Vector4i(4, 5, 6, 2),Vector4i(4, 0, 6, 5),Vector4i(6, 1, 3, 0),Vector4i(4, 1, 6, 0)
				},
				{
					Vector4i(4, 5, 6, 2),Vector4i(4, 0, 6, 5),Vector4i(6, 4, 3, 0),Vector4i(4, 1, 3, 0)
				},
				{
					Vector4i(4, 5, 6, 2),Vector4i(1, 0, 6, 5),Vector4i(6, 1, 3, 0),Vector4i(4, 1, 6, 5)
				},
				{
					Vector4i(4, 5, 6, 2),Vector4i(4, 0, 3, 5),Vector4i(6, 4, 3, 5),Vector4i(4, 1, 3, 0)
				},
				{
					Vector4i(4, 5, 6, 2),Vector4i(1, 0, 3, 5),Vector4i(1, 3, 6, 5),Vector4i(4, 1, 6, 5)
				},
				{
					Vector4i(4, 5, 6, 2),Vector4i(1, 0, 3, 5),Vector4i(4, 3, 6, 5),Vector4i(4, 1, 3, 5)
				},
				{
					Vector4i(4, 5, 6, 2),Vector4i(1, 7, 3, 4),Vector4i(3, 6, 4, 7),Vector4i(5, 6, 0, 7),Vector4i(3, 7, 0, 6),Vector4i(1, 5, 0, 7),Vector4i(4, 5, 1, 7),Vector4i(1, 0, 3, 7),Vector4i(4, 7, 6, 5)
				},
				{
					Vector4i(4, 5, 6, 2),Vector4i(1, 7, 6, 4),Vector4i(1, 7, 3, 6),Vector4i(5, 6, 3, 7),Vector4i(3, 7, 0, 5),Vector4i(4, 0, 1, 7),Vector4i(4, 5, 0, 7),Vector4i(1, 0, 3, 7),Vector4i(4, 7, 6, 5)
				}
			},
			{},
			{

				{
					Vector4i(3, 4, 5, 1),Vector4i(4, 5, 1, 0),Vector4i(5, 2, 1, 0)
				},
				{
					Vector4i(3, 4, 5, 1),Vector4i(4, 5, 1, 2),Vector4i(4, 2, 1, 0)
				}
			},
			{

				{
					Vector4i(0, 7, 5, 4),Vector4i(3, 5, 6, 7),Vector4i(1, 2, 4, 7),Vector4i(2, 0, 4, 7),Vector4i(1, 7, 6, 2),Vector4i(1, 3, 6, 7),Vector4i(1, 4, 5, 7),Vector4i(1, 5, 3, 7),Vector4i(0, 7, 6, 5),Vector4i(0, 7, 2, 6)
				},
				{
					Vector4i(0, 7, 5, 4),Vector4i(3, 5, 6, 7),Vector4i(1, 2, 4, 7),Vector4i(2, 0, 4, 7),Vector4i(1, 7, 6, 2),Vector4i(1, 3, 6, 7),Vector4i(1, 4, 3, 7),Vector4i(4, 5, 3, 7),Vector4i(0, 7, 6, 5),Vector4i(0, 7, 2, 6)
				},
				{
					Vector4i(0, 7, 5, 4),Vector4i(3, 5, 6, 7),Vector4i(1, 2, 4, 7),Vector4i(2, 0, 4, 7),Vector4i(1, 7, 6, 2),Vector4i(1, 3, 6, 7),Vector4i(1, 4, 5, 7),Vector4i(1, 5, 3, 7),Vector4i(0, 7, 2, 5),Vector4i(5, 7, 2, 6)
				},
				{
					Vector4i(0, 7, 5, 4),Vector4i(3, 5, 6, 7),Vector4i(1, 2, 4, 7),Vector4i(2, 0, 4, 7),Vector4i(1, 7, 6, 2),Vector4i(1, 3, 6, 7),Vector4i(1, 4, 3, 7),Vector4i(4, 5, 3, 7),Vector4i(0, 7, 2, 5),Vector4i(5, 7, 2, 6)
				}
			},
			{

				{
					Vector4i(2, 4, 6, 7),Vector4i(3, 5, 6, 7),Vector4i(0, 1, 4, 7),Vector4i(0, 7, 4, 2),Vector4i(0, 7, 5, 1),Vector4i(1, 7, 5, 3),Vector4i(1, 7, 6, 4),Vector4i(1, 7, 3, 6),Vector4i(0, 7, 6, 5),Vector4i(0, 7, 2, 6)
				},
				{
					Vector4i(2, 4, 6, 7),Vector4i(3, 5, 6, 7),Vector4i(0, 1, 4, 7),Vector4i(0, 7, 4, 2),Vector4i(0, 7, 5, 1),Vector4i(1, 7, 5, 3),Vector4i(1, 7, 3, 4),Vector4i(4, 7, 3, 6),Vector4i(0, 7, 6, 5),Vector4i(0, 7, 2, 6)
				},
				{
					Vector4i(2, 4, 6, 7),Vector4i(3, 5, 6, 7),Vector4i(0, 1, 4, 7),Vector4i(0, 7, 4, 2),Vector4i(0, 7, 5, 1),Vector4i(1, 7, 5, 3),Vector4i(1, 7, 6, 4),Vector4i(1, 7, 3, 6),Vector4i(0, 7, 2, 5),Vector4i(5, 7, 2, 6)
				},
				{
					Vector4i(2, 4, 6, 7),Vector4i(3, 5, 6, 7),Vector4i(0, 1, 4, 7),Vector4i(0, 7, 4, 2),Vector4i(0, 7, 5, 1),Vector4i(1, 7, 5, 3),Vector4i(1, 7, 3, 4),Vector4i(4, 7, 3, 6),Vector4i(0, 7, 2, 5),Vector4i(5, 7, 2, 6)
				}
			},
			{

				{
					Vector4i(2, 5, 7, 4),Vector4i(7, 4, 6, 0),Vector4i(2, 4, 7, 0),Vector4i(3, 6, 7, 4),Vector4i(7, 4, 5, 1),Vector4i(3, 4, 7, 1)
				},
				{
					Vector4i(2, 5, 7, 4),Vector4i(7, 4, 6, 0),Vector4i(2, 4, 7, 0),Vector4i(3, 6, 7, 4),Vector4i(7, 4, 5, 3),Vector4i(3, 4, 5, 1)
				},
				{
					Vector4i(2, 5, 7, 4),Vector4i(7, 4, 6, 0),Vector4i(2, 4, 7, 0),Vector4i(1, 6, 7, 4),Vector4i(7, 4, 5, 1),Vector4i(3, 6, 7, 1)
				},
				{
					Vector4i(2, 5, 7, 4),Vector4i(7, 4, 6, 0),Vector4i(2, 4, 7, 0),Vector4i(1, 3, 5, 8),Vector4i(5, 8, 3, 7),Vector4i(6, 8, 4, 7),Vector4i(5, 7, 4, 8),Vector4i(1, 8, 4, 6),Vector4i(3, 8, 1, 6),Vector4i(1, 8, 5, 4),Vector4i(3, 6, 7, 8)
				},
				{
					Vector4i(2, 5, 7, 4),Vector4i(7, 4, 6, 2),Vector4i(2, 4, 6, 0),Vector4i(3, 6, 7, 4),Vector4i(7, 4, 5, 1),Vector4i(3, 4, 7, 1)
				},
				{
					Vector4i(2, 5, 7, 4),Vector4i(7, 4, 6, 2),Vector4i(2, 4, 6, 0),Vector4i(3, 6, 7, 4),Vector4i(7, 4, 5, 3),Vector4i(3, 4, 5, 1)
				},
				{
					Vector4i(2, 5, 7, 4),Vector4i(7, 4, 6, 2),Vector4i(2, 4, 6, 0),Vector4i(1, 6, 7, 4),Vector4i(7, 4, 5, 1),Vector4i(3, 6, 7, 1)
				},
				{
					Vector4i(2, 5, 7, 4),Vector4i(7, 4, 6, 2),Vector4i(2, 4, 6, 0),Vector4i(1, 3, 5, 8),Vector4i(5, 8, 3, 7),Vector4i(6, 8, 4, 7),Vector4i(5, 7, 4, 8),Vector4i(1, 8, 4, 6),Vector4i(3, 8, 1, 6),Vector4i(1, 8, 5, 4),Vector4i(3, 6, 7, 8)
				},
				{
					Vector4i(0, 5, 7, 4),Vector4i(7, 4, 6, 0),Vector4i(2, 5, 7, 0),Vector4i(3, 6, 7, 4),Vector4i(7, 4, 5, 1),Vector4i(3, 4, 7, 1)
				},
				{
					Vector4i(0, 5, 7, 4),Vector4i(7, 4, 6, 0),Vector4i(2, 5, 7, 0),Vector4i(3, 6, 7, 4),Vector4i(7, 4, 5, 3),Vector4i(3, 4, 5, 1)
				},
				{
					Vector4i(0, 5, 7, 4),Vector4i(7, 4, 6, 0),Vector4i(2, 5, 7, 0),Vector4i(1, 6, 7, 4),Vector4i(7, 4, 5, 1),Vector4i(3, 6, 7, 1)
				},
				{
					Vector4i(0, 5, 7, 4),Vector4i(7, 4, 6, 0),Vector4i(2, 5, 7, 0),Vector4i(1, 3, 5, 8),Vector4i(5, 8, 3, 7),Vector4i(6, 8, 4, 7),Vector4i(5, 7, 4, 8),Vector4i(1, 8, 4, 6),Vector4i(3, 8, 1, 6),Vector4i(1, 8, 5, 4),Vector4i(3, 6, 7, 8)
				},
				{
					Vector4i(2, 5, 6, 4),Vector4i(7, 5, 6, 2),Vector4i(2, 4, 6, 0),Vector4i(3, 6, 5, 4),Vector4i(7, 6, 5, 3),Vector4i(3, 4, 5, 1)
				},
				{
					Vector4i(2, 5, 6, 4),Vector4i(7, 5, 6, 2),Vector4i(2, 4, 6, 0),Vector4i(1, 6, 5, 4),Vector4i(1, 6, 7, 5),Vector4i(3, 6, 7, 1)
				},
				{
					Vector4i(2, 5, 6, 4),Vector4i(7, 5, 6, 2),Vector4i(2, 4, 6, 0),Vector4i(1, 6, 5, 4),Vector4i(3, 6, 7, 5),Vector4i(3, 6, 5, 1)
				},
				{
					Vector4i(2, 5, 6, 4),Vector4i(7, 5, 6, 2),Vector4i(2, 4, 6, 0),Vector4i(1, 3, 7, 8),Vector4i(1, 7, 5, 8),Vector4i(6, 8, 5, 7),Vector4i(5, 6, 4, 8),Vector4i(3, 8, 1, 4),Vector4i(3, 8, 4, 6),Vector4i(1, 8, 5, 4),Vector4i(3, 6, 7, 8)
				},
				{
					Vector4i(0, 5, 6, 4),Vector4i(0, 5, 7, 6),Vector4i(2, 5, 7, 0),Vector4i(3, 6, 5, 4),Vector4i(7, 6, 5, 3),Vector4i(3, 4, 5, 1)
				},
				{
					Vector4i(0, 5, 6, 4),Vector4i(0, 5, 7, 6),Vector4i(2, 5, 7, 0),Vector4i(1, 6, 5, 4),Vector4i(1, 6, 7, 5),Vector4i(3, 6, 7, 1)
				},
				{
					Vector4i(0, 5, 6, 4),Vector4i(0, 5, 7, 6),Vector4i(2, 5, 7, 0),Vector4i(1, 6, 5, 4),Vector4i(3, 6, 7, 5),Vector4i(3, 6, 5, 1)
				},
				{
					Vector4i(0, 5, 6, 4),Vector4i(0, 5, 7, 6),Vector4i(2, 5, 7, 0),Vector4i(1, 3, 7, 8),Vector4i(1, 7, 5, 8),Vector4i(6, 8, 5, 7),Vector4i(5, 6, 4, 8),Vector4i(3, 8, 1, 4),Vector4i(3, 8, 4, 6),Vector4i(1, 8, 5, 4),Vector4i(3, 6, 7, 8)
				},
				{
					Vector4i(0, 5, 6, 4),Vector4i(2, 5, 7, 6),Vector4i(2, 5, 6, 0),Vector4i(3, 6, 5, 4),Vector4i(7, 6, 5, 3),Vector4i(3, 4, 5, 1)
				},
				{
					Vector4i(0, 5, 6, 4),Vector4i(2, 5, 7, 6),Vector4i(2, 5, 6, 0),Vector4i(1, 6, 5, 4),Vector4i(1, 6, 7, 5),Vector4i(3, 6, 7, 1)
				},
				{
					Vector4i(0, 5, 6, 4),Vector4i(2, 5, 7, 6),Vector4i(2, 5, 6, 0),Vector4i(1, 6, 5, 4),Vector4i(3, 6, 7, 5),Vector4i(3, 6, 5, 1)
				},
				{
					Vector4i(0, 5, 6, 4),Vector4i(2, 5, 7, 6),Vector4i(2, 5, 6, 0),Vector4i(1, 3, 7, 8),Vector4i(1, 7, 5, 8),Vector4i(6, 8, 5, 7),Vector4i(5, 6, 4, 8),Vector4i(3, 8, 1, 4),Vector4i(3, 8, 4, 6),Vector4i(1, 8, 5, 4),Vector4i(3, 6, 7, 8)
				},
				{
					Vector4i(0, 2, 6, 8),Vector4i(6, 8, 2, 7),Vector4i(5, 8, 4, 7),Vector4i(6, 7, 4, 8),Vector4i(0, 8, 4, 5),Vector4i(2, 8, 0, 5),Vector4i(0, 8, 6, 4),Vector4i(2, 5, 7, 8),Vector4i(3, 6, 7, 4),Vector4i(7, 4, 5, 1),Vector4i(3, 4, 7, 1)
				},
				{
					Vector4i(0, 2, 6, 8),Vector4i(6, 8, 2, 7),Vector4i(5, 8, 4, 7),Vector4i(6, 7, 4, 8),Vector4i(0, 8, 4, 5),Vector4i(2, 8, 0, 5),Vector4i(0, 8, 6, 4),Vector4i(2, 5, 7, 8),Vector4i(3, 6, 7, 4),Vector4i(7, 4, 5, 3),Vector4i(3, 4, 5, 1)
				},
				{
					Vector4i(0, 2, 6, 8),Vector4i(6, 8, 2, 7),Vector4i(5, 8, 4, 7),Vector4i(6, 7, 4, 8),Vector4i(0, 8, 4, 5),Vector4i(2, 8, 0, 5),Vector4i(0, 8, 6, 4),Vector4i(2, 5, 7, 8),Vector4i(1, 6, 7, 4),Vector4i(7, 4, 5, 1),Vector4i(3, 6, 7, 1)
				},
				{
					Vector4i(0, 2, 6, 8),Vector4i(6, 8, 2, 7),Vector4i(5, 8, 4, 7),Vector4i(6, 7, 4, 8),Vector4i(0, 8, 4, 5),Vector4i(2, 8, 0, 5),Vector4i(0, 8, 6, 4),Vector4i(2, 5, 7, 8),Vector4i(1, 3, 5, 9),Vector4i(5, 9, 3, 7),Vector4i(6, 9, 4, 7),Vector4i(5, 7, 4, 9),Vector4i(1, 9, 4, 6),Vector4i(3, 9, 1, 6),Vector4i(1, 9, 5, 4),Vector4i(3, 6, 7, 9)
				},
				{
					Vector4i(0, 2, 7, 8),Vector4i(0, 7, 6, 8),Vector4i(5, 8, 6, 7),Vector4i(6, 5, 4, 8),Vector4i(2, 8, 0, 4),Vector4i(2, 8, 4, 5),Vector4i(0, 8, 6, 4),Vector4i(2, 5, 7, 8),Vector4i(3, 6, 5, 4),Vector4i(7, 6, 5, 3),Vector4i(3, 4, 5, 1)
				},
				{
					Vector4i(0, 2, 7, 8),Vector4i(0, 7, 6, 8),Vector4i(5, 8, 6, 7),Vector4i(6, 5, 4, 8),Vector4i(2, 8, 0, 4),Vector4i(2, 8, 4, 5),Vector4i(0, 8, 6, 4),Vector4i(2, 5, 7, 8),Vector4i(1, 6, 5, 4),Vector4i(1, 6, 7, 5),Vector4i(3, 6, 7, 1)
				},
				{
					Vector4i(0, 2, 7, 8),Vector4i(0, 7, 6, 8),Vector4i(5, 8, 6, 7),Vector4i(6, 5, 4, 8),Vector4i(2, 8, 0, 4),Vector4i(2, 8, 4, 5),Vector4i(0, 8, 6, 4),Vector4i(2, 5, 7, 8),Vector4i(1, 6, 5, 4),Vector4i(3, 6, 7, 5),Vector4i(3, 6, 5, 1)
				},
				{
					Vector4i(0, 2, 7, 8),Vector4i(0, 7, 6, 8),Vector4i(5, 8, 6, 7),Vector4i(6, 5, 4, 8),Vector4i(2, 8, 0, 4),Vector4i(2, 8, 4, 5),Vector4i(0, 8, 6, 4),Vector4i(2, 5, 7, 8),Vector4i(1, 3, 7, 9),Vector4i(1, 7, 5, 9),Vector4i(6, 9, 5, 7),Vector4i(5, 6, 4, 9),Vector4i(3, 9, 1, 4),Vector4i(3, 9, 4, 6),Vector4i(1, 9, 5, 4),Vector4i(3, 6, 7, 9)
				}
			},
			{},
			{},
			{},
			{},
			{

				{
					Vector4i(3, 0, 5, 4),Vector4i(4, 1, 0, 5),Vector4i(5, 1, 0, 2)
				},
				{
					Vector4i(3, 0, 5, 4),Vector4i(4, 2, 0, 5),Vector4i(4, 1, 0, 2)
				}
			},
			{

				{
					Vector4i(1, 4, 5, 7),Vector4i(3, 7, 6, 5),Vector4i(0, 7, 4, 2),Vector4i(2, 7, 4, 1),Vector4i(0, 2, 6, 7),Vector4i(0, 7, 6, 3),Vector4i(0, 7, 5, 4),Vector4i(0, 7, 3, 5),Vector4i(1, 5, 6, 7),Vector4i(1, 6, 2, 7)
				},
				{
					Vector4i(1, 4, 5, 7),Vector4i(3, 7, 6, 5),Vector4i(0, 7, 4, 2),Vector4i(2, 7, 4, 1),Vector4i(0, 2, 6, 7),Vector4i(0, 7, 6, 3),Vector4i(0, 7, 3, 4),Vector4i(4, 7, 3, 5),Vector4i(1, 5, 6, 7),Vector4i(1, 6, 2, 7)
				},
				{
					Vector4i(1, 4, 5, 7),Vector4i(3, 7, 6, 5),Vector4i(0, 7, 4, 2),Vector4i(2, 7, 4, 1),Vector4i(0, 2, 6, 7),Vector4i(0, 7, 6, 3),Vector4i(0, 7, 5, 4),Vector4i(0, 7, 3, 5),Vector4i(1, 5, 2, 7),Vector4i(5, 6, 2, 7)
				},
				{
					Vector4i(1, 4, 5, 7),Vector4i(3, 7, 6, 5),Vector4i(0, 7, 4, 2),Vector4i(2, 7, 4, 1),Vector4i(0, 2, 6, 7),Vector4i(0, 7, 6, 3),Vector4i(0, 7, 3, 4),Vector4i(4, 7, 3, 5),Vector4i(1, 5, 2, 7),Vector4i(5, 6, 2, 7)
				}
			},
			{},
			{},
			{

				{
					Vector4i(2, 7, 6, 4),Vector4i(3, 7, 6, 5),Vector4i(0, 1, 4, 7),Vector4i(1, 2, 4, 7),Vector4i(0, 7, 5, 1),Vector4i(0, 3, 5, 7),Vector4i(0, 4, 6, 7),Vector4i(0, 6, 3, 7),Vector4i(1, 5, 6, 7),Vector4i(1, 6, 2, 7)
				},
				{
					Vector4i(2, 7, 6, 4),Vector4i(3, 7, 6, 5),Vector4i(0, 1, 4, 7),Vector4i(1, 2, 4, 7),Vector4i(0, 7, 5, 1),Vector4i(0, 3, 5, 7),Vector4i(0, 4, 3, 7),Vector4i(4, 6, 3, 7),Vector4i(1, 5, 6, 7),Vector4i(1, 6, 2, 7)
				},
				{
					Vector4i(2, 7, 6, 4),Vector4i(3, 7, 6, 5),Vector4i(0, 1, 4, 7),Vector4i(1, 2, 4, 7),Vector4i(0, 7, 5, 1),Vector4i(0, 3, 5, 7),Vector4i(0, 4, 6, 7),Vector4i(0, 6, 3, 7),Vector4i(1, 5, 2, 7),Vector4i(5, 6, 2, 7)
				},
				{
					Vector4i(2, 7, 6, 4),Vector4i(3, 7, 6, 5),Vector4i(0, 1, 4, 7),Vector4i(1, 2, 4, 7),Vector4i(0, 7, 5, 1),Vector4i(0, 3, 5, 7),Vector4i(0, 4, 3, 7),Vector4i(4, 6, 3, 7),Vector4i(1, 5, 2, 7),Vector4i(5, 6, 2, 7)
				}
			},
			{

				{
					Vector4i(2, 4, 7, 5),Vector4i(7, 1, 6, 4),Vector4i(2, 1, 7, 4),Vector4i(3, 4, 7, 6),Vector4i(7, 0, 5, 4),Vector4i(3, 0, 7, 4)
				},
				{
					Vector4i(2, 4, 7, 5),Vector4i(7, 1, 6, 4),Vector4i(2, 1, 7, 4),Vector4i(3, 4, 7, 6),Vector4i(7, 3, 5, 4),Vector4i(3, 0, 5, 4)
				},
				{
					Vector4i(2, 4, 7, 5),Vector4i(7, 1, 6, 4),Vector4i(2, 1, 7, 4),Vector4i(0, 4, 7, 6),Vector4i(7, 0, 5, 4),Vector4i(3, 0, 7, 6)
				},
				{
					Vector4i(2, 4, 7, 5),Vector4i(7, 1, 6, 4),Vector4i(2, 1, 7, 4),Vector4i(0, 8, 5, 3),Vector4i(5, 7, 3, 8),Vector4i(6, 7, 4, 8),Vector4i(5, 8, 4, 7),Vector4i(0, 6, 4, 8),Vector4i(3, 6, 0, 8),Vector4i(0, 4, 5, 8),Vector4i(3, 8, 7, 6)
				},
				{
					Vector4i(2, 4, 7, 5),Vector4i(7, 2, 6, 4),Vector4i(2, 1, 6, 4),Vector4i(3, 4, 7, 6),Vector4i(7, 0, 5, 4),Vector4i(3, 0, 7, 4)
				},
				{
					Vector4i(2, 4, 7, 5),Vector4i(7, 2, 6, 4),Vector4i(2, 1, 6, 4),Vector4i(3, 4, 7, 6),Vector4i(7, 3, 5, 4),Vector4i(3, 0, 5, 4)
				},
				{
					Vector4i(2, 4, 7, 5),Vector4i(7, 2, 6, 4),Vector4i(2, 1, 6, 4),Vector4i(0, 4, 7, 6),Vector4i(7, 0, 5, 4),Vector4i(3, 0, 7, 6)
				},
				{
					Vector4i(2, 4, 7, 5),Vector4i(7, 2, 6, 4),Vector4i(2, 1, 6, 4),Vector4i(0, 8, 5, 3),Vector4i(5, 7, 3, 8),Vector4i(6, 7, 4, 8),Vector4i(5, 8, 4, 7),Vector4i(0, 6, 4, 8),Vector4i(3, 6, 0, 8),Vector4i(0, 4, 5, 8),Vector4i(3, 8, 7, 6)
				},
				{
					Vector4i(1, 4, 7, 5),Vector4i(7, 1, 6, 4),Vector4i(2, 1, 7, 5),Vector4i(3, 4, 7, 6),Vector4i(7, 0, 5, 4),Vector4i(3, 0, 7, 4)
				},
				{
					Vector4i(1, 4, 7, 5),Vector4i(7, 1, 6, 4),Vector4i(2, 1, 7, 5),Vector4i(3, 4, 7, 6),Vector4i(7, 3, 5, 4),Vector4i(3, 0, 5, 4)
				},
				{
					Vector4i(1, 4, 7, 5),Vector4i(7, 1, 6, 4),Vector4i(2, 1, 7, 5),Vector4i(0, 4, 7, 6),Vector4i(7, 0, 5, 4),Vector4i(3, 0, 7, 6)
				},
				{
					Vector4i(1, 4, 7, 5),Vector4i(7, 1, 6, 4),Vector4i(2, 1, 7, 5),Vector4i(0, 8, 5, 3),Vector4i(5, 7, 3, 8),Vector4i(6, 7, 4, 8),Vector4i(5, 8, 4, 7),Vector4i(0, 6, 4, 8),Vector4i(3, 6, 0, 8),Vector4i(0, 4, 5, 8),Vector4i(3, 8, 7, 6)
				},
				{
					Vector4i(2, 4, 6, 5),Vector4i(7, 2, 6, 5),Vector4i(2, 1, 6, 4),Vector4i(3, 4, 5, 6),Vector4i(7, 3, 5, 6),Vector4i(3, 0, 5, 4)
				},
				{
					Vector4i(2, 4, 6, 5),Vector4i(7, 2, 6, 5),Vector4i(2, 1, 6, 4),Vector4i(0, 4, 5, 6),Vector4i(0, 5, 7, 6),Vector4i(3, 0, 7, 6)
				},
				{
					Vector4i(2, 4, 6, 5),Vector4i(7, 2, 6, 5),Vector4i(2, 1, 6, 4),Vector4i(0, 4, 5, 6),Vector4i(3, 5, 7, 6),Vector4i(3, 0, 5, 6)
				},
				{
					Vector4i(2, 4, 6, 5),Vector4i(7, 2, 6, 5),Vector4i(2, 1, 6, 4),Vector4i(0, 8, 7, 3),Vector4i(0, 8, 5, 7),Vector4i(6, 7, 5, 8),Vector4i(5, 8, 4, 6),Vector4i(3, 4, 0, 8),Vector4i(3, 6, 4, 8),Vector4i(0, 4, 5, 8),Vector4i(3, 8, 7, 6)
				},
				{
					Vector4i(1, 4, 6, 5),Vector4i(1, 6, 7, 5),Vector4i(2, 1, 7, 5),Vector4i(3, 4, 5, 6),Vector4i(7, 3, 5, 6),Vector4i(3, 0, 5, 4)
				},
				{
					Vector4i(1, 4, 6, 5),Vector4i(1, 6, 7, 5),Vector4i(2, 1, 7, 5),Vector4i(0, 4, 5, 6),Vector4i(0, 5, 7, 6),Vector4i(3, 0, 7, 6)
				},
				{
					Vector4i(1, 4, 6, 5),Vector4i(1, 6, 7, 5),Vector4i(2, 1, 7, 5),Vector4i(0, 4, 5, 6),Vector4i(3, 5, 7, 6),Vector4i(3, 0, 5, 6)
				},
				{
					Vector4i(1, 4, 6, 5),Vector4i(1, 6, 7, 5),Vector4i(2, 1, 7, 5),Vector4i(0, 8, 7, 3),Vector4i(0, 8, 5, 7),Vector4i(6, 7, 5, 8),Vector4i(5, 8, 4, 6),Vector4i(3, 4, 0, 8),Vector4i(3, 6, 4, 8),Vector4i(0, 4, 5, 8),Vector4i(3, 8, 7, 6)
				},
				{
					Vector4i(1, 4, 6, 5),Vector4i(2, 6, 7, 5),Vector4i(2, 1, 6, 5),Vector4i(3, 4, 5, 6),Vector4i(7, 3, 5, 6),Vector4i(3, 0, 5, 4)
				},
				{
					Vector4i(1, 4, 6, 5),Vector4i(2, 6, 7, 5),Vector4i(2, 1, 6, 5),Vector4i(0, 4, 5, 6),Vector4i(0, 5, 7, 6),Vector4i(3, 0, 7, 6)
				},
				{
					Vector4i(1, 4, 6, 5),Vector4i(2, 6, 7, 5),Vector4i(2, 1, 6, 5),Vector4i(0, 4, 5, 6),Vector4i(3, 5, 7, 6),Vector4i(3, 0, 5, 6)
				},
				{
					Vector4i(1, 4, 6, 5),Vector4i(2, 6, 7, 5),Vector4i(2, 1, 6, 5),Vector4i(0, 8, 7, 3),Vector4i(0, 8, 5, 7),Vector4i(6, 7, 5, 8),Vector4i(5, 8, 4, 6),Vector4i(3, 4, 0, 8),Vector4i(3, 6, 4, 8),Vector4i(0, 4, 5, 8),Vector4i(3, 8, 7, 6)
				},
				{
					Vector4i(1, 8, 6, 2),Vector4i(6, 7, 2, 8),Vector4i(5, 7, 4, 8),Vector4i(6, 8, 4, 7),Vector4i(1, 5, 4, 8),Vector4i(2, 5, 1, 8),Vector4i(1, 4, 6, 8),Vector4i(2, 8, 7, 5),Vector4i(3, 4, 7, 6),Vector4i(7, 0, 5, 4),Vector4i(3, 0, 7, 4)
				},
				{
					Vector4i(1, 8, 6, 2),Vector4i(6, 7, 2, 8),Vector4i(5, 7, 4, 8),Vector4i(6, 8, 4, 7),Vector4i(1, 5, 4, 8),Vector4i(2, 5, 1, 8),Vector4i(1, 4, 6, 8),Vector4i(2, 8, 7, 5),Vector4i(3, 4, 7, 6),Vector4i(7, 3, 5, 4),Vector4i(3, 0, 5, 4)
				},
				{
					Vector4i(1, 8, 6, 2),Vector4i(6, 7, 2, 8),Vector4i(5, 7, 4, 8),Vector4i(6, 8, 4, 7),Vector4i(1, 5, 4, 8),Vector4i(2, 5, 1, 8),Vector4i(1, 4, 6, 8),Vector4i(2, 8, 7, 5),Vector4i(0, 4, 7, 6),Vector4i(7, 0, 5, 4),Vector4i(3, 0, 7, 6)
				},
				{
					Vector4i(1, 8, 6, 2),Vector4i(6, 7, 2, 8),Vector4i(5, 7, 4, 8),Vector4i(6, 8, 4, 7),Vector4i(1, 5, 4, 8),Vector4i(2, 5, 1, 8),Vector4i(1, 4, 6, 8),Vector4i(2, 8, 7, 5),Vector4i(0, 9, 5, 3),Vector4i(5, 7, 3, 9),Vector4i(6, 7, 4, 9),Vector4i(5, 9, 4, 7),Vector4i(0, 6, 4, 9),Vector4i(3, 6, 0, 9),Vector4i(0, 4, 5, 9),Vector4i(3, 9, 7, 6)
				},
				{
					Vector4i(1, 8, 7, 2),Vector4i(1, 8, 6, 7),Vector4i(5, 7, 6, 8),Vector4i(6, 8, 4, 5),Vector4i(2, 4, 1, 8),Vector4i(2, 5, 4, 8),Vector4i(1, 4, 6, 8),Vector4i(2, 8, 7, 5),Vector4i(3, 4, 5, 6),Vector4i(7, 3, 5, 6),Vector4i(3, 0, 5, 4)
				},
				{
					Vector4i(1, 8, 7, 2),Vector4i(1, 8, 6, 7),Vector4i(5, 7, 6, 8),Vector4i(6, 8, 4, 5),Vector4i(2, 4, 1, 8),Vector4i(2, 5, 4, 8),Vector4i(1, 4, 6, 8),Vector4i(2, 8, 7, 5),Vector4i(0, 4, 5, 6),Vector4i(0, 5, 7, 6),Vector4i(3, 0, 7, 6)
				},
				{
					Vector4i(1, 8, 7, 2),Vector4i(1, 8, 6, 7),Vector4i(5, 7, 6, 8),Vector4i(6, 8, 4, 5),Vector4i(2, 4, 1, 8),Vector4i(2, 5, 4, 8),Vector4i(1, 4, 6, 8),Vector4i(2, 8, 7, 5),Vector4i(0, 4, 5, 6),Vector4i(3, 5, 7, 6),Vector4i(3, 0, 5, 6)
				},
				{
					Vector4i(1, 8, 7, 2),Vector4i(1, 8, 6, 7),Vector4i(5, 7, 6, 8),Vector4i(6, 8, 4, 5),Vector4i(2, 4, 1, 8),Vector4i(2, 5, 4, 8),Vector4i(1, 4, 6, 8),Vector4i(2, 8, 7, 5),Vector4i(0, 9, 7, 3),Vector4i(0, 9, 5, 7),Vector4i(6, 7, 5, 9),Vector4i(5, 9, 4, 6),Vector4i(3, 4, 0, 9),Vector4i(3, 6, 4, 9),Vector4i(0, 4, 5, 9),Vector4i(3, 9, 7, 6)
				}
			},
			{},
			{},
			{

				{
					Vector4i(4, 5, 6, 3),Vector4i(4, 1, 6, 5),Vector4i(6, 0, 2, 1),Vector4i(4, 0, 6, 1)
				},
				{
					Vector4i(4, 5, 6, 3),Vector4i(4, 1, 6, 5),Vector4i(6, 4, 2, 1),Vector4i(4, 0, 2, 1)
				},
				{
					Vector4i(4, 5, 6, 3),Vector4i(0, 1, 6, 5),Vector4i(6, 0, 2, 1),Vector4i(4, 0, 6, 5)
				},
				{
					Vector4i(4, 5, 6, 3),Vector4i(4, 1, 2, 5),Vector4i(6, 4, 2, 5),Vector4i(4, 0, 2, 1)
				},
				{
					Vector4i(4, 5, 6, 3),Vector4i(0, 1, 2, 5),Vector4i(0, 2, 6, 5),Vector4i(4, 0, 6, 5)
				},
				{
					Vector4i(4, 5, 6, 3),Vector4i(0, 1, 2, 5),Vector4i(4, 2, 6, 5),Vector4i(4, 0, 2, 5)
				},
				{
					Vector4i(4, 5, 6, 3),Vector4i(0, 7, 2, 4),Vector4i(2, 6, 4, 7),Vector4i(5, 6, 1, 7),Vector4i(2, 7, 1, 6),Vector4i(0, 5, 1, 7),Vector4i(4, 5, 0, 7),Vector4i(0, 1, 2, 7),Vector4i(4, 7, 6, 5)
				},
				{
					Vector4i(4, 5, 6, 3),Vector4i(0, 7, 6, 4),Vector4i(0, 7, 2, 6),Vector4i(5, 6, 2, 7),Vector4i(2, 7, 1, 5),Vector4i(4, 1, 0, 7),Vector4i(4, 5, 1, 7),Vector4i(0, 1, 2, 7),Vector4i(4, 7, 6, 5)
				}
			},
			{},
			{},
			{},
			{},
			{},
			{},
			{}
		}};

		assert(!table[idx].empty());
		return table[idx];
	}


	const std::vector<std::vector<Vector2i>>& CutTable::get_diag_confs(const int idx) {
		static const std::array<std::vector<std::vector<Vector2i>>, 64> table= {{

			{

				{
				
				}
			},
			{

				{
				
				}
			},
			{

				{
				
				}
			},
			{

				{
					Vector2i(0, 5)
				},
				{
					Vector2i(2, 4)
				}
			},
			{

				{
				
				}
			},
			{

				{
					Vector2i(1, 5)
				},
				{
					Vector2i(2, 4)
				}
			},
			{

				{
					Vector2i(1, 5)
				},
				{
					Vector2i(0, 4)
				}
			},
			{},
			{

				{
				
				}
			},
			{

				{
					Vector2i(1, 5)
				},
				{
					Vector2i(3, 4)
				}
			},
			{

				{
				
				}
			},
			{

				{
					Vector2i(0, 5),Vector2i(1, 6)
				},
				{
					Vector2i(1, 6),Vector2i(2, 4)
				},
				{
					Vector2i(0, 5),Vector2i(3, 4)
				},
				{
					Vector2i(2, 4),Vector2i(3, 4)
				}
			},
			{

				{
					Vector2i(2, 5)
				},
				{
					Vector2i(3, 4)
				}
			},
			{

				{
					Vector2i(1, 6),Vector2i(2, 4),Vector2i(2, 6)
				},
				{
					Vector2i(2, 4),Vector2i(2, 6),Vector2i(3, 4)
				},
				{
					Vector2i(1, 5),Vector2i(1, 6),Vector2i(2, 6)
				},
				{
					Vector2i(2, 4),Vector2i(3, 4),Vector2i(3, 5)
				},
				{
					Vector2i(1, 5),Vector2i(1, 6),Vector2i(3, 5)
				},
				{
					Vector2i(1, 5),Vector2i(3, 4),Vector2i(3, 5)
				},
				{
					Vector2i(1, 5),Vector2i(2, 6),Vector2i(3, 4)
				},
				{
					Vector2i(1, 6),Vector2i(2, 4),Vector2i(3, 5)
				}
			},
			{

				{
					Vector2i(1, 5),Vector2i(2, 6)
				},
				{
					Vector2i(0, 4),Vector2i(2, 6)
				},
				{
					Vector2i(1, 5),Vector2i(3, 5)
				},
				{
					Vector2i(0, 4),Vector2i(3, 5)
				}
			},
			{},
			{

				{
				
				}
			},
			{

				{
					Vector2i(0, 5)
				},
				{
					Vector2i(3, 4)
				}
			},
			{

				{
					Vector2i(2, 5)
				},
				{
					Vector2i(3, 4)
				}
			},
			{

				{
					Vector2i(0, 6),Vector2i(2, 4),Vector2i(2, 6)
				},
				{
					Vector2i(2, 4),Vector2i(2, 6),Vector2i(3, 4)
				},
				{
					Vector2i(0, 5),Vector2i(0, 6),Vector2i(2, 6)
				},
				{
					Vector2i(2, 4),Vector2i(3, 4),Vector2i(3, 5)
				},
				{
					Vector2i(0, 5),Vector2i(0, 6),Vector2i(3, 5)
				},
				{
					Vector2i(0, 5),Vector2i(3, 4),Vector2i(3, 5)
				},
				{
					Vector2i(0, 5),Vector2i(2, 6),Vector2i(3, 4)
				},
				{
					Vector2i(0, 6),Vector2i(2, 4),Vector2i(3, 5)
				}
			},
			{

				{
				
				}
			},
			{

				{
					Vector2i(0, 6),Vector2i(1, 5)
				},
				{
					Vector2i(0, 6),Vector2i(2, 4)
				},
				{
					Vector2i(1, 5),Vector2i(3, 4)
				},
				{
					Vector2i(2, 4),Vector2i(3, 4)
				}
			},
			{

				{
					Vector2i(1, 5),Vector2i(2, 6)
				},
				{
					Vector2i(0, 4),Vector2i(2, 6)
				},
				{
					Vector2i(1, 5),Vector2i(3, 4)
				},
				{
					Vector2i(0, 4),Vector2i(3, 4)
				}
			},
			{},
			{

				{
					Vector2i(0, 5)
				},
				{
					Vector2i(1, 4)
				}
			},
			{},
			{

				{
					Vector2i(0, 6),Vector2i(2, 6)
				},
				{
					Vector2i(0, 6),Vector2i(3, 4)
				},
				{
					Vector2i(1, 5),Vector2i(2, 6)
				},
				{
					Vector2i(1, 5),Vector2i(3, 4)
				}
			},
			{},
			{

				{
					Vector2i(0, 6),Vector2i(2, 5)
				},
				{
					Vector2i(0, 6),Vector2i(3, 4)
				},
				{
					Vector2i(1, 5),Vector2i(2, 5)
				},
				{
					Vector2i(1, 5),Vector2i(3, 4)
				}
			},
			{},
			{

				{
					Vector2i(0, 7),Vector2i(1, 5),Vector2i(2, 7),Vector2i(3, 5)
				},
				{
					Vector2i(0, 7),Vector2i(1, 5),Vector2i(3, 4),Vector2i(3, 5)
				},
				{
					Vector2i(0, 7),Vector2i(1, 5),Vector2i(2, 6),Vector2i(2, 7)
				},
				{
					Vector2i(0, 7),Vector2i(1, 5),Vector2i(2, 6),Vector2i(3, 4)
				},
				{
					Vector2i(1, 5),Vector2i(1, 6),Vector2i(2, 7),Vector2i(3, 5)
				},
				{
					Vector2i(1, 5),Vector2i(1, 6),Vector2i(3, 4),Vector2i(3, 5)
				},
				{
					Vector2i(1, 5),Vector2i(1, 6),Vector2i(2, 6),Vector2i(2, 7)
				},
				{
					Vector2i(1, 5),Vector2i(1, 6),Vector2i(2, 6),Vector2i(3, 4)
				},
				{
					Vector2i(0, 4),Vector2i(0, 7),Vector2i(2, 7),Vector2i(3, 5)
				},
				{
					Vector2i(0, 4),Vector2i(0, 7),Vector2i(3, 4),Vector2i(3, 5)
				},
				{
					Vector2i(0, 4),Vector2i(0, 7),Vector2i(2, 6),Vector2i(2, 7)
				},
				{
					Vector2i(0, 4),Vector2i(0, 7),Vector2i(2, 6),Vector2i(3, 4)
				},
				{
					Vector2i(1, 5),Vector2i(1, 6),Vector2i(3, 4),Vector2i(3, 5)
				},
				{
					Vector2i(1, 5),Vector2i(1, 6),Vector2i(2, 6),Vector2i(2, 7)
				},
				{
					Vector2i(1, 5),Vector2i(1, 6),Vector2i(2, 6),Vector2i(3, 4)
				},
				{
					Vector2i(1, 5),Vector2i(1, 6),Vector2i(2, 7),Vector2i(3, 5)
				},
				{
					Vector2i(0, 4),Vector2i(0, 7),Vector2i(3, 4),Vector2i(3, 5)
				},
				{
					Vector2i(0, 4),Vector2i(0, 7),Vector2i(2, 6),Vector2i(2, 7)
				},
				{
					Vector2i(0, 4),Vector2i(0, 7),Vector2i(2, 6),Vector2i(3, 4)
				},
				{
					Vector2i(0, 4),Vector2i(0, 7),Vector2i(2, 7),Vector2i(3, 5)
				},
				{
					Vector2i(0, 4),Vector2i(1, 6),Vector2i(3, 4),Vector2i(3, 5)
				},
				{
					Vector2i(0, 4),Vector2i(1, 6),Vector2i(2, 6),Vector2i(2, 7)
				},
				{
					Vector2i(0, 4),Vector2i(1, 6),Vector2i(2, 6),Vector2i(3, 4)
				},
				{
					Vector2i(0, 4),Vector2i(1, 6),Vector2i(2, 7),Vector2i(3, 5)
				},
				{
					Vector2i(0, 4),Vector2i(1, 6),Vector2i(2, 7),Vector2i(3, 5)
				},
				{
					Vector2i(0, 4),Vector2i(1, 6),Vector2i(3, 4),Vector2i(3, 5)
				},
				{
					Vector2i(0, 4),Vector2i(1, 6),Vector2i(2, 6),Vector2i(2, 7)
				},
				{
					Vector2i(0, 4),Vector2i(1, 6),Vector2i(2, 6),Vector2i(3, 4)
				},
				{
					Vector2i(0, 7),Vector2i(1, 5),Vector2i(3, 4),Vector2i(3, 5)
				},
				{
					Vector2i(0, 7),Vector2i(1, 5),Vector2i(2, 6),Vector2i(2, 7)
				},
				{
					Vector2i(0, 7),Vector2i(1, 5),Vector2i(2, 6),Vector2i(3, 4)
				},
				{
					Vector2i(0, 7),Vector2i(1, 5),Vector2i(2, 7),Vector2i(3, 5)
				}
			},
			{},
			{

				{
				
				}
			},
			{

				{
				
				}
			},
			{

				{
					Vector2i(1, 5)
				},
				{
					Vector2i(3, 4)
				}
			},
			{

				{
					Vector2i(0, 5),Vector2i(1, 6)
				},
				{
					Vector2i(1, 6),Vector2i(2, 4)
				},
				{
					Vector2i(0, 5),Vector2i(3, 5)
				},
				{
					Vector2i(2, 4),Vector2i(3, 5)
				}
			},
			{

				{
					Vector2i(0, 5)
				},
				{
					Vector2i(3, 4)
				}
			},
			{

				{
					Vector2i(0, 6),Vector2i(1, 5)
				},
				{
					Vector2i(0, 6),Vector2i(2, 4)
				},
				{
					Vector2i(1, 5),Vector2i(3, 5)
				},
				{
					Vector2i(2, 4),Vector2i(3, 5)
				}
			},
			{

				{
					Vector2i(0, 4),Vector2i(0, 6),Vector2i(1, 6)
				},
				{
					Vector2i(0, 4),Vector2i(0, 6),Vector2i(3, 4)
				},
				{
					Vector2i(0, 6),Vector2i(1, 5),Vector2i(1, 6)
				},
				{
					Vector2i(0, 4),Vector2i(3, 4),Vector2i(3, 5)
				},
				{
					Vector2i(1, 5),Vector2i(1, 6),Vector2i(3, 5)
				},
				{
					Vector2i(1, 5),Vector2i(3, 4),Vector2i(3, 5)
				},
				{
					Vector2i(0, 6),Vector2i(1, 5),Vector2i(3, 4)
				},
				{
					Vector2i(0, 4),Vector2i(1, 6),Vector2i(3, 5)
				}
			},
			{},
			{

				{
					Vector2i(0, 5)
				},
				{
					Vector2i(2, 4)
				}
			},
			{

				{
					Vector2i(0, 6),Vector2i(1, 5)
				},
				{
					Vector2i(0, 6),Vector2i(3, 4)
				},
				{
					Vector2i(1, 5),Vector2i(2, 5)
				},
				{
					Vector2i(2, 5),Vector2i(3, 4)
				}
			},
			{

				{
					Vector2i(0, 6),Vector2i(1, 6)
				},
				{
					Vector2i(0, 6),Vector2i(3, 4)
				},
				{
					Vector2i(1, 6),Vector2i(2, 5)
				},
				{
					Vector2i(2, 5),Vector2i(3, 4)
				}
			},
			{

				{
					Vector2i(0, 7),Vector2i(1, 7),Vector2i(2, 4),Vector2i(3, 4)
				},
				{
					Vector2i(0, 7),Vector2i(2, 4),Vector2i(3, 4),Vector2i(3, 5)
				},
				{
					Vector2i(0, 7),Vector2i(1, 6),Vector2i(1, 7),Vector2i(2, 4)
				},
				{
					Vector2i(0, 7),Vector2i(1, 6),Vector2i(2, 4),Vector2i(3, 5)
				},
				{
					Vector2i(1, 7),Vector2i(2, 4),Vector2i(2, 6),Vector2i(3, 4)
				},
				{
					Vector2i(2, 4),Vector2i(2, 6),Vector2i(3, 4),Vector2i(3, 5)
				},
				{
					Vector2i(1, 6),Vector2i(1, 7),Vector2i(2, 4),Vector2i(2, 6)
				},
				{
					Vector2i(1, 6),Vector2i(2, 4),Vector2i(2, 6),Vector2i(3, 5)
				},
				{
					Vector2i(0, 5),Vector2i(0, 7),Vector2i(1, 7),Vector2i(3, 4)
				},
				{
					Vector2i(0, 5),Vector2i(0, 7),Vector2i(3, 4),Vector2i(3, 5)
				},
				{
					Vector2i(0, 5),Vector2i(0, 7),Vector2i(1, 6),Vector2i(1, 7)
				},
				{
					Vector2i(0, 5),Vector2i(0, 7),Vector2i(1, 6),Vector2i(3, 5)
				},
				{
					Vector2i(2, 4),Vector2i(2, 6),Vector2i(3, 4),Vector2i(3, 5)
				},
				{
					Vector2i(1, 6),Vector2i(1, 7),Vector2i(2, 4),Vector2i(2, 6)
				},
				{
					Vector2i(1, 6),Vector2i(2, 4),Vector2i(2, 6),Vector2i(3, 5)
				},
				{
					Vector2i(1, 7),Vector2i(2, 4),Vector2i(2, 6),Vector2i(3, 4)
				},
				{
					Vector2i(0, 5),Vector2i(0, 7),Vector2i(3, 4),Vector2i(3, 5)
				},
				{
					Vector2i(0, 5),Vector2i(0, 7),Vector2i(1, 6),Vector2i(1, 7)
				},
				{
					Vector2i(0, 5),Vector2i(0, 7),Vector2i(1, 6),Vector2i(3, 5)
				},
				{
					Vector2i(0, 5),Vector2i(0, 7),Vector2i(1, 7),Vector2i(3, 4)
				},
				{
					Vector2i(0, 5),Vector2i(2, 6),Vector2i(3, 4),Vector2i(3, 5)
				},
				{
					Vector2i(0, 5),Vector2i(1, 6),Vector2i(1, 7),Vector2i(2, 6)
				},
				{
					Vector2i(0, 5),Vector2i(1, 6),Vector2i(2, 6),Vector2i(3, 5)
				},
				{
					Vector2i(0, 5),Vector2i(1, 7),Vector2i(2, 6),Vector2i(3, 4)
				},
				{
					Vector2i(0, 5),Vector2i(1, 7),Vector2i(2, 6),Vector2i(3, 4)
				},
				{
					Vector2i(0, 5),Vector2i(2, 6),Vector2i(3, 4),Vector2i(3, 5)
				},
				{
					Vector2i(0, 5),Vector2i(1, 6),Vector2i(1, 7),Vector2i(2, 6)
				},
				{
					Vector2i(0, 5),Vector2i(1, 6),Vector2i(2, 6),Vector2i(3, 5)
				},
				{
					Vector2i(0, 7),Vector2i(2, 4),Vector2i(3, 4),Vector2i(3, 5)
				},
				{
					Vector2i(0, 7),Vector2i(1, 6),Vector2i(1, 7),Vector2i(2, 4)
				},
				{
					Vector2i(0, 7),Vector2i(1, 6),Vector2i(2, 4),Vector2i(3, 5)
				},
				{
					Vector2i(0, 7),Vector2i(1, 7),Vector2i(2, 4),Vector2i(3, 4)
				}
			},
			{},
			{},
			{},
			{},
			{

				{
					Vector2i(1, 5)
				},
				{
					Vector2i(2, 4)
				}
			},
			{

				{
					Vector2i(0, 5),Vector2i(1, 6)
				},
				{
					Vector2i(1, 6),Vector2i(3, 4)
				},
				{
					Vector2i(0, 5),Vector2i(2, 5)
				},
				{
					Vector2i(2, 5),Vector2i(3, 4)
				}
			},
			{},
			{},
			{

				{
					Vector2i(0, 6),Vector2i(1, 6)
				},
				{
					Vector2i(1, 6),Vector2i(3, 4)
				},
				{
					Vector2i(0, 6),Vector2i(2, 5)
				},
				{
					Vector2i(2, 5),Vector2i(3, 4)
				}
			},
			{

				{
					Vector2i(0, 7),Vector2i(1, 7),Vector2i(2, 4),Vector2i(3, 4)
				},
				{
					Vector2i(1, 7),Vector2i(2, 4),Vector2i(3, 4),Vector2i(3, 5)
				},
				{
					Vector2i(0, 6),Vector2i(0, 7),Vector2i(1, 7),Vector2i(2, 4)
				},
				{
					Vector2i(0, 6),Vector2i(1, 7),Vector2i(2, 4),Vector2i(3, 5)
				},
				{
					Vector2i(0, 7),Vector2i(2, 4),Vector2i(2, 6),Vector2i(3, 4)
				},
				{
					Vector2i(2, 4),Vector2i(2, 6),Vector2i(3, 4),Vector2i(3, 5)
				},
				{
					Vector2i(0, 6),Vector2i(0, 7),Vector2i(2, 4),Vector2i(2, 6)
				},
				{
					Vector2i(0, 6),Vector2i(2, 4),Vector2i(2, 6),Vector2i(3, 5)
				},
				{
					Vector2i(0, 7),Vector2i(1, 5),Vector2i(1, 7),Vector2i(3, 4)
				},
				{
					Vector2i(1, 5),Vector2i(1, 7),Vector2i(3, 4),Vector2i(3, 5)
				},
				{
					Vector2i(0, 6),Vector2i(0, 7),Vector2i(1, 5),Vector2i(1, 7)
				},
				{
					Vector2i(0, 6),Vector2i(1, 5),Vector2i(1, 7),Vector2i(3, 5)
				},
				{
					Vector2i(2, 4),Vector2i(2, 6),Vector2i(3, 4),Vector2i(3, 5)
				},
				{
					Vector2i(0, 6),Vector2i(0, 7),Vector2i(2, 4),Vector2i(2, 6)
				},
				{
					Vector2i(0, 6),Vector2i(2, 4),Vector2i(2, 6),Vector2i(3, 5)
				},
				{
					Vector2i(0, 7),Vector2i(2, 4),Vector2i(2, 6),Vector2i(3, 4)
				},
				{
					Vector2i(1, 5),Vector2i(1, 7),Vector2i(3, 4),Vector2i(3, 5)
				},
				{
					Vector2i(0, 6),Vector2i(0, 7),Vector2i(1, 5),Vector2i(1, 7)
				},
				{
					Vector2i(0, 6),Vector2i(1, 5),Vector2i(1, 7),Vector2i(3, 5)
				},
				{
					Vector2i(0, 7),Vector2i(1, 5),Vector2i(1, 7),Vector2i(3, 4)
				},
				{
					Vector2i(1, 5),Vector2i(2, 6),Vector2i(3, 4),Vector2i(3, 5)
				},
				{
					Vector2i(0, 6),Vector2i(0, 7),Vector2i(1, 5),Vector2i(2, 6)
				},
				{
					Vector2i(0, 6),Vector2i(1, 5),Vector2i(2, 6),Vector2i(3, 5)
				},
				{
					Vector2i(0, 7),Vector2i(1, 5),Vector2i(2, 6),Vector2i(3, 4)
				},
				{
					Vector2i(0, 7),Vector2i(1, 5),Vector2i(2, 6),Vector2i(3, 4)
				},
				{
					Vector2i(1, 5),Vector2i(2, 6),Vector2i(3, 4),Vector2i(3, 5)
				},
				{
					Vector2i(0, 6),Vector2i(0, 7),Vector2i(1, 5),Vector2i(2, 6)
				},
				{
					Vector2i(0, 6),Vector2i(1, 5),Vector2i(2, 6),Vector2i(3, 5)
				},
				{
					Vector2i(1, 7),Vector2i(2, 4),Vector2i(3, 4),Vector2i(3, 5)
				},
				{
					Vector2i(0, 6),Vector2i(0, 7),Vector2i(1, 7),Vector2i(2, 4)
				},
				{
					Vector2i(0, 6),Vector2i(1, 7),Vector2i(2, 4),Vector2i(3, 5)
				},
				{
					Vector2i(0, 7),Vector2i(1, 7),Vector2i(2, 4),Vector2i(3, 4)
				}
			},
			{},
			{},
			{

				{
					Vector2i(0, 6),Vector2i(1, 4),Vector2i(1, 6)
				},
				{
					Vector2i(1, 4),Vector2i(1, 6),Vector2i(2, 4)
				},
				{
					Vector2i(0, 5),Vector2i(0, 6),Vector2i(1, 6)
				},
				{
					Vector2i(1, 4),Vector2i(2, 4),Vector2i(2, 5)
				},
				{
					Vector2i(0, 5),Vector2i(0, 6),Vector2i(2, 5)
				},
				{
					Vector2i(0, 5),Vector2i(2, 4),Vector2i(2, 5)
				},
				{
					Vector2i(0, 5),Vector2i(1, 6),Vector2i(2, 4)
				},
				{
					Vector2i(0, 6),Vector2i(1, 4),Vector2i(2, 5)
				}
			},
			{},
			{},
			{},
			{},
			{},
			{},
			{}
		}};

		//assert(!table[idx].empty());
		return table[idx];
	}


	const std::vector<std::vector<std::array<bool, 4>>>& CutTable::get_surface_conf(const int idx) {
		static const std::array<std::vector<std::vector<std::array<bool, 4>>>, 64> table= {{

			{

				{
					{{false, false, false, false}}
				}
			},
			{

				{
					{{false, false, false, true}},{{false, true, false, false}}
				}
			},
			{

				{
					{{false, false, false, true}},{{false, true, false, false}}
				}
			},
			{

				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				}
			},
			{

				{
					{{false, false, false, true}},{{false, true, false, false}}
				}
			},
			{

				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				}
			},
			{

				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				}
			},
			{},
			{

				{
					{{false, false, false, true}},{{false, true, false, false}}
				}
			},
			{

				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				}
			},
			{

				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				}
			},
			{

				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				}
			},
			{

				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				}
			},
			{

				{
					{{false, true, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, true, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}}
				},
				{
					{{false, true, false, false}},{{false, false, false, false}},{{false, false, true, false}},{{false, false, false, false}}
				},
				{
					{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}}
				},
				{
					{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}}
				},
				{
					{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}}
				}
			},
			{

				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				}
			},
			{},
			{

				{
					{{false, false, false, true}},{{false, true, false, false}}
				}
			},
			{

				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				}
			},
			{

				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				}
			},
			{

				{
					{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}}
				},
				{
					{{false, false, false, true}},{{false, false, false, false}},{{false, false, true, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}}
				},
				{
					{{false, false, false, true}},{{false, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}}
				},
				{
					{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}}
				}
			},
			{

				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				}
			},
			{

				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				}
			},
			{

				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				}
			},
			{},
			{

				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				}
			},
			{},
			{

				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				}
			},
			{},
			{

				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				}
			},
			{},
			{

				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				}
			},
			{},
			{

				{
					{{false, false, false, true}},{{false, true, false, false}}
				}
			},
			{

				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				}
			},
			{

				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				}
			},
			{

				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				}
			},
			{

				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				}
			},
			{

				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				}
			},
			{

				{
					{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}}
				},
				{
					{{false, false, false, true}},{{false, false, false, false}},{{false, false, true, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}}
				},
				{
					{{false, false, false, true}},{{false, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}}
				},
				{
					{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}}
				}
			},
			{},
			{

				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				}
			},
			{

				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				}
			},
			{

				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				}
			},
			{

				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, false, false, true}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}},{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				}
			},
			{},
			{},
			{},
			{},
			{

				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				}
			},
			{

				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				}
			},
			{},
			{},
			{

				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				}
			},
			{

				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{true, false, false, false}},{{true, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				}
			},
			{},
			{},
			{

				{
					{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, true}},{{false, true, false, false}},{{false, false, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}}
				},
				{
					{{false, false, false, true}},{{false, false, false, false}},{{false, false, true, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}}
				},
				{
					{{false, false, false, true}},{{false, false, false, false}},{{false, true, false, false}},{{false, false, false, false}}
				},
				{
					{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}}
				},
				{
					{{false, false, false, true}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, false, false, false}},{{false, true, false, false}}
				}
			},
			{},
			{},
			{},
			{},
			{},
			{},
			{}
		}};

		assert(!table[idx].empty());
		return table[idx];
	}


	const std::vector<std::vector<Vector4i>>& CutTable::get_face_id_conf(const int idx) {
		static const std::array<std::vector<std::vector<Vector4i>>, 64> table= {{

			{

				{
					Vector4i(0, 1, 2, 3)
				}
			},
			{

				{
					Vector4i(1, 3, 2, -1),Vector4i(0, -1, 2, 3)
				}
			},
			{

				{
					Vector4i(2, 3, 0, -1),Vector4i(1, -1, 0, 3)
				}
			},
			{

				{
					Vector4i(-1, 3, 2, 0),Vector4i(-1, -1, 3, 2),Vector4i(1, 0, 3, -1)
				},
				{
					Vector4i(-1, 3, 2, 0),Vector4i(0, -1, 3, -1),Vector4i(1, -1, 3, 2)
				}
			},
			{

				{
					Vector4i(0, 3, 1, -1),Vector4i(2, -1, 1, 3)
				}
			},
			{

				{
					Vector4i(-1, 1, 2, 3),Vector4i(-1, 2, 3, -1),Vector4i(0, -1, 3, 1)
				},
				{
					Vector4i(-1, 1, 2, 3),Vector4i(1, -1, 3, -1),Vector4i(0, 2, 3, -1)
				}
			},
			{

				{
					Vector4i(-1, 3, 0, 1),Vector4i(-1, -1, 3, 0),Vector4i(2, 1, 3, -1)
				},
				{
					Vector4i(-1, 3, 0, 1),Vector4i(1, -1, 3, -1),Vector4i(2, -1, 3, 0)
				}
			},
			{},
			{

				{
					Vector4i(3, 2, 1, -1),Vector4i(0, -1, 1, 2)
				}
			},
			{

				{
					Vector4i(-1, 2, 3, 1),Vector4i(-1, -1, 2, 3),Vector4i(0, 1, 2, -1)
				},
				{
					Vector4i(-1, 2, 3, 1),Vector4i(1, -1, 2, -1),Vector4i(0, -1, 2, 3)
				}
			},
			{

				{
					Vector4i(2, 3, -1, -1),Vector4i(1, -1, -1, 3),Vector4i(2, -1, -1, 0),Vector4i(1, 0, -1, -1)
				}
			},
			{

				{
					Vector4i(-1, 3, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 1),Vector4i(-1, 1, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 2)
				},
				{
					Vector4i(-1, 3, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 1),Vector4i(-1, 1, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 2)
				},
				{
					Vector4i(-1, 3, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 1),Vector4i(-1, 1, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 2)
				},
				{
					Vector4i(-1, 3, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 1),Vector4i(-1, 1, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 2)
				}
			},
			{

				{
					Vector4i(-1, 2, 3, 1),Vector4i(-1, 3, 1, -1),Vector4i(0, -1, 1, 2)
				},
				{
					Vector4i(-1, 2, 3, 1),Vector4i(2, -1, 1, -1),Vector4i(0, 3, 1, -1)
				}
			},
			{

				{
					Vector4i(1, -1, 3, 2),Vector4i(1, -1, 3, -1),Vector4i(0, 2, -1, 1),Vector4i(-1, 2, 3, -1)
				},
				{
					Vector4i(1, -1, 3, 2),Vector4i(1, -1, 3, -1),Vector4i(-1, 2, -1, 1),Vector4i(0, 2, 3, -1)
				},
				{
					Vector4i(1, -1, 3, 2),Vector4i(1, -1, 3, -1),Vector4i(0, 2, -1, 1),Vector4i(-1, 2, 3, -1)
				},
				{
					Vector4i(1, -1, 3, 2),Vector4i(1, -1, 3, -1),Vector4i(-1, 2, -1, 1),Vector4i(0, 2, 3, -1)
				},
				{
					Vector4i(1, -1, 3, 2),Vector4i(1, 0, 3, -1),Vector4i(1, 2, -1, -1),Vector4i(-1, 2, 3, -1)
				},
				{
					Vector4i(1, -1, 3, 2),Vector4i(1, 0, 3, -1),Vector4i(1, 2, -1, -1),Vector4i(-1, 2, 3, -1)
				},
				{
					Vector4i(1, -1, 3, 2),Vector4i(-1, -1, -1, 2),Vector4i(-1, 2, -1, -1),Vector4i(-1, 1, -1, -1),Vector4i(-1, -1, -1, 1),Vector4i(-1, 3, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, -1)
				},
				{
					Vector4i(1, -1, 3, 2),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 2),Vector4i(-1, 1, -1, -1),Vector4i(-1, -1, -1, 1),Vector4i(-1, 3, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, -1)
				}
			},
			{

				{
					Vector4i(-1, 3, -1, -1),Vector4i(-1, -1, -1, 1),Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, 1, -1, -1),Vector4i(-1, 1, -1, -1)
				},
				{
					Vector4i(-1, 3, -1, -1),Vector4i(-1, -1, -1, 1),Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, 1, -1, -1),Vector4i(-1, 1, -1, -1)
				},
				{
					Vector4i(-1, 3, -1, -1),Vector4i(-1, -1, -1, 1),Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, 1, -1, -1),Vector4i(-1, 1, -1, -1)
				},
				{
					Vector4i(-1, 3, -1, -1),Vector4i(-1, -1, -1, 1),Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, 1, -1, -1),Vector4i(-1, 1, -1, -1)
				}
			},
			{},
			{

				{
					Vector4i(3, 0, 2, -1),Vector4i(1, -1, 2, 0)
				}
			},
			{

				{
					Vector4i(-1, 0, 3, 2),Vector4i(-1, 3, 2, -1),Vector4i(1, -1, 2, 0)
				},
				{
					Vector4i(-1, 0, 3, 2),Vector4i(0, -1, 2, -1),Vector4i(1, 3, 2, -1)
				}
			},
			{

				{
					Vector4i(-1, 0, 3, 2),Vector4i(-1, -1, 0, 3),Vector4i(1, 2, 0, -1)
				},
				{
					Vector4i(-1, 0, 3, 2),Vector4i(2, -1, 0, -1),Vector4i(1, -1, 0, 3)
				}
			},
			{

				{
					Vector4i(0, 2, 3, -1),Vector4i(0, -1, 3, -1),Vector4i(1, 0, -1, 2),Vector4i(-1, -1, 3, 2)
				},
				{
					Vector4i(0, 2, 3, -1),Vector4i(0, -1, 3, -1),Vector4i(-1, 0, -1, 2),Vector4i(1, -1, 3, 2)
				},
				{
					Vector4i(0, 2, 3, -1),Vector4i(0, -1, 3, -1),Vector4i(1, 0, -1, 2),Vector4i(-1, -1, 3, 2)
				},
				{
					Vector4i(0, 2, 3, -1),Vector4i(0, -1, 3, -1),Vector4i(-1, 0, -1, 2),Vector4i(1, -1, 3, 2)
				},
				{
					Vector4i(0, 2, 3, -1),Vector4i(0, -1, 3, 1),Vector4i(0, -1, -1, 2),Vector4i(-1, -1, 3, 2)
				},
				{
					Vector4i(0, 2, 3, -1),Vector4i(0, -1, 3, 1),Vector4i(0, -1, -1, 2),Vector4i(-1, -1, 3, 2)
				},
				{
					Vector4i(0, 2, 3, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 0),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, -1)
				},
				{
					Vector4i(0, 2, 3, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, -1)
				}
			},
			{

				{
					Vector4i(0, 3, -1, -1),Vector4i(2, -1, -1, 3),Vector4i(0, -1, -1, 1),Vector4i(2, 1, -1, -1)
				}
			},
			{

				{
					Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 1),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, 3, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 2, -1, -1)
				},
				{
					Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 1),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, 3, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 2, -1, -1)
				},
				{
					Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 1),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, 3, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 2, -1, -1)
				},
				{
					Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 1),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, 3, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 2, -1, -1)
				}
			},
			{

				{
					Vector4i(-1, 3, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, 2),Vector4i(-1, 2, -1, -1),Vector4i(-1, 1, -1, -1),Vector4i(-1, 1, -1, -1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 0)
				},
				{
					Vector4i(-1, 3, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, 2),Vector4i(-1, 2, -1, -1),Vector4i(-1, 1, -1, -1),Vector4i(-1, 1, -1, -1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 0)
				},
				{
					Vector4i(-1, 3, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, 2),Vector4i(-1, 2, -1, -1),Vector4i(-1, 1, -1, -1),Vector4i(-1, 1, -1, -1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 0)
				},
				{
					Vector4i(-1, 3, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, 2),Vector4i(-1, 2, -1, -1),Vector4i(-1, 1, -1, -1),Vector4i(-1, 1, -1, -1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 0)
				}
			},
			{},
			{

				{
					Vector4i(-1, 2, 1, 0),Vector4i(-1, -1, 2, 1),Vector4i(3, 0, 2, -1)
				},
				{
					Vector4i(-1, 2, 1, 0),Vector4i(0, -1, 2, -1),Vector4i(3, -1, 2, 1)
				}
			},
			{},
			{

				{
					Vector4i(-1, 0, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 2)
				},
				{
					Vector4i(-1, 0, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 2)
				},
				{
					Vector4i(-1, 0, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 2)
				},
				{
					Vector4i(-1, 0, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 2)
				}
			},
			{},
			{

				{
					Vector4i(-1, -1, -1, 1),Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 3),Vector4i(-1, 3, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, 1, -1, -1),Vector4i(-1, 1, -1, -1),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 2)
				},
				{
					Vector4i(-1, -1, -1, 1),Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 3),Vector4i(-1, 3, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, 1, -1, -1),Vector4i(-1, 1, -1, -1),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 2)
				},
				{
					Vector4i(-1, -1, -1, 1),Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 3),Vector4i(-1, 3, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, 1, -1, -1),Vector4i(-1, 1, -1, -1),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 2)
				},
				{
					Vector4i(-1, -1, -1, 1),Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 3),Vector4i(-1, 3, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, 1, -1, -1),Vector4i(-1, 1, -1, -1),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 2)
				}
			},
			{},
			{

				{
					Vector4i(-1, 0, 3, -1),Vector4i(1, -1, -1, 2),Vector4i(-1, -1, 3, 2),Vector4i(-1, 2, 1, -1),Vector4i(3, -1, -1, 0),Vector4i(-1, -1, 1, 0)
				},
				{
					Vector4i(-1, 0, 3, -1),Vector4i(1, -1, -1, 2),Vector4i(-1, -1, 3, 2),Vector4i(-1, 2, 1, -1),Vector4i(-1, -1, -1, 0),Vector4i(3, -1, 1, 0)
				},
				{
					Vector4i(-1, 0, 3, -1),Vector4i(1, -1, -1, 2),Vector4i(-1, -1, 3, 2),Vector4i(-1, -1, 1, -1),Vector4i(3, -1, -1, 0),Vector4i(-1, 2, 1, 0)
				},
				{
					Vector4i(-1, 0, 3, -1),Vector4i(1, -1, -1, 2),Vector4i(-1, -1, 3, 2),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 3),Vector4i(-1, 2, -1, -1)
				},
				{
					Vector4i(-1, 0, 3, -1),Vector4i(-1, -1, -1, 2),Vector4i(1, -1, 3, 2),Vector4i(-1, 2, 1, -1),Vector4i(3, -1, -1, 0),Vector4i(-1, -1, 1, 0)
				},
				{
					Vector4i(-1, 0, 3, -1),Vector4i(-1, -1, -1, 2),Vector4i(1, -1, 3, 2),Vector4i(-1, 2, 1, -1),Vector4i(-1, -1, -1, 0),Vector4i(3, -1, 1, 0)
				},
				{
					Vector4i(-1, 0, 3, -1),Vector4i(-1, -1, -1, 2),Vector4i(1, -1, 3, 2),Vector4i(-1, -1, 1, -1),Vector4i(3, -1, -1, 0),Vector4i(-1, 2, 1, 0)
				},
				{
					Vector4i(-1, 0, 3, -1),Vector4i(-1, -1, -1, 2),Vector4i(1, -1, 3, 2),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 3),Vector4i(-1, 2, -1, -1)
				},
				{
					Vector4i(-1, -1, 3, -1),Vector4i(1, -1, -1, 2),Vector4i(-1, 0, 3, 2),Vector4i(-1, 2, 1, -1),Vector4i(3, -1, -1, 0),Vector4i(-1, -1, 1, 0)
				},
				{
					Vector4i(-1, -1, 3, -1),Vector4i(1, -1, -1, 2),Vector4i(-1, 0, 3, 2),Vector4i(-1, 2, 1, -1),Vector4i(-1, -1, -1, 0),Vector4i(3, -1, 1, 0)
				},
				{
					Vector4i(-1, -1, 3, -1),Vector4i(1, -1, -1, 2),Vector4i(-1, 0, 3, 2),Vector4i(-1, -1, 1, -1),Vector4i(3, -1, -1, 0),Vector4i(-1, 2, 1, 0)
				},
				{
					Vector4i(-1, -1, 3, -1),Vector4i(1, -1, -1, 2),Vector4i(-1, 0, 3, 2),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 3),Vector4i(-1, 2, -1, -1)
				},
				{
					Vector4i(-1, -1, 3, -1),Vector4i(-1, -1, 0, 2),Vector4i(1, -1, 3, 2),Vector4i(-1, -1, 1, -1),Vector4i(-1, -1, 2, 0),Vector4i(3, -1, 1, 0)
				},
				{
					Vector4i(-1, -1, 3, -1),Vector4i(-1, -1, 0, 2),Vector4i(1, -1, 3, 2),Vector4i(-1, -1, 1, 3),Vector4i(-1, -1, -1, 0),Vector4i(-1, 2, 1, 0)
				},
				{
					Vector4i(-1, -1, 3, -1),Vector4i(-1, -1, 0, 2),Vector4i(1, -1, 3, 2),Vector4i(-1, -1, 1, 3),Vector4i(-1, 2, -1, 0),Vector4i(-1, -1, 1, 0)
				},
				{
					Vector4i(-1, -1, 3, -1),Vector4i(-1, -1, 0, 2),Vector4i(1, -1, 3, 2),Vector4i(-1, 0, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 3),Vector4i(-1, 2, -1, -1)
				},
				{
					Vector4i(-1, -1, 3, 1),Vector4i(-1, -1, -1, 2),Vector4i(-1, 0, 3, 2),Vector4i(-1, -1, 1, -1),Vector4i(-1, -1, 2, 0),Vector4i(3, -1, 1, 0)
				},
				{
					Vector4i(-1, -1, 3, 1),Vector4i(-1, -1, -1, 2),Vector4i(-1, 0, 3, 2),Vector4i(-1, -1, 1, 3),Vector4i(-1, -1, -1, 0),Vector4i(-1, 2, 1, 0)
				},
				{
					Vector4i(-1, -1, 3, 1),Vector4i(-1, -1, -1, 2),Vector4i(-1, 0, 3, 2),Vector4i(-1, -1, 1, 3),Vector4i(-1, 2, -1, 0),Vector4i(-1, -1, 1, 0)
				},
				{
					Vector4i(-1, -1, 3, 1),Vector4i(-1, -1, -1, 2),Vector4i(-1, 0, 3, 2),Vector4i(-1, 0, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 3),Vector4i(-1, 2, -1, -1)
				},
				{
					Vector4i(-1, -1, 3, 1),Vector4i(-1, 0, -1, 2),Vector4i(-1, -1, 3, 2),Vector4i(-1, -1, 1, -1),Vector4i(-1, -1, 2, 0),Vector4i(3, -1, 1, 0)
				},
				{
					Vector4i(-1, -1, 3, 1),Vector4i(-1, 0, -1, 2),Vector4i(-1, -1, 3, 2),Vector4i(-1, -1, 1, 3),Vector4i(-1, -1, -1, 0),Vector4i(-1, 2, 1, 0)
				},
				{
					Vector4i(-1, -1, 3, 1),Vector4i(-1, 0, -1, 2),Vector4i(-1, -1, 3, 2),Vector4i(-1, -1, 1, 3),Vector4i(-1, 2, -1, 0),Vector4i(-1, -1, 1, 0)
				},
				{
					Vector4i(-1, -1, 3, 1),Vector4i(-1, 0, -1, 2),Vector4i(-1, -1, 3, 2),Vector4i(-1, 0, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 3),Vector4i(-1, 2, -1, -1)
				},
				{
					Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 1),Vector4i(-1, 0, -1, -1),Vector4i(-1, 2, 1, -1),Vector4i(3, -1, -1, 0),Vector4i(-1, -1, 1, 0)
				},
				{
					Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 1),Vector4i(-1, 0, -1, -1),Vector4i(-1, 2, 1, -1),Vector4i(-1, -1, -1, 0),Vector4i(3, -1, 1, 0)
				},
				{
					Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 1),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, 1, -1),Vector4i(3, -1, -1, 0),Vector4i(-1, 2, 1, 0)
				},
				{
					Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 1),Vector4i(-1, 0, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 3),Vector4i(-1, 2, -1, -1)
				},
				{
					Vector4i(-1, 2, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 1),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, 1, -1),Vector4i(-1, -1, 2, 0),Vector4i(3, -1, 1, 0)
				},
				{
					Vector4i(-1, 2, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 1),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, 1, 3),Vector4i(-1, -1, -1, 0),Vector4i(-1, 2, 1, 0)
				},
				{
					Vector4i(-1, 2, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 1),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, 1, 3),Vector4i(-1, 2, -1, 0),Vector4i(-1, -1, 1, 0)
				},
				{
					Vector4i(-1, 2, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 1),Vector4i(-1, 0, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 3),Vector4i(-1, 2, -1, -1)
				}
			},
			{},
			{

				{
					Vector4i(3, 1, 0, -1),Vector4i(2, -1, 0, 1)
				}
			},
			{

				{
					Vector4i(1, 3, -1, -1),Vector4i(0, -1, -1, 3),Vector4i(1, -1, -1, 2),Vector4i(0, 2, -1, -1)
				}
			},
			{

				{
					Vector4i(-1, 1, 3, 0),Vector4i(-1, 3, 0, -1),Vector4i(2, -1, 0, 1)
				},
				{
					Vector4i(-1, 1, 3, 0),Vector4i(1, -1, 0, -1),Vector4i(2, 3, 0, -1)
				}
			},
			{

				{
					Vector4i(-1, 3, -1, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 2),Vector4i(-1, 1, -1, -1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, 0, -1, -1),Vector4i(-1, 0, -1, -1)
				},
				{
					Vector4i(-1, 3, -1, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 2),Vector4i(-1, 1, -1, -1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, 0, -1, -1),Vector4i(-1, 0, -1, -1)
				},
				{
					Vector4i(-1, 3, -1, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 2),Vector4i(-1, 1, -1, -1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, 0, -1, -1),Vector4i(-1, 0, -1, -1)
				},
				{
					Vector4i(-1, 3, -1, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 2),Vector4i(-1, 1, -1, -1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, 0, -1, -1),Vector4i(-1, 0, -1, -1)
				}
			},
			{

				{
					Vector4i(-1, 1, 3, 0),Vector4i(-1, -1, 1, 3),Vector4i(2, 0, 1, -1)
				},
				{
					Vector4i(-1, 1, 3, 0),Vector4i(0, -1, 1, -1),Vector4i(2, -1, 1, 3)
				}
			},
			{

				{
					Vector4i(-1, -1, -1, 3),Vector4i(-1, 1, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, 0, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 1)
				},
				{
					Vector4i(-1, -1, -1, 3),Vector4i(-1, 1, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, 0, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 1)
				},
				{
					Vector4i(-1, -1, -1, 3),Vector4i(-1, 1, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, 0, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 1)
				},
				{
					Vector4i(-1, -1, -1, 3),Vector4i(-1, 1, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, 0, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 1)
				}
			},
			{

				{
					Vector4i(1, 0, 3, -1),Vector4i(1, -1, 3, -1),Vector4i(2, 1, -1, 0),Vector4i(-1, -1, 3, 0)
				},
				{
					Vector4i(1, 0, 3, -1),Vector4i(1, -1, 3, -1),Vector4i(-1, 1, -1, 0),Vector4i(2, -1, 3, 0)
				},
				{
					Vector4i(1, 0, 3, -1),Vector4i(1, -1, 3, -1),Vector4i(2, 1, -1, 0),Vector4i(-1, -1, 3, 0)
				},
				{
					Vector4i(1, 0, 3, -1),Vector4i(1, -1, 3, -1),Vector4i(-1, 1, -1, 0),Vector4i(2, -1, 3, 0)
				},
				{
					Vector4i(1, 0, 3, -1),Vector4i(1, -1, 3, 2),Vector4i(1, -1, -1, 0),Vector4i(-1, -1, 3, 0)
				},
				{
					Vector4i(1, 0, 3, -1),Vector4i(1, -1, 3, 2),Vector4i(1, -1, -1, 0),Vector4i(-1, -1, 3, 0)
				},
				{
					Vector4i(1, 0, 3, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 1),Vector4i(-1, 1, -1, -1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, -1)
				},
				{
					Vector4i(1, 0, 3, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, 1),Vector4i(-1, 1, -1, -1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, -1)
				}
			},
			{},
			{

				{
					Vector4i(-1, 0, 2, 1),Vector4i(-1, 2, 1, -1),Vector4i(3, -1, 1, 0)
				},
				{
					Vector4i(-1, 0, 2, 1),Vector4i(0, -1, 1, -1),Vector4i(3, 2, 1, -1)
				}
			},
			{

				{
					Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 2),Vector4i(-1, 1, -1, -1),Vector4i(-1, 1, -1, -1)
				},
				{
					Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 2),Vector4i(-1, 1, -1, -1),Vector4i(-1, 1, -1, -1)
				},
				{
					Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 2),Vector4i(-1, 1, -1, -1),Vector4i(-1, 1, -1, -1)
				},
				{
					Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 2),Vector4i(-1, 1, -1, -1),Vector4i(-1, 1, -1, -1)
				}
			},
			{

				{
					Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 3),Vector4i(-1, 3, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, 1, -1, -1),Vector4i(-1, 1, -1, -1)
				},
				{
					Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 3),Vector4i(-1, 3, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, 1, -1, -1),Vector4i(-1, 1, -1, -1)
				},
				{
					Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 3),Vector4i(-1, 3, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, 1, -1, -1),Vector4i(-1, 1, -1, -1)
				},
				{
					Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 3),Vector4i(-1, 3, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, 1, -1, -1),Vector4i(-1, 1, -1, -1)
				}
			},
			{

				{
					Vector4i(-1, -1, 3, 0),Vector4i(2, 1, -1, -1),Vector4i(-1, 1, 3, -1),Vector4i(-1, -1, 2, 1),Vector4i(3, 0, -1, -1),Vector4i(-1, 0, 2, -1)
				},
				{
					Vector4i(-1, -1, 3, 0),Vector4i(2, 1, -1, -1),Vector4i(-1, 1, 3, -1),Vector4i(-1, -1, 2, 1),Vector4i(-1, 0, -1, -1),Vector4i(3, 0, 2, -1)
				},
				{
					Vector4i(-1, -1, 3, 0),Vector4i(2, 1, -1, -1),Vector4i(-1, 1, 3, -1),Vector4i(-1, -1, 2, -1),Vector4i(3, 0, -1, -1),Vector4i(-1, 0, 2, 1)
				},
				{
					Vector4i(-1, -1, 3, 0),Vector4i(2, 1, -1, -1),Vector4i(-1, 1, 3, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, -1, -1, 1)
				},
				{
					Vector4i(-1, -1, 3, 0),Vector4i(-1, 1, -1, -1),Vector4i(2, 1, 3, -1),Vector4i(-1, -1, 2, 1),Vector4i(3, 0, -1, -1),Vector4i(-1, 0, 2, -1)
				},
				{
					Vector4i(-1, -1, 3, 0),Vector4i(-1, 1, -1, -1),Vector4i(2, 1, 3, -1),Vector4i(-1, -1, 2, 1),Vector4i(-1, 0, -1, -1),Vector4i(3, 0, 2, -1)
				},
				{
					Vector4i(-1, -1, 3, 0),Vector4i(-1, 1, -1, -1),Vector4i(2, 1, 3, -1),Vector4i(-1, -1, 2, -1),Vector4i(3, 0, -1, -1),Vector4i(-1, 0, 2, 1)
				},
				{
					Vector4i(-1, -1, 3, 0),Vector4i(-1, 1, -1, -1),Vector4i(2, 1, 3, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, -1, -1, 1)
				},
				{
					Vector4i(-1, -1, 3, -1),Vector4i(2, 1, -1, -1),Vector4i(-1, 1, 3, 0),Vector4i(-1, -1, 2, 1),Vector4i(3, 0, -1, -1),Vector4i(-1, 0, 2, -1)
				},
				{
					Vector4i(-1, -1, 3, -1),Vector4i(2, 1, -1, -1),Vector4i(-1, 1, 3, 0),Vector4i(-1, -1, 2, 1),Vector4i(-1, 0, -1, -1),Vector4i(3, 0, 2, -1)
				},
				{
					Vector4i(-1, -1, 3, -1),Vector4i(2, 1, -1, -1),Vector4i(-1, 1, 3, 0),Vector4i(-1, -1, 2, -1),Vector4i(3, 0, -1, -1),Vector4i(-1, 0, 2, 1)
				},
				{
					Vector4i(-1, -1, 3, -1),Vector4i(2, 1, -1, -1),Vector4i(-1, 1, 3, 0),Vector4i(-1, -1, -1, 0),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, -1, -1, 1)
				},
				{
					Vector4i(-1, -1, 3, -1),Vector4i(-1, 1, 0, -1),Vector4i(2, 1, 3, -1),Vector4i(-1, -1, 2, -1),Vector4i(-1, 0, 1, -1),Vector4i(3, 0, 2, -1)
				},
				{
					Vector4i(-1, -1, 3, -1),Vector4i(-1, 1, 0, -1),Vector4i(2, 1, 3, -1),Vector4i(-1, 3, 2, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, 0, 2, 1)
				},
				{
					Vector4i(-1, -1, 3, -1),Vector4i(-1, 1, 0, -1),Vector4i(2, 1, 3, -1),Vector4i(-1, 3, 2, -1),Vector4i(-1, 0, -1, 1),Vector4i(-1, 0, 2, -1)
				},
				{
					Vector4i(-1, -1, 3, -1),Vector4i(-1, 1, 0, -1),Vector4i(2, 1, 3, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, -1, -1, 1)
				},
				{
					Vector4i(-1, 2, 3, -1),Vector4i(-1, 1, -1, -1),Vector4i(-1, 1, 3, 0),Vector4i(-1, -1, 2, -1),Vector4i(-1, 0, 1, -1),Vector4i(3, 0, 2, -1)
				},
				{
					Vector4i(-1, 2, 3, -1),Vector4i(-1, 1, -1, -1),Vector4i(-1, 1, 3, 0),Vector4i(-1, 3, 2, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, 0, 2, 1)
				},
				{
					Vector4i(-1, 2, 3, -1),Vector4i(-1, 1, -1, -1),Vector4i(-1, 1, 3, 0),Vector4i(-1, 3, 2, -1),Vector4i(-1, 0, -1, 1),Vector4i(-1, 0, 2, -1)
				},
				{
					Vector4i(-1, 2, 3, -1),Vector4i(-1, 1, -1, -1),Vector4i(-1, 1, 3, 0),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, -1, -1, 1)
				},
				{
					Vector4i(-1, 2, 3, -1),Vector4i(-1, 1, -1, 0),Vector4i(-1, 1, 3, -1),Vector4i(-1, -1, 2, -1),Vector4i(-1, 0, 1, -1),Vector4i(3, 0, 2, -1)
				},
				{
					Vector4i(-1, 2, 3, -1),Vector4i(-1, 1, -1, 0),Vector4i(-1, 1, 3, -1),Vector4i(-1, 3, 2, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, 0, 2, 1)
				},
				{
					Vector4i(-1, 2, 3, -1),Vector4i(-1, 1, -1, 0),Vector4i(-1, 1, 3, -1),Vector4i(-1, 3, 2, -1),Vector4i(-1, 0, -1, 1),Vector4i(-1, 0, 2, -1)
				},
				{
					Vector4i(-1, 2, 3, -1),Vector4i(-1, 1, -1, 0),Vector4i(-1, 1, 3, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, -1, -1, 1)
				},
				{
					Vector4i(-1, -1, -1, 1),Vector4i(-1, 1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, 2, 1),Vector4i(3, 0, -1, -1),Vector4i(-1, 0, 2, -1)
				},
				{
					Vector4i(-1, -1, -1, 1),Vector4i(-1, 1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, 2, 1),Vector4i(-1, 0, -1, -1),Vector4i(3, 0, 2, -1)
				},
				{
					Vector4i(-1, -1, -1, 1),Vector4i(-1, 1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, 2, -1),Vector4i(3, 0, -1, -1),Vector4i(-1, 0, 2, 1)
				},
				{
					Vector4i(-1, -1, -1, 1),Vector4i(-1, 1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 0),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, -1, -1, 1)
				},
				{
					Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, 2, -1),Vector4i(-1, 0, 1, -1),Vector4i(3, 0, 2, -1)
				},
				{
					Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, 3, 2, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, 0, 2, 1)
				},
				{
					Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, 3, 2, -1),Vector4i(-1, 0, -1, 1),Vector4i(-1, 0, 2, -1)
				},
				{
					Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, -1, -1, 1)
				}
			},
			{},
			{},
			{},
			{},
			{

				{
					Vector4i(-1, 0, 2, 1),Vector4i(-1, -1, 0, 2),Vector4i(3, 1, 0, -1)
				},
				{
					Vector4i(-1, 0, 2, 1),Vector4i(1, -1, 0, -1),Vector4i(3, -1, 0, 2)
				}
			},
			{

				{
					Vector4i(-1, -1, -1, 2),Vector4i(-1, 0, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, -1, -1, 1),Vector4i(-1, 1, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 0)
				},
				{
					Vector4i(-1, -1, -1, 2),Vector4i(-1, 0, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, -1, -1, 1),Vector4i(-1, 1, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 0)
				},
				{
					Vector4i(-1, -1, -1, 2),Vector4i(-1, 0, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, -1, -1, 1),Vector4i(-1, 1, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 0)
				},
				{
					Vector4i(-1, -1, -1, 2),Vector4i(-1, 0, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, 3, -1, -1),Vector4i(-1, -1, -1, 1),Vector4i(-1, 1, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 0)
				}
			},
			{},
			{},
			{

				{
					Vector4i(-1, 1, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 0)
				},
				{
					Vector4i(-1, 1, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 0)
				},
				{
					Vector4i(-1, 1, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 0)
				},
				{
					Vector4i(-1, 1, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, 2, -1, -1),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, 0)
				}
			},
			{

				{
					Vector4i(-1, 1, 3, -1),Vector4i(2, -1, -1, 0),Vector4i(-1, -1, 3, 0),Vector4i(-1, 0, 2, -1),Vector4i(3, -1, -1, 1),Vector4i(-1, -1, 2, 1)
				},
				{
					Vector4i(-1, 1, 3, -1),Vector4i(2, -1, -1, 0),Vector4i(-1, -1, 3, 0),Vector4i(-1, 0, 2, -1),Vector4i(-1, -1, -1, 1),Vector4i(3, -1, 2, 1)
				},
				{
					Vector4i(-1, 1, 3, -1),Vector4i(2, -1, -1, 0),Vector4i(-1, -1, 3, 0),Vector4i(-1, -1, 2, -1),Vector4i(3, -1, -1, 1),Vector4i(-1, 0, 2, 1)
				},
				{
					Vector4i(-1, 1, 3, -1),Vector4i(2, -1, -1, 0),Vector4i(-1, -1, 3, 0),Vector4i(-1, 1, -1, -1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 3),Vector4i(-1, 0, -1, -1)
				},
				{
					Vector4i(-1, 1, 3, -1),Vector4i(-1, -1, -1, 0),Vector4i(2, -1, 3, 0),Vector4i(-1, 0, 2, -1),Vector4i(3, -1, -1, 1),Vector4i(-1, -1, 2, 1)
				},
				{
					Vector4i(-1, 1, 3, -1),Vector4i(-1, -1, -1, 0),Vector4i(2, -1, 3, 0),Vector4i(-1, 0, 2, -1),Vector4i(-1, -1, -1, 1),Vector4i(3, -1, 2, 1)
				},
				{
					Vector4i(-1, 1, 3, -1),Vector4i(-1, -1, -1, 0),Vector4i(2, -1, 3, 0),Vector4i(-1, -1, 2, -1),Vector4i(3, -1, -1, 1),Vector4i(-1, 0, 2, 1)
				},
				{
					Vector4i(-1, 1, 3, -1),Vector4i(-1, -1, -1, 0),Vector4i(2, -1, 3, 0),Vector4i(-1, 1, -1, -1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 3),Vector4i(-1, 0, -1, -1)
				},
				{
					Vector4i(-1, -1, 3, -1),Vector4i(2, -1, -1, 0),Vector4i(-1, 1, 3, 0),Vector4i(-1, 0, 2, -1),Vector4i(3, -1, -1, 1),Vector4i(-1, -1, 2, 1)
				},
				{
					Vector4i(-1, -1, 3, -1),Vector4i(2, -1, -1, 0),Vector4i(-1, 1, 3, 0),Vector4i(-1, 0, 2, -1),Vector4i(-1, -1, -1, 1),Vector4i(3, -1, 2, 1)
				},
				{
					Vector4i(-1, -1, 3, -1),Vector4i(2, -1, -1, 0),Vector4i(-1, 1, 3, 0),Vector4i(-1, -1, 2, -1),Vector4i(3, -1, -1, 1),Vector4i(-1, 0, 2, 1)
				},
				{
					Vector4i(-1, -1, 3, -1),Vector4i(2, -1, -1, 0),Vector4i(-1, 1, 3, 0),Vector4i(-1, 1, -1, -1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 3),Vector4i(-1, 0, -1, -1)
				},
				{
					Vector4i(-1, -1, 3, -1),Vector4i(-1, -1, 1, 0),Vector4i(2, -1, 3, 0),Vector4i(-1, -1, 2, -1),Vector4i(-1, -1, 0, 1),Vector4i(3, -1, 2, 1)
				},
				{
					Vector4i(-1, -1, 3, -1),Vector4i(-1, -1, 1, 0),Vector4i(2, -1, 3, 0),Vector4i(-1, -1, 2, 3),Vector4i(-1, -1, -1, 1),Vector4i(-1, 0, 2, 1)
				},
				{
					Vector4i(-1, -1, 3, -1),Vector4i(-1, -1, 1, 0),Vector4i(2, -1, 3, 0),Vector4i(-1, -1, 2, 3),Vector4i(-1, 0, -1, 1),Vector4i(-1, -1, 2, 1)
				},
				{
					Vector4i(-1, -1, 3, -1),Vector4i(-1, -1, 1, 0),Vector4i(2, -1, 3, 0),Vector4i(-1, 1, -1, -1),Vector4i(-1, 1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 3),Vector4i(-1, 0, -1, -1)
				},
				{
					Vector4i(-1, -1, 3, 2),Vector4i(-1, -1, -1, 0),Vector4i(-1, 1, 3, 0),Vector4i(-1, -1, 2, -1),Vector4i(-1, -1, 0, 1),Vector4i(3, -1, 2, 1)
				},
				{
					Vector4i(-1, -1, 3, 2),Vector4i(-1, -1, -1, 0),Vector4i(-1, 1, 3, 0),Vector4i(-1, -1, 2, 3),Vector4i(-1, -1, -1, 1),Vector4i(-1, 0, 2, 1)
				},
				{
					Vector4i(-1, -1, 3, 2),Vector4i(-1, -1, -1, 0),Vector4i(-1, 1, 3, 0),Vector4i(-1, -1, 2, 3),Vector4i(-1, 0, -1, 1),Vector4i(-1, -1, 2, 1)
				},
				{
					Vector4i(-1, -1, 3, 2),Vector4i(-1, -1, -1, 0),Vector4i(-1, 1, 3, 0),Vector4i(-1, 1, -1, -1),Vector4i(-1, 1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 3),Vector4i(-1, 0, -1, -1)
				},
				{
					Vector4i(-1, -1, 3, 2),Vector4i(-1, 1, -1, 0),Vector4i(-1, -1, 3, 0),Vector4i(-1, -1, 2, -1),Vector4i(-1, -1, 0, 1),Vector4i(3, -1, 2, 1)
				},
				{
					Vector4i(-1, -1, 3, 2),Vector4i(-1, 1, -1, 0),Vector4i(-1, -1, 3, 0),Vector4i(-1, -1, 2, 3),Vector4i(-1, -1, -1, 1),Vector4i(-1, 0, 2, 1)
				},
				{
					Vector4i(-1, -1, 3, 2),Vector4i(-1, 1, -1, 0),Vector4i(-1, -1, 3, 0),Vector4i(-1, -1, 2, 3),Vector4i(-1, 0, -1, 1),Vector4i(-1, -1, 2, 1)
				},
				{
					Vector4i(-1, -1, 3, 2),Vector4i(-1, 1, -1, 0),Vector4i(-1, -1, 3, 0),Vector4i(-1, 1, -1, -1),Vector4i(-1, 1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 3),Vector4i(-1, 0, -1, -1)
				},
				{
					Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 2),Vector4i(-1, 1, -1, -1),Vector4i(-1, 0, 2, -1),Vector4i(3, -1, -1, 1),Vector4i(-1, -1, 2, 1)
				},
				{
					Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 2),Vector4i(-1, 1, -1, -1),Vector4i(-1, 0, 2, -1),Vector4i(-1, -1, -1, 1),Vector4i(3, -1, 2, 1)
				},
				{
					Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 2),Vector4i(-1, 1, -1, -1),Vector4i(-1, -1, 2, -1),Vector4i(3, -1, -1, 1),Vector4i(-1, 0, 2, 1)
				},
				{
					Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 2),Vector4i(-1, 1, -1, -1),Vector4i(-1, 1, -1, -1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 3),Vector4i(-1, 0, -1, -1)
				},
				{
					Vector4i(-1, 0, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 2),Vector4i(-1, 1, -1, -1),Vector4i(-1, -1, 2, -1),Vector4i(-1, -1, 0, 1),Vector4i(3, -1, 2, 1)
				},
				{
					Vector4i(-1, 0, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 2),Vector4i(-1, 1, -1, -1),Vector4i(-1, -1, 2, 3),Vector4i(-1, -1, -1, 1),Vector4i(-1, 0, 2, 1)
				},
				{
					Vector4i(-1, 0, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 2),Vector4i(-1, 1, -1, -1),Vector4i(-1, -1, 2, 3),Vector4i(-1, 0, -1, 1),Vector4i(-1, -1, 2, 1)
				},
				{
					Vector4i(-1, 0, -1, -1),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, 2),Vector4i(-1, 1, -1, -1),Vector4i(-1, 1, -1, -1),Vector4i(-1, 1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, -1),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 3),Vector4i(-1, 0, -1, -1)
				}
			},
			{},
			{},
			{

				{
					Vector4i(0, 1, 2, -1),Vector4i(0, -1, 2, -1),Vector4i(3, 0, -1, 1),Vector4i(-1, -1, 2, 1)
				},
				{
					Vector4i(0, 1, 2, -1),Vector4i(0, -1, 2, -1),Vector4i(-1, 0, -1, 1),Vector4i(3, -1, 2, 1)
				},
				{
					Vector4i(0, 1, 2, -1),Vector4i(0, -1, 2, -1),Vector4i(3, 0, -1, 1),Vector4i(-1, -1, 2, 1)
				},
				{
					Vector4i(0, 1, 2, -1),Vector4i(0, -1, 2, -1),Vector4i(-1, 0, -1, 1),Vector4i(3, -1, 2, 1)
				},
				{
					Vector4i(0, 1, 2, -1),Vector4i(0, -1, 2, 3),Vector4i(0, -1, -1, 1),Vector4i(-1, -1, 2, 1)
				},
				{
					Vector4i(0, 1, 2, -1),Vector4i(0, -1, 2, 3),Vector4i(0, -1, -1, 1),Vector4i(-1, -1, 2, 1)
				},
				{
					Vector4i(0, 1, 2, -1),Vector4i(-1, 1, -1, -1),Vector4i(-1, -1, -1, 1),Vector4i(-1, -1, -1, 0),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, -1)
				},
				{
					Vector4i(0, 1, 2, -1),Vector4i(-1, 1, -1, -1),Vector4i(-1, 1, -1, -1),Vector4i(-1, -1, -1, 0),Vector4i(-1, 0, -1, -1),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 2),Vector4i(-1, -1, -1, 3),Vector4i(-1, -1, -1, -1)
				}
			},
			{},
			{},
			{},
			{},
			{},
			{},
			{}
		}};

		assert(!table[idx].empty());
		return table[idx];
	}


}
