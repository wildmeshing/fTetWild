// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include <floattetwild/Mesh.hpp>
#include <floattetwild/LocalOperations.h>

#include <geogram/mesh/mesh_partition.h>


#include <numeric>

namespace floatTetWild {

	void Mesh::one_ring_vertex_coloring(std::vector<Scalar> &colors) const
	{
		colors.resize(tet_vertices.size());
		std::fill(colors.begin(), colors.end(), -1);
		colors[0] = 0;

		std::vector<bool> available(tet_vertices.size(), true);

		std::vector<int> ring;

		for(int i = 1; i < tet_vertices.size(); ++i)
		{
			const auto &v = tet_vertices[i];
			if(v.is_removed)
				continue;

			ring.clear();
			for(const auto &t : v.conn_tets)
			{
				for(int j = 0; j < 4; ++j)
					ring.push_back(tets[t][j]);

			}
			vector_unique(ring);

			for(const auto n : ring)
			{
				if (colors[n] != -1)
					available[colors[n]] = false;
			}

			int first_available_col;
			for (first_available_col = 0; first_available_col < available.size(); first_available_col++){
				if (available[first_available_col])
					break;
			}

			assert(available[first_available_col]);

			colors[i] = first_available_col;

			for(const auto n : ring)
			{
				if (colors[n] != -1)
					available[colors[n]] = true;
			}

		}
	}

	void Mesh::one_ring_vertex_sets(const int threshold, std::vector<std::vector<int>> &concurrent_sets, std::vector<int> &serial_set) const
	{
		std::vector<Scalar> colors;
		one_ring_vertex_coloring(colors);
		int max_c = -1;
		for(const auto c : colors)
			max_c = std::max(max_c, int(c));

		concurrent_sets.clear();
		concurrent_sets.resize(max_c+1);
		serial_set.clear();

		for(size_t i = 0; i < colors.size(); ++i)
		{
			const int col = colors[i];
			//removed vertex
			if(col < 0)
				serial_set.push_back(i);
			else
				concurrent_sets[col].push_back(i);
		}



		for(int i = concurrent_sets.size() -1; i >= 0; --i)
		{
			if(concurrent_sets[i].size() < threshold)
			{
				serial_set.insert(serial_set.end(), concurrent_sets[i].begin(), concurrent_sets[i].end());
				concurrent_sets.erase(concurrent_sets.begin()+i);
			}
		}
	}


	void Mesh::one_ring_edge_set(const std::vector<std::array<int, 2>> &edges, const std::vector<bool>& v_is_removed, const std::vector<bool>& f_is_removed,
	        const std::vector<std::unordered_set<int>>& conn_fs, const std::vector<Vector3>& input_vertices, std::vector<int> &safe_set)
	{
		// std::vector<int> indices(edges.size());
		// std::vector<Scalar> weights(edges.size());

		// for(size_t e_id = 0; e_id < edges.size(); ++e_id)
		// {
		// 	indices[e_id] = e_id;
		// 	weights[e_id] = (input_vertices[edges[e_id][0]] - input_vertices[edges[e_id][1]]).squaredNorm();
		// }

		// std::sort(indices.begin(), indices.end(), [&](const int i, const int j){ return weights[i] < weights[j]; });
		std::vector<int> indices(edges.size());
		std::iota(std::begin(indices), std::end(indices), 0);
		floatTetWild::Random::shuffle(indices);

		std::vector<bool> unsafe_face(f_is_removed.size(), false);
		safe_set.clear();
		for(const int e_id : indices)
		{
			const auto e = edges[e_id];
			if(v_is_removed[e[0]] || v_is_removed[e[1]])
				continue;

			bool ok = true;


			for(const int f : conn_fs[e[0]])
			{
				if(f_is_removed[f])
					continue;

				if(unsafe_face[f])
				{
					ok=false;
					break;
				}
			}
			if(!ok)
				continue;
			for(const int f : conn_fs[e[1]])
			{
				if(f_is_removed[f])
					continue;

				if(unsafe_face[f])
				{
					ok=false;
					break;
				}
			}
			if(!ok)
				continue;

			safe_set.push_back(e_id);

			for(const int f : conn_fs[e[0]])
			{
				if(f_is_removed[f])
					continue;

				assert(!unsafe_face[f]);
				unsafe_face[f]=true;
			}
			for(const int f : conn_fs[e[1]])
			{
				if(f_is_removed[f])
					continue;

				// assert(!unsafe_face[f]);
				unsafe_face[f]=true;
			}
		}
	}


	// void Mesh::one_ring_edge_set(const std::vector<std::array<int, 2>> &edges, std::vector<int> &safe_set) const
	// {
	// 	// std::vector<int> indices(edges.size());
	// 	// std::vector<Scalar> weights(edges.size());

	// 	// for(size_t e_id = 0; e_id < edges.size(); ++e_id)
	// 	// {
	// 	// 	indices[e_id] = e_id;
	// 	// 	weights[e_id] = (tet_vertices[edges[e_id][0]].pos - tet_vertices[edges[e_id][1]].pos).squaredNorm();
	// 	// }
	// 	// std::sort(indices.begin(), indices.end(), [&](const int i, const int j){ return weights[i] < weights[j]; });

	// 	std::vector<int> indices(edges.size());
	// 	std::iota(std::begin(indices), std::end(indices), 0);
	// 	std::random_shuffle(indices.begin(), indices.end());

	// 	std::vector<bool> unsafe_tets(tets.size(), false);
	// 	safe_set.clear();
	// 	for(const int e_id : indices)
	// 	{
	// 		const auto &e = edges[e_id];
	// 		if(tet_vertices[e[0]].is_removed || tet_vertices[e[1]].is_removed)
	// 			continue;

	// 		bool ok = true;

	// 		for(const auto &t : tet_vertices[e[0]].conn_tets)
	// 		{
	// 			if(tets[t].is_removed)
	// 				continue;

	// 			if(unsafe_tets[t])
	// 			{
	// 				ok=false;
	// 				break;
	// 			}
	// 		}
	// 		if(!ok)
	// 			continue;
	// 		for(const auto &t : tet_vertices[e[1]].conn_tets)
	// 		{
	// 			if(tets[t].is_removed)
	// 				continue;

	// 			if(unsafe_tets[t])
	// 			{
	// 				ok=false;
	// 				break;
	// 			}
	// 		}
	// 		if(!ok)
	// 			continue;

	// 		safe_set.push_back(e_id);

	// 		for(const auto &t : tet_vertices[e[0]].conn_tets)
	// 		{
	// 			if(tets[t].is_removed)
	// 				continue;

	// 			assert(!unsafe_tets[t]);
	// 			unsafe_tets[t]=true;
	// 		}
	// 		for(const auto &t : tet_vertices[e[1]].conn_tets)
	// 		{
	// 			if(tets[t].is_removed)
	// 				continue;

	// 			// assert(!unsafe_tets[t]);
	// 			unsafe_tets[t]=true;
	// 		}
	// 	}
	// }

	void Mesh::partition(const int n_parts, std::vector<std::vector<int>> &tets_id) const
	{
		GEO::Mesh M;
		// Setup vertices
		M.vertices.create_vertices((int) tet_vertices.size());
		for (int i = 0; i < (int) M.vertices.nb(); ++i) {
			GEO::vec3 &p = M.vertices.point(i);
			p[0] = tet_vertices[i][0];
			p[1] = tet_vertices[i][1];
			p[2] = tet_vertices[i][2];
		}

		M.cells.create_tets((int) get_t_num());

		int c = 0;
		for (int index = 0; index < (int) tets.size(); ++index) {
			if(tets[index].is_removed)
				continue;

			for (int lv = 0; lv < 4; ++lv) {
				M.cells.set_vertex(c, lv, tets[index][lv]);
			}

			++c;
		}
		M.cells.connect();

		{
			c = 0;
			GEO::Attribute<int> original_indices(M.cells.attributes(), "indices");
			for (int index = 0; index < (int) tets.size(); ++index) {
				if(tets[index].is_removed)
					continue;
				original_indices[c] = index;
				++c;
			}
		}

		GEO::vector<GEO::index_t> facet_ptr, tet_ptr;
		GEO::mesh_partition(M, GEO::MESH_PARTITION_HILBERT, facet_ptr, tet_ptr, n_parts);
		//GEO::MESH_PARTITION_CONNECTED_COMPONENTS

		GEO::Attribute<int> original_indices(M.cells.attributes(), "indices");

		tets_id.clear();
		tets_id.resize(tet_ptr.size()-1);
		for(int i = 0; i < tet_ptr.size()-1; ++i)
		{
			for(int j = tet_ptr[i]; j < tet_ptr[i+1]; ++j)
				tets_id[i].push_back(original_indices[j]);
		}


		// static const int UNVISITED = -1;

	    //     std::vector<int> new_index(tets.size(), UNVISITED);
	    //     std::stack<int> S;
	    //     int new_cur_index = 0;
	    //     for(int t = 0; t < tets.size(); ++t) {
	    //         if(new_index[t] == UNVISITED) {
	    //             // tet_ptr.push_back(new_cur_index);
	    //             new_index[t] = new_cur_index;
	    //             new_cur_index++;
	    //             S.push(t);
	    //         }
	    //         while(!S.empty()) {
	    //             int t1 = S.top();
	    //             S.pop();

	    //             for(int lf = 0; lf < 4; lf++) {
	    //                 int t2 = tets[t1].opp_t_ids[lf];
	    //                 if(t2 != -1 && new_index[t2] == UNVISITED) {
	    //                     new_index[t2] = new_cur_index;
	    //                     new_cur_index++;
	    //                     S.push(int(t2));
	    //                 }
	    //             }
	    //         }
	    //     }
	    //     assert(new_cur_index == tets.size());
	    //     // tet_ptr.push_back(new_cur_index);
	    //     // Permutation::invert(new_index);
	    //     // M.cells.permute_elements(new_index);
	}

}