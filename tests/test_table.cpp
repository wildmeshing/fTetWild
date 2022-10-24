////////////////////////////////////////////////////////////////////////////////
#include <floattetwild/auto_table.hpp>
#include <floattetwild/Predicates.hpp>
#include <catch2/catch_all.hpp>
#include <iostream>
////////////////////////////////////////////////////////////////////////////////

using namespace floatTetWild;

static const std::vector<int> valid_indices = { 0, 1, 2, 3, 4, 5, 6, 8, 9, 12, 13, 16, 17, 18, 19, 24, 30, 32, 34, 36, 38, 40, 43, 48, 53, 56};

// TEST_CASE("orientation", "[table]") {
// 	// std::vector<Vector3> vertices = {Vector3(0, 0, 0), Vector3(1, 0, 0), Vector3(0, 1, 0), Vector3(0, 0, 1)};
// 	// for(const auto idx : valid_indices)
// 	// {
// 	// 	const auto &tet_confs = CutTable::get_tet_confs(idx);

// 	// 	for(const auto &tets : tet_confs)
// 	// 	{
// 	// 		for(const auto &tet: tets)
// 	// 		{
// 	// 			std::cout<<tet.transpose()<<std::endl;
// 	// 			const auto orientation = Predicates::orient_3d(vertices[tet(0)], vertices[tet(1)], vertices[tet(2)], vertices[tet(3)]);
// 	// 			REQUIRE(orientation == Predicates::ORI_POSITIVE);
// 	// 		}
// 	// 	}
// 	// }
// }

TEST_CASE("table", "[table]") {

	for(const auto idx : valid_indices)
	{
		const auto &tet_confs = CutTable::get_tet_confs(idx);
		const auto &diag_confs = CutTable::get_diag_confs(idx);

		REQUIRE(tet_confs.size() == diag_confs.size());

		if(tet_confs.size() == 1)
		{
			REQUIRE(diag_confs.front().empty());
		}
		else if(tet_confs.size() == 2)
		{
			for(const auto &ec : diag_confs)
				REQUIRE(ec.size() == 1);
		}
		else if(tet_confs.size() == 8)
		{
			for(const auto &ec : diag_confs)
				REQUIRE(ec.size() == 3);
		}
		else
		{
			// for(const auto &ec : diag_confs){
			// 	REQUIRE(ec.size() == 5);
			// }
		}
	}
}
