set(SOURCES
        AABBWrapper.h
        AABBWrapper.cpp


        auto_table.hpp
        auto_table.cpp

	LocalOperations.h
	LocalOperations.cpp
	EdgeSplitting.h
	EdgeSplitting.cpp
	EdgeCollapsing.h
	EdgeCollapsing.cpp
        EdgeSwapping.h
        EdgeSwapping.cpp
        VertexSmoothing.h
        VertexSmoothing.cpp

        MeshImprovement.h
        MeshImprovement.cpp



        FloatTetDelaunay.h
        FloatTetDelaunay.cpp
	FloatTetCutting.h
	FloatTetCutting.cpp
	FloatTetCuttingCheck.h
	FloatTetCuttingCheck.cpp
	FloatTetCuttingParallel.h
	FloatTetCuttingParallel.cpp

	intersections.h
	intersections.cpp

	Simplification.h
	Simplification.cpp

        Mesh.cpp
        Mesh.hpp

        MeshIO.hpp
        MeshIO.cpp

        Logger.hpp
        Logger.cpp

        Parameters.h

        Types.hpp

        Statistics.h
)

prepend_current_path(SOURCES)
float_tetwild_copy_headers(${SOURCES})
float_tetwild_set_source_group(${SOURCES})
target_sources(FloatTetwild PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################

add_subdirectory(external)
