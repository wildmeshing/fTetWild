set(SOURCES
	predicates.c
	Predicates.cpp
	Predicates.hpp

	mesh_AABB.h
	mesh_AABB.cpp


	MshSaver.h
	MshSaver.cpp
	Exception.h

	triangle_triangle_intersection.cpp

	getRSS.c
	get_mem.h
	get_mem.cpp

	bfs_orient.h
	bfs_orient.cpp

	FastWindingNumber.hpp
	FastWindingNumber.cpp

		Rational.h
		)

prepend_current_path(SOURCES)
float_tetwild_copy_headers(${SOURCES})
float_tetwild_set_source_group(${SOURCES})
target_sources(FloatTetwild PRIVATE ${SOURCES})
