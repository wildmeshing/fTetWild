if(FLOAT_TETWILD_WITH_EXACT_ENVELOPE )
	set(SOURCES
	#	predicates.c
		Predicates.cpp
		Predicates.hpp

		mesh_AABB.h
		mesh_AABB.cpp


		MshSaver.h
		MshSaver.cpp
		MshLoader.h
		MshLoader.cpp
		Exception.h

		triangle_triangle_intersection.cpp

		getRSS.c
		get_mem.h
		get_mem.cpp

		bfs_orient.h
		bfs_orient.cpp

#			WindingNumber.h
#		FastWindingNumber.hpp
#		FastWindingNumber.cpp

		Rational.h
	)
else()
	set(SOURCES
#			predicates.c
			Predicates.cpp
			Predicates.hpp

			mesh_AABB.h
			mesh_AABB.cpp


			MshSaver.h
			MshSaver.cpp
			MshLoader.h
			MshLoader.cpp
			Exception.h

			triangle_triangle_intersection.cpp

			getRSS.c
			get_mem.h
			get_mem.cpp

			bfs_orient.h
			bfs_orient.cpp

#			WindingNumber.h
#			FastWindingNumber.hpp
#			FastWindingNumber.cpp

			Rational.h
			)
endif()

prepend_current_path(SOURCES)
float_tetwild_copy_headers(${SOURCES})
float_tetwild_set_source_group(${SOURCES})
target_sources(FloatTetwild PRIVATE ${SOURCES})
