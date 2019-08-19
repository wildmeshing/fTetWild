cmake_minimum_required(VERSION 3.1)
################################################################################

set(test_sources
	main.cpp
	tests.cpp
	test_predicates.cpp
	test_table.cpp
	# test_subdivision.cpp
	# test_tree.cpp
    # test_envelope.cpp
)

add_executable(unit_tests ${test_sources})

foreach(source IN ITEMS ${test_sources})
    source_group("tests" FILES "${source}")
endforeach()

target_link_libraries(unit_tests PUBLIC FloatTetwild catch warnings::all)

# set(DATA_DIR "${THIRD_PARTY_DIR}/data/")
# target_compile_definitions(unit_tests PUBLIC -DDATA_DIR=\"${DATA_DIR}\")


if(FLOAT_TETWILD_WITH_SANITIZERS)
	add_sanitizers(unit_tests)
endif()

# Register tests
set(PARSE_CATCH_TESTS_ADD_TO_CONFIGURE_DEPENDS ON)
include(Catch)
catch_discover_tests(unit_tests)
