# Download Catch2 unit test framework using FetchContent
include(FetchContent)
FetchContent_Declare(
    Catch2
    URL https://github.com/catchorg/Catch2/archive/refs/tags/v3.5.3.tar.gz
    URL_HASH MD5=1f51d817ce81d54b12e87d06e159305f
)
FetchContent_MakeAvailable(Catch2)

# Add Catch2 to module path for test discovery
list(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/extras)

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

target_link_libraries(unit_tests PUBLIC FloatTetwild Catch2::Catch2WithMain warnings::all)

# set(DATA_DIR "${THIRD_PARTY_DIR}/data/")
# target_compile_definitions(unit_tests PUBLIC -DDATA_DIR=\"${DATA_DIR}\")

if(FLOAT_TETWILD_WITH_SANITIZERS)
	add_sanitizers(unit_tests)
endif()

# Register tests
set(PARSE_CATCH_TESTS_ADD_TO_CONFIGURE_DEPENDS ON)
include(Catch)
catch_discover_tests(unit_tests)
