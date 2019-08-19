################################################################################
cmake_minimum_required(VERSION 3.1)
project(FloatTetwild)
################################################################################

# cmake_policy(SET CMP0063 NEW)
# set(CMAKE_CXX_VISIBILITY_PRESET hidden)

# Detects whether this is a top-level project
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
	set(FLOAT_TETWILD_TOPLEVEL_PROJECT ON)
else()
	set(FLOAT_TETWILD_TOPLEVEL_PROJECT OFF)
endif()

if(INPUT_THIRD_PARTY_DIR)
	set(THIRD_PARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/${INPUT_THIRD_PARTY_DIR}/)
else()
	set(THIRD_PARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/)
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH ${THIRD_PARTY_DIR}/Catch2/contrib)

# Color output
include(UseColors)

# Prepend function
include(PrependCurrentPath)

# FloatTetwild utils
include(FloatTetwildUtils)

# Extra warnings
include(Warnings)

# Use C++11/14
include(CXXFeatures)

# Sort projects inside the solution
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Generate position independent code by default
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

################################################################################

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	foreach(config ${CMAKE_CONFIGURATION_TYPES})
		string(TOUPPER ${config} config)
		string(REPLACE /MD /MT CMAKE_C_FLAGS_${config} "${CMAKE_C_FLAGS_${config}}")
		string(REPLACE /MD /MT CMAKE_CXX_FLAGS_${config} "${CMAKE_CXX_FLAGS_${config}}")
	endforeach()
endif()

################################################################################

# FloatTetwild options
option(FLOAT_TETWILD_ENABLE_TBB         "Enable TBB"                    ON)
option(FLOAT_TETWILD_USE_FLOAT          "Use floats instead of double"  OFF)


# Sanitizer options
option(SANITIZE_ADDRESS   "Sanitize Address"       OFF)
option(SANITIZE_MEMORY    "Sanitize Memory"        OFF)
option(SANITIZE_THREAD    "Sanitize Thread"        OFF)
option(SANITIZE_UNDEFINED "Sanitize Undefined"     OFF)

# Options for libigl modules
option(LIBIGL_USE_STATIC_LIBRARY          "Use libigl as static library" OFF)
option(LIBIGL_WITH_COMISO            "Use CoMiso"                   OFF)
option(LIBIGL_WITH_EMBREE            "Use Embree"                   OFF)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"                   OFF)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"                     OFF)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"                    OFF)
option(LIBIGL_WITH_PNG               "Use PNG"                      OFF)
option(LIBIGL_WITH_TETGEN            "Use Tetgen"                   OFF)
option(LIBIGL_WITH_TRIANGLE          "Use Triangle"                 OFF)
option(LIBIGL_WITH_PREDICATES        "Use exact predicates"         OFF)
option(LIBIGL_WITH_XML               "Use XML"                      OFF)
#option(LIBIGL_WITH_PYTHON            "Use Python"                   ${LIBIGL_BUILD_PYTHON})


# Sanitizers
if(FLOAT_TETWILD_WITH_SANITIZERS)
	list(APPEND CMAKE_MODULE_PATH ${THIRD_PARTY_DIR}/sanitizers-cmake/cmake)
endif()

# Setup dependencies
include(FloatTetwildDependencies)

################################################################################
# FloatTetwild library
################################################################################

# add_library() can only be called without any source since CMake 3.11 ...
add_library(${PROJECT_NAME} src/Logger.cpp)

# Public include directory for FloatTetwild
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_BINARY_DIR}/include)

# set(MESH_PATH "${CMAKE_CURRENT_SOURCE_DIR}/tests/")
# target_compile_definitions(FloatTetwild PUBLIC -DFLOAT_TETWILD_MESH_PATH=\"${MESH_PATH}\")

# Dependencies
# target_link_libraries(FloatTetwild PUBLIC clipper::clipper)

# Extra warnings
target_link_libraries(${PROJECT_NAME} PRIVATE warnings::all)

# Use C++14
target_compile_features(${PROJECT_NAME} PUBLIC ${CXX14_FEATURES})
# target_compile_definitions(FloatTetwild PUBLIC -DEIGEN_STACK_ALLOCATION_LIMIT=8388608)

if(FLOAT_TETWILD_WITH_SANITIZERS)
	add_sanitizers(${PROJECT_NAME})
endif()

################################################################################
# FloatTetwild binary
################################################################################

# Main executable
if(FLOAT_TETWILD_TOPLEVEL_PROJECT)
	add_executable(${PROJECT_NAME}_bin src/main.cpp)
	target_compile_features(${PROJECT_NAME}_bin PUBLIC ${CXX14_FEATURES})

	target_link_libraries(${PROJECT_NAME}_bin
		PUBLIC
			${PROJECT_NAME}
			warnings::all
	)

	if(FLOAT_TETWILD_WITH_SANITIZERS)
		add_sanitizers(${PROJECT_NAME}_bin)
	endif()

	if(NOT (${CMAKE_VERSION} VERSION_LESS "3.6.0"))
		set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${PROJECT_NAME}_bin)
	endif()
endif()

################################################################################
# Required libraries
################################################################################



# TBB
if(FLOAT_TETWILD_ENABLE_TBB)
	float_tetwild_download_tbb()
	set(TBB_BUILD_STATIC ON CACHE BOOL " " FORCE)
	set(TBB_BUILD_SHARED OFF CACHE BOOL " " FORCE)
	set(TBB_BUILD_TBBMALLOC OFF CACHE BOOL " " FORCE)
	set(TBB_BUILD_TBBMALLOC_PROXY OFF CACHE BOOL " " FORCE)
	set(TBB_BUILD_TESTS OFF CACHE BOOL " " FORCE)

	add_subdirectory(${THIRD_PARTY_DIR}/tbb tbb)
	set_property(TARGET tbb_static tbb_def_files PROPERTY FOLDER "dependencies")

	target_compile_definitions(tbb_static PUBLIC -DUSE_TBB)
	target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${THIRD_PARTY_DIR}/tbb/include)
	target_link_libraries(${PROJECT_NAME} PUBLIC tbb_static)
endif()


if(${FLOAT_TETWILD_USE_FLOAT})
	target_compile_definitions(${PROJECT_NAME} PUBLIC -DFLOAT_TETWILD_USE_FLOAT)
endif()

# libigl
float_tetwild_download_libigl()
find_package(LIBIGL REQUIRED)
target_link_libraries(${PROJECT_NAME} PUBLIC igl::core)

# Geogram
float_tetwild_download_geogram()
include(geogram)
target_link_libraries(${PROJECT_NAME} PUBLIC geogram)



#spdlog
target_link_libraries(${PROJECT_NAME} PUBLIC spdlog::spdlog)
target_link_libraries(${PROJECT_NAME} PUBLIC CLI11::CLI11)

# Json
# float_tetwild_download_json()
# add_library(json INTERFACE)
# target_include_directories(json SYSTEM INTERFACE ${THIRD_PARTY_DIR}/json/include)
# target_link_libraries(FloatTetwild PUBLIC json)


# Catch2
float_tetwild_download_catch2()
add_library(catch INTERFACE)
target_include_directories(catch SYSTEM INTERFACE ${THIRD_PARTY_DIR}/Catch2/single_include/catch2)


#aabbcc
# float_tetwild_download_aabbcc()
# add_library(aabbcc ${THIRD_PARTY_DIR}/aabbcc/src/AABB.h ${THIRD_PARTY_DIR}/aabbcc/src/AABB.cc)
# target_include_directories(aabbcc PUBLIC ${THIRD_PARTY_DIR}/aabbcc/src)
# target_link_libraries(${PROJECT_NAME} PUBLIC aabbcc)

if(LIBIGL_WITH_TETGEN)
	target_compile_definitions(${PROJECT_NAME} PUBLIC -DLIBIGL_WITH_TETGEN)
	target_link_libraries(${PROJECT_NAME} PUBLIC igl::tetgen)
endif()

################################################################################
# Optional libraries
################################################################################


################################################################################
# Clutter management
################################################################################

################################################################################
# Subdirectories
################################################################################

# Other sources
add_subdirectory(src)

# Compile extras only if this is a top-level project
if(FLOAT_TETWILD_TOPLEVEL_PROJECT)
	# Auxiliary binaries

	# Unit tests
	include(CTest)
	enable_testing()
	add_subdirectory(tests)
endif()
