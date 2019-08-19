# Prepare dependencies
#
# For each third-party library, if the appropriate target doesn't exist yet,
# download it via external project, and add_subdirectory to build it alongside
# this project.

### Configuration
set(FLOAT_TETWILD_ROOT     "${CMAKE_CURRENT_LIST_DIR}/..")
set(FLOAT_TETWILD_EXTERNAL ${THIRD_PARTY_DIR})

# Download and update 3rdparty libraries
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})
list(REMOVE_DUPLICATES CMAKE_MODULE_PATH)
include(FloatTetwildDownloadExternal)

################################################################################
# Required libraries
################################################################################

# Sanitizers
if(FLOAT_TETWILD_WITH_SANITIZERS)
    float_tetwild_download_sanitizers()
    find_package(Sanitizers)
endif()

# CL11
if(NOT TARGET CLI11::CLI11)
    float_tetwild_download_cli11()
    add_subdirectory(${FLOAT_TETWILD_EXTERNAL}/cli11)
endif()

# spdlog
if(NOT TARGET spdlog::spdlog)
	float_tetwild_download_spdlog()
	add_subdirectory(${FLOAT_TETWILD_EXTERNAL}/spdlog)
endif()

