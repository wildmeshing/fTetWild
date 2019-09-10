################################################################################

include(DownloadProject)

# With CMake 3.8 and above, we can hide warnings about git being in a
# detached head by passing an extra GIT_CONFIG option
if(NOT (${CMAKE_VERSION} VERSION_LESS "3.8.0"))
    set(FLOAT_TETWILD_EXTRA_OPTIONS "GIT_CONFIG advice.detachedHead=false")
else()
    set(FLOAT_TETWILD_EXTRA_OPTIONS "")
endif()

# Shortcut function
function(float_tetwild_download_project name)
    download_project(
        PROJ         ${name}
        SOURCE_DIR   ${FLOAT_TETWILD_EXTERNAL}/${name}
        DOWNLOAD_DIR ${FLOAT_TETWILD_EXTERNAL}/.cache/${name}
        QUIET
        ${FLOAT_TETWILD_EXTRA_OPTIONS}
        ${ARGN}
    )
endfunction()

################################################################################

## libigl
function(float_tetwild_download_libigl)
    float_tetwild_download_project(libigl
        GIT_REPOSITORY https://github.com/libigl/libigl.git
        GIT_TAG        b0d7740e0b7e887a7e93601c4c557ecf762b389b
    )
endfunction()

## Json
function(float_tetwild_download_json)
    float_tetwild_download_project(json
        GIT_REPOSITORY https://github.com/jdumas/json
        GIT_TAG        0901d33bf6e7dfe6f70fd9d142c8f5c6695c6c5b
    )
endfunction()

## Catch2
function(float_tetwild_download_catch2)
    float_tetwild_download_project(Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v2.4.2
    )
endfunction()

## CLI11
function(float_tetwild_download_cli11)
    float_tetwild_download_project(cli11
        URL     https://github.com/CLIUtils/CLI11/archive/v1.6.1.tar.gz
        URL_MD5 48ef97262adb0b47a2f0a7edbda6e2aa
    )
endfunction()

## tbb
function(float_tetwild_download_tbb)
    float_tetwild_download_project(tbb
        GIT_REPOSITORY https://github.com/wjakob/tbb.git
        GIT_TAG        08b4341a1893a72656467e96137f1f99d0112547
    )
endfunction()

## Sanitizers
function(float_tetwild_download_sanitizers)
    float_tetwild_download_project(sanitizers-cmake
        GIT_REPOSITORY https://github.com/arsenm/sanitizers-cmake.git
        GIT_TAG        6947cff3a9c9305eb9c16135dd81da3feb4bf87f
    )
endfunction()

## fmt
function(float_tetwild_download_fmt)
    float_tetwild_download_project(fmt
        GIT_REPOSITORY https://github.com/fmtlib/fmt
        GIT_TAG 5.3.0
    )
endfunction()

## spdlog
function(float_tetwild_download_spdlog)
    float_tetwild_download_project(spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG         v1.3.1
    )
endfunction()

## Geogram LGPL
function(float_tetwild_download_geogram)
    float_tetwild_download_project(geogram
        GIT_REPOSITORY https://github.com/alicevision/geogram.git
        GIT_TAG        v1.6.8
    )
endfunction()

## aabbcc
function(float_tetwild_download_aabbcc)
    float_tetwild_download_project(aabbcc
            GIT_REPOSITORY https://github.com/lohedges/aabbcc.git
            GIT_TAG        0c85e61362d384d70c71946826bfed0fb24a74ba
            )
endfunction()
