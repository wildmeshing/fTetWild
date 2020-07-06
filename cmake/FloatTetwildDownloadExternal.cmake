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
        URL          https://github.com/catchorg/Catch2/archive/v2.4.2.tar.gz
        URL_MD5      26927b878b1f42633f15a9ef1c4bd8e7
    )
endfunction()

## CLI11
function(float_tetwild_download_cli11)
    float_tetwild_download_project(cli11
            URL     https://github.com/CLIUtils/CLI11/archive/v1.8.0.tar.gz
            URL_MD5 5e5470abcb76422360409297bfc446ac
    )
endfunction()

## tbb
function(float_tetwild_download_tbb)
    float_tetwild_download_project(tbb
        GIT_REPOSITORY https://github.com/wjakob/tbb.git
        GIT_TAG        ddbe45cd3ad89df9a84cd77013d5898fc48b8e89
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
        URL      https://github.com/fmtlib/fmt/archive/5.3.0.tar.gz
        URL_MD5  1015bf3ff2a140dfe03de50ee2469401
    )
endfunction()

## spdlog
function(float_tetwild_download_spdlog)
    float_tetwild_download_project(spdlog
        URL         https://github.com/gabime/spdlog/archive/v1.3.1.tar.gz
        URL_MD5     3c17dd6983de2a66eca8b5a0b213d29f
    )
endfunction()

## Geogram LGPL
function(float_tetwild_download_geogram)
    float_tetwild_download_project(geogram
        URL         https://github.com/alicevision/geogram/archive/v1.6.8.tar.gz
        URL_MD5     9be60608f984de909fcf06ba421f970a
    )
endfunction()

## aabbcc
function(float_tetwild_download_aabbcc)
    float_tetwild_download_project(aabbcc
            GIT_REPOSITORY https://github.com/lohedges/aabbcc.git
            GIT_TAG        0c85e61362d384d70c71946826bfed0fb24a74ba
            )
endfunction()

## winding number
function(float_tetwild_download_windingnumber)
    float_tetwild_download_project(windingnumber
            GIT_REPOSITORY https://github.com/alecjacobson/WindingNumber.git
            GIT_TAG        1e6081e52905575d8e98fb8b7c0921274a18752f
            )
endfunction()
