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
        GIT_TAG        45cfc79fede992ea3923ded9de3c21d1c4faced1
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
        URL          https://github.com/catchorg/Catch2/archive/refs/tags/v3.1.1.tar.gz
        URL_MD5      1f3e0d8c3297252f77d643ff06d058cb
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
#        GIT_REPOSITORY   https://github.com/polyfem/geogram.git
#        GIT_TAG          e6b9612f1146370e40deaa341b4dd7ef90502102
            GIT_REPOSITORY https://github.com/Yixin-Hu/geogram
            GIT_TAG        b613750341a6cdd31ae8df80ecfc26ac7ca1a6ad
    )
endfunction()

## aabbcc
function(float_tetwild_download_aabbcc)
    float_tetwild_download_project(aabbcc
            GIT_REPOSITORY https://github.com/lohedges/aabbcc.git
            GIT_TAG        0c85e61362d384d70c71946826bfed0fb24a74ba
            )
endfunction()

### winding number
#function(float_tetwild_download_windingnumber)
#    float_tetwild_download_project(windingnumber
#            GIT_REPOSITORY https://github.com/alecjacobson/WindingNumber.git
#            GIT_TAG        bde8780ec848fa71c1294a0af50347e968b19493
#            )
#endfunction()


## exact envelope
function(float_tetwild_download_exact_envelope)
    float_tetwild_download_project(exact_envelope
            GIT_REPOSITORY https://github.com/wangbolun300/fast-envelope
            GIT_TAG        520ee04b6c69a802db31d1fd3a3e6e382d10ef98
            )
endfunction()