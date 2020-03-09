# Try to find the GMP librairies
# GMP_FOUND - system has GMP lib
# GMP_INCLUDE_DIRS - the GMP include directory
# GMP_LIBRARIES - Libraries needed to use GMP

if (GMP_INCLUDE_DIRS AND GMP_LIBRARIES)
        # Already in cache, be silent
        set(GMP_FIND_QUIETLY TRUE)
endif (GMP_INCLUDE_DIRS AND GMP_LIBRARIES)



#if(WIN32)
#	if(CYGWIN)
#		triwild_download_gmp_cygwin()
	#elseif(MINGW)
#		triwild_download_gmp_mingw()
	#else()
#		triwild_download_gmp_vc()
	#endif()

	#SET(GMP_WINDOWS_PATH ${THIRD_PARTY_DIR}/gmp)
#endif()


find_path(GMP_INCLUDE_DIRS NAMES gmp.h PATHS $ENV{GMP_INC} ${GMP_WINDOWS_PATH})
find_library(GMP_LIBRARIES NAMES gmp libgmp PATHS $ENV{GMP_LIB} ${GMP_WINDOWS_PATH})
find_library(GMPXX_LIBRARIES NAMES gmpxx libgmpxx PATHS $ENV{GMP_LIB} ${GMP_WINDOWS_PATH})
#MESSAGE(STATUS "GMP libs: " ${GMP_LIBRARIES} " " ${GMPXX_LIBRARIES} )

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(GMP DEFAULT_MSG GMP_INCLUDE_DIRS GMP_LIBRARIES)

mark_as_advanced(GMP_INCLUDE_DIRS GMP_LIBRARIES)
MESSAGE(STATUS "GMP libs: " ${GMP_LIBRARIES} " " ${GMP_INCLUDE_DIRS} )