# Try FFTW_HOME first
if(DEFINED ENV{FFTW_HOME})
    set(_FFTW_HINT "$ENV{FFTW_HOME}")
endif()

# Try pkg-config if available
find_package(PkgConfig QUIET)

if(PkgConfig_FOUND)
    pkg_check_modules(PC_FFTW3 QUIET fftw3)
endif()

# Find include dir
find_path(FFTW3_INCLUDE_DIR
    NAMES fftw3.h
    HINTS
        ${_FFTW_HINT}/include
        ${PC_FFTW3_INCLUDE_DIRS}
)

# Find library
find_library(FFTW3_LIBRARY
    NAMES fftw3
    HINTS
        ${_FFTW_HINT}/lib
        ${PC_FFTW3_LIBRARY_DIRS}
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW3
    REQUIRED_VARS FFTW3_INCLUDE_DIR FFTW3_LIBRARY
)

# Target
if(FFTW3_FOUND AND NOT TARGET FFTW3)
    add_library(FFTW3 UNKNOWN IMPORTED)

    set_target_properties(FFTW3 PROPERTIES
        IMPORTED_LOCATION "${FFTW3_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIR}"
        INTERFACE_COMPILE_DEFINITIONS "FFTW3"
    )
endif()
