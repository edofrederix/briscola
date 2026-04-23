# This module requires the following environment variables to be set:
#
#   WM_PROJECT
#   WM_PROJECT_DIR
#   WM_PROJECT_VERSION
#   WM_COMPILER
#   WM_COMPILE_OPTION
#   WM_ARCH
#   WM_ARCH_OPTION
#   WM_LABEL_SIZE
#   WM_LABEL_OPTION
#   WM_PRECISION_OPTION
#
# These variables are set by the OpenFOAM bashrc or cshrc scripts.

# Locate OpenFOAM installation via environment variables
if(NOT DEFINED ENV{WM_PROJECT_DIR})
    message(FATAL_ERROR
        "OpenFOAM environment not set. "
        "Load OpenFOAM before running CMake.")
endif()

# Required include paths from OpenFOAM
set(OPENFOAM_INCLUDE_DIR
    $ENV{WM_PROJECT_DIR}/src/OpenFOAM/lnInclude
    $ENV{WM_PROJECT_DIR}/src/Pstream/mpi/lnInclude
    $ENV{WM_PROJECT_DIR}/src/OSspecific/POSIX/lnInclude
)

# Reconstruct FOAM environment variables

set(WM_OPTIONS
    $ENV{WM_ARCH}$ENV{WM_COMPILER}$ENV{WM_PRECISION_OPTION}$ENV{WM_LABEL_OPTION}$ENV{WM_COMPILE_OPTION}
)

set(FOAM_LIBBIN
    $ENV{WM_PROJECT_DIR}/platforms/${WM_OPTIONS}/lib
)

set(FOAM_APPBIN
    $ENV{WM_PROJECT_DIR}/platforms/{$WM_OPTIONS}/bin
)

set(FOAM_USER_LIBBIN
    $ENV{HOME}/$ENV{WM_PROJECT}/$ENV{USER}-$ENV{WM_PROJECT_VERSION}/platforms/${WM_OPTIONS}/lib
)

set(FOAM_USER_APPBIN
    $ENV{HOME}/$ENV{WM_PROJECT}/$ENV{USER}-$ENV{WM_PROJECT_VERSION}/platforms/${WM_OPTIONS}/bin
)

# Core required OF libs
find_library(OPENFOAM_LIB OpenFOAM HINTS ${FOAM_LIBBIN})

# Compile definitions
set(OPENFOAM_COMPILE_DEFINITIONS
    WM_ARCH_OPTION=$ENV{WM_ARCH_OPTION}
    WM_$ENV{WM_PRECISION_OPTION}
    WM_LABEL_SIZE=$ENV{WM_LABEL_SIZE}
    $ENV{WM_ARCH}
    NoRepository
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OpenFOAM DEFAULT_MSG OPENFOAM_LIB)

# Target
if(OPENFOAM_FOUND AND NOT TARGET OpenFOAM)
    add_library(OpenFOAM UNKNOWN IMPORTED)

    set_target_properties(OpenFOAM PROPERTIES
        IMPORTED_LOCATION "${OPENFOAM_LIB}"
        INTERFACE_INCLUDE_DIRECTORIES "${OPENFOAM_INCLUDE_DIR}"
    )
endif()
