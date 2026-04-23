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

# Core required OF libs
find_library(OPENFOAM_LIB OpenFOAM HINTS $ENV{FOAM_LIBBIN})

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
