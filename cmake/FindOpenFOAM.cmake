# Locate OpenFOAM installation via environment variables
if(NOT DEFINED ENV{WM_PROJECT_DIR})
    message(FATAL_ERROR
        "OpenFOAM environment not set. "
        "Load OpenFOAM before running CMake.")
endif()

# Required include paths from OpenFOAM
set(OPENFOAM_INCLUDE_DIRS
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
    OMPI_SKIP_MPICXX
)

# Compile options
set(OPENFOAM_COMPILE_OPTIONS

    # Default OpenFOAM warning flags
    -Wall
    -Wextra
    -Wold-style-cast
    -Wnon-virtual-dtor
    -Wno-unused-parameter
    -Wno-invalid-offsetof
    -Wno-attributes
    -ftemplate-depth-100

    # Custom warning flags
    -Wshadow

    # Build type flags
    $<$<CONFIG:Release>:-O3>
    $<$<CONFIG:Debug>:-g -O0 -DFULLDEBUG>
    $<$<CONFIG:Profile>:-g -O2>
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OpenFOAM DEFAULT_MSG OPENFOAM_LIB)

# Set OpenFOAM target
if(OPENFOAM_FOUND AND NOT TARGET OpenFOAM)
    add_library(OpenFOAM UNKNOWN IMPORTED)

    set_target_properties(OpenFOAM PROPERTIES
        IMPORTED_LOCATION "${OPENFOAM_LIB}"
        INTERFACE_INCLUDE_DIRECTORIES "${OPENFOAM_INCLUDE_DIRS}"
        INTERFACE_COMPILE_DEFINITIONS "${OPENFOAM_COMPILE_DEFINITIONS}"
        INTERFACE_COMPILE_OPTIONS "${OPENFOAM_COMPILE_OPTIONS}"
    )
endif()
