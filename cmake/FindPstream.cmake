# Locate OpenFOAM installation via environment variables
if(NOT DEFINED ENV{WM_PROJECT_DIR})
    message(FATAL_ERROR
        "OpenFOAM environment not set. "
        "Load OpenFOAM before running CMake.")
endif()

# Required include paths from Pstream
set(PSTREAM_INCLUDE_DIRS
    $ENV{WM_PROJECT_DIR}/src/Pstream/mpi/lnInclude
)

# Core required OF libs
find_library(PSTREAM_LIB Pstream HINTS
    "${FOAM_LIBBIN}/openmpi-system"
    "${FOAM_LIBBIN}/mpi-system"
    "${FOAM_LIBBIN}/mpich-gm"
    "${FOAM_LIBBIN}/mvapich2"
    "${FOAM_LIBBIN}/hpmpi"
    "${FOAM_LIBBIN}/mpi"
    "${FOAM_LIBBIN}/fjmpi"
    "${FOAM_LIBBIN}/qsmpi"
    "${FOAM_LIBBIN}/dummy"
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Pstream DEFAULT_MSG PSTREAM_LIB)

# Target
if(PSTREAM_FOUND AND NOT TARGET Pstream)
    add_library(Pstream UNKNOWN IMPORTED)

    set_target_properties(Pstream PROPERTIES
        IMPORTED_LOCATION "${PSTREAM_LIB}"
        INTERFACE_INCLUDE_DIRECTORIES "${PSTREAM_INCLUDE_DIRS}"
    )
endif()
