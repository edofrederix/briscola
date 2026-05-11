# Try PETSC_HOME first
if(DEFINED ENV{PETSC_HOME})
    set(_PETSC_HINT "$ENV{PETSC_HOME}")
endif()

# pkg-config fallback
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
    pkg_check_modules(PC_PETSC QUIET petsc)
endif()

# Include dir
find_path(PETSC_INCLUDE_DIR
    NAMES petsc.h
    HINTS
        ${_PETSC_HINT}/include
        ${PC_PETSC_INCLUDE_DIRS}
)

# Library
find_library(PETSC_LIBRARY
    NAMES petsc
    HINTS
        ${_PETSC_HINT}/lib
        ${PC_PETSC_LIBRARY_DIRS}
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PETSc
    REQUIRED_VARS PETSC_INCLUDE_DIR PETSC_LIBRARY
)

# Imported target
if(PETSC_FOUND AND NOT TARGET PETSc)
    add_library(PETSc UNKNOWN IMPORTED)

    set_target_properties(PETSc PROPERTIES
        IMPORTED_LOCATION "${PETSC_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${PETSC_INCLUDE_DIR}"
        INTERFACE_COMPILE_DEFINITIONS "PETSC"
    )
endif()
