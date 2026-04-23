# Try EIGEN_HOME first
if(DEFINED ENV{EIGEN_HOME})
    set(_EIGEN_HINT "$ENV{EIGEN_HOME}")
endif()

# pkg-config fallback
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
    pkg_check_modules(PC_EIGEN QUIET eigen3)
endif()

# Include dir
find_path(EIGEN_INCLUDE_DIR
    NAMES Eigen/Dense
    HINTS
        ${_EIGEN_HINT}/include/eigen3
        ${PC_EIGEN_INCLUDE_DIRS}
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Eigen3 REQUIRED_VARS EIGEN_INCLUDE_DIR)

# Imported target
if(EIGEN3_FOUND AND NOT TARGET Eigen3)
    add_library(Eigen3 INTERFACE IMPORTED)

    set_target_properties(Eigen3 PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${EIGEN_INCLUDE_DIR}"
        INTERFACE_COMPILE_DEFINITIONS "EIGEN3"
    )
endif()
