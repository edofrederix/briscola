cmake_minimum_required(VERSION 3.22)
project(Briscola CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Build types

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type")
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
  "Debug" "Release" "Profile"
)

# Store built libraries and applications in lib and bin in the build directory

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Helper modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Load OpenFOAM from environment

find_package(OpenFOAM REQUIRED)
find_package(Pstream REQUIRED)
include(OpenFOAMConfig)

# Dependencies

find_package(MPI REQUIRED)
find_package(FFTW3 REQUIRED)
find_package(PETSc)
find_package(Eigen3)

# Find Python3 and check for numpy availability

find_package(Python3 COMPONENTS Interpreter)

if(Python3_FOUND)
    execute_process(
        COMMAND ${Python3_EXECUTABLE} -c "import numpy"
        RESULT_VARIABLE NUMPY_RESULT
        OUTPUT_QUIET
        ERROR_QUIET
    )
endif()

# Default prefix points to OpenFOAM's user install location

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    cmake_path(SET install_path NORMALIZE "$ENV{FOAM_USER_LIBBIN}/../")
    set(CMAKE_INSTALL_PREFIX ${install_path} CACHE PATH "Install prefix" FORCE)
    message(STATUS "Using OpenFOAM's user install location as prefix:
        ${install_path}")
endif()

# Add sub-directories

add_subdirectory(src)
add_subdirectory(applications)
