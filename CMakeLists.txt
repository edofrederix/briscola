cmake_minimum_required(VERSION 3.15)
project(Briscola CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Build types

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type")
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
  "Debug" "Release" "Profile"
)

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

# Add sub-directories

add_subdirectory(src)
add_subdirectory(applications)
