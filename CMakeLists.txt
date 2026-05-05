cmake_minimum_required(VERSION 3.22)
project(Briscola CXX)

# Include main Briscola CMake file

include(cmake/Briscola.cmake)

# Require OpenFOAM

find_package(OpenFOAM REQUIRED)
find_package(Pstream REQUIRED)

# Other dependencies

find_package(MPI REQUIRED)
find_package(FFTW3 REQUIRED)
find_package(PETSc)
find_package(Eigen3)

# Enable testing

enable_testing()

# Add sub-directories

add_subdirectory(src)
add_subdirectory(applications)
