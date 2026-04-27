# Define executable
add_executable(matrixMultiplicationTiming Test-matrixMultiplicationTiming.C)

# Dependencies
target_link_libraries(matrixMultiplicationTiming PRIVATE briscolaFiniteVolume)

# OpenFOAM configuration
openfoam_configure_target(matrixMultiplicationTiming)

# Add to tests target
add_dependencies(tests matrixMultiplicationTiming)

# Add test
add_test(NAME matrixMultiplicationTiming
    COMMAND matrixMultiplicationTiming
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
