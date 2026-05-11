# Define executable
add_executable(matrixMultiplicationTiming Test-matrixMultiplicationTiming.C)

# Dependencies
target_link_libraries(matrixMultiplicationTiming PRIVATE briscolaFiniteVolume)

# Add to tests target
add_dependencies(tests matrixMultiplicationTiming)

# Add test
add_test(NAME matrixMultiplicationTiming
    COMMAND matrixMultiplicationTiming
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
