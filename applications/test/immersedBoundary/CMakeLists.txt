# Define executable
add_executable(immersedBoundary Test-immersedBoundary.C)

# Dependencies
target_link_libraries(immersedBoundary PRIVATE briscolaFiniteVolume)

# Add to tests target
add_dependencies(tests immersedBoundary)

# Add test
add_test(NAME immersedBoundary
    COMMAND $<TARGET_FILE:immersedBoundary>
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
