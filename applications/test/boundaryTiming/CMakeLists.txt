# Define executable
add_executable(boundaryTiming Test-boundaryTiming.C)

# Dependencies
target_link_libraries(boundaryTiming PRIVATE briscolaFiniteVolume)

# OpenFOAM configuration
openfoam_configure_target(boundaryTiming)

# Add to tests target
add_dependencies(tests boundaryTiming)

# Add test
add_test(NAME boundaryTiming
    COMMAND
        ${MPIEXEC_EXECUTABLE}
        ${MPIEXEC_NUMPROC_FLAG} 64 --oversubscribe
        $<TARGET_FILE:boundaryTiming> -parallel
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

set_tests_properties(boundaryTiming PROPERTIES PROCESSORS 64)
