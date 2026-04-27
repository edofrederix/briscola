# Define executable
add_executable(externalSolverTiming Test-externalSolverTiming.C)

# Dependencies
target_link_libraries(externalSolverTiming PRIVATE briscolaFiniteVolume)

# OpenFOAM configuration
openfoam_configure_target(externalSolverTiming)

# Add to tests target
add_dependencies(tests externalSolverTiming)

# Add test
add_test(NAME externalSolverTiming
    COMMAND
        ${MPIEXEC_EXECUTABLE}
        ${MPIEXEC_NUMPROC_FLAG} 8 --oversubscribe
        $<TARGET_FILE:externalSolverTiming> -parallel
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

set_tests_properties(externalSolverTiming PROPERTIES PROCESSORS 8)
