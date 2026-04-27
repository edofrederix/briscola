# Define executable
add_executable(parallelBoundaryCondition Test-parallelBoundaryCondition.C)

# Dependencies
target_link_libraries(parallelBoundaryCondition PRIVATE briscolaFiniteVolume)

# OpenFOAM configuration
openfoam_configure_target(parallelBoundaryCondition)

# Add to tests target
add_dependencies(tests parallelBoundaryCondition)

# Add test
add_test(NAME parallelBoundaryCondition
    COMMAND
        ${MPIEXEC_EXECUTABLE}
        ${MPIEXEC_NUMPROC_FLAG} 27 --oversubscribe
        $<TARGET_FILE:parallelBoundaryCondition> -parallel
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

set_tests_properties(parallelBoundaryCondition PROPERTIES PROCESSORS 27)
