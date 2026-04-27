# Define executable
add_executable(prolongation Test-prolongation.C)

# Dependencies
target_link_libraries(prolongation PRIVATE briscolaFiniteVolume)

# OpenFOAM configuration
openfoam_configure_target(prolongation)

# Add to tests target
add_dependencies(tests prolongation)

# Add test
add_test(NAME prolongation
    COMMAND
        ${MPIEXEC_EXECUTABLE}
        ${MPIEXEC_NUMPROC_FLAG} 4 --oversubscribe
        $<TARGET_FILE:prolongation> -parallel
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

set_tests_properties(prolongation PROPERTIES PROCESSORS 4)
