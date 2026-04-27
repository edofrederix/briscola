# Define executable
add_executable(restriction Test-restriction.C)

# Dependencies
target_link_libraries(restriction PRIVATE briscolaFiniteVolume)

# OpenFOAM configuration
openfoam_configure_target(restriction)

# Add to tests target
add_dependencies(tests restriction)

# Add test
add_test(NAME restriction
    COMMAND
        ${MPIEXEC_EXECUTABLE}
        ${MPIEXEC_NUMPROC_FLAG} 4 --oversubscribe
        $<TARGET_FILE:restriction> -parallel
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

set_tests_properties(restriction PROPERTIES PROCESSORS 4)
