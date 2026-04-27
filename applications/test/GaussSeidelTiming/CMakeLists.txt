# Define executable
add_executable(GaussSeidelTiming Test-GaussSeidelTiming.C)

# Dependencies
target_link_libraries(GaussSeidelTiming PRIVATE briscolaFiniteVolume)

# OpenFOAM configuration
openfoam_configure_target(GaussSeidelTiming)

# Add to tests target
add_dependencies(tests GaussSeidelTiming)

# Add test
add_test(NAME GaussSeidelTiming
    COMMAND
        ${MPIEXEC_EXECUTABLE}
        ${MPIEXEC_NUMPROC_FLAG} 8 --oversubscribe
        $<TARGET_FILE:GaussSeidelTiming> -parallel
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

set_tests_properties(GaussSeidelTiming PROPERTIES PROCESSORS 8)
