# Define executable
add_executable(metrics Test-metrics.C)

# Dependencies
target_link_libraries(metrics PRIVATE briscolaFiniteVolume)

# Add to tests target
add_dependencies(tests metrics)

# Add test
add_test(NAME metrics
    COMMAND
        ${MPIEXEC_EXECUTABLE}
        ${MPIEXEC_NUMPROC_FLAG} 4 --oversubscribe
        $<TARGET_FILE:metrics> -parallel
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

set_tests_properties(metrics PROPERTIES PROCESSORS 4)
