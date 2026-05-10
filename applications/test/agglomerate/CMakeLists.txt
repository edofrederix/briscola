# Define executable
add_executable(agglomerate Test-agglomerate.C)

# Dependencies
target_link_libraries(agglomerate PRIVATE briscolaFiniteVolume)

# Add to tests target
add_dependencies(tests agglomerate)

# Add test
add_test(NAME agglomerate
    COMMAND
        ${MPIEXEC_EXECUTABLE}
        ${MPIEXEC_NUMPROC_FLAG} 192 --oversubscribe
        $<TARGET_FILE:agglomerate> -parallel
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

set_tests_properties(agglomerate PROPERTIES PROCESSORS 192)
