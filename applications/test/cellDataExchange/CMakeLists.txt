# Define executable
add_executable(cellDataExchange Test-cellDataExchange.C)

# Dependencies
target_link_libraries(cellDataExchange PRIVATE briscolaFiniteVolume)

# Add to tests target
add_dependencies(tests cellDataExchange)

# Add test
add_test(NAME cellDataExchange
    COMMAND
        ${MPIEXEC_EXECUTABLE}
        ${MPIEXEC_NUMPROC_FLAG} 4 --oversubscribe
        $<TARGET_FILE:cellDataExchange> -parallel
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

set_tests_properties(cellDataExchange PROPERTIES PROCESSORS 4)
