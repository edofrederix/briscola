# Define executable
add_executable(cellDataExchangeUnstructured Test-cellDataExchangeUnstructured.C)

# Dependencies
target_link_libraries(cellDataExchangeUnstructured PRIVATE briscolaFiniteVolume)

# Add to tests target
add_dependencies(tests cellDataExchangeUnstructured)

# Add test
add_test(NAME cellDataExchangeUnstructured
    COMMAND
        ${MPIEXEC_EXECUTABLE}
        ${MPIEXEC_NUMPROC_FLAG} 5 --oversubscribe
        $<TARGET_FILE:cellDataExchangeUnstructured> -parallel
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

set_tests_properties(cellDataExchangeUnstructured PROPERTIES PROCESSORS 5)
