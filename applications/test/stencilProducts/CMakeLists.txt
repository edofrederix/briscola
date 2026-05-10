# Define executable
add_executable(stencilProducts Test-stencilProducts.C)

# Dependencies
target_link_libraries(stencilProducts PRIVATE briscolaFiniteVolume)

# Add to tests target
add_dependencies(tests stencilProducts)

# Add test
add_test(NAME stencilProducts
    COMMAND
        ${MPIEXEC_EXECUTABLE}
        ${MPIEXEC_NUMPROC_FLAG} 4 --oversubscribe
        $<TARGET_FILE:stencilProducts> -parallel
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

set_tests_properties(stencilProducts PROPERTIES PROCESSORS 4)
