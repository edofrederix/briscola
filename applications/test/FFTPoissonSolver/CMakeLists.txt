# Define executable
add_executable(FFTPoissonSolver Test-FFTPoissonSolver.C)

# Dependencies
target_link_libraries(FFTPoissonSolver PRIVATE briscolaFiniteVolume FFTW3)

# Add to tests target
add_dependencies(tests FFTPoissonSolver)

# Add test
add_test(NAME FFTPoissonSolver
    COMMAND
        ${MPIEXEC_EXECUTABLE}
        ${MPIEXEC_NUMPROC_FLAG} 8 --oversubscribe
        $<TARGET_FILE:FFTPoissonSolver> -parallel
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

set_tests_properties(FFTPoissonSolver PROPERTIES PROCESSORS 8)
