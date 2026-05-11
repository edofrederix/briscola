# Define executable
add_executable(FFTPoissonSolver2D Test-FFTPoissonSolver2D.C)

# Dependencies
target_link_libraries(FFTPoissonSolver2D PRIVATE briscolaFiniteVolume FFTW3)

# Add to tests target
add_dependencies(tests FFTPoissonSolver2D)

# Add test
add_test(NAME FFTPoissonSolver2D
    COMMAND
        ${MPIEXEC_EXECUTABLE}
        ${MPIEXEC_NUMPROC_FLAG} 4 --oversubscribe
        $<TARGET_FILE:FFTPoissonSolver2D> -parallel
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

set_tests_properties(FFTPoissonSolver2D PROPERTIES PROCESSORS 4)
