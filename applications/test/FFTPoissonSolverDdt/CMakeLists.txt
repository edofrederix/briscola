# Define executable
add_executable(FFTPoissonSolverDdt Test-FFTPoissonSolverDdt.C)

# Dependencies
target_link_libraries(FFTPoissonSolverDdt PRIVATE briscolaFiniteVolume FFTW3)

# OpenFOAM configuration
openfoam_configure_target(FFTPoissonSolverDdt)

# Add to tests target
add_dependencies(tests FFTPoissonSolverDdt)

# Add test
add_test(NAME FFTPoissonSolverDdt
    COMMAND $<TARGET_FILE:FFTPoissonSolverDdt>
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
