# Define executable
add_executable(briscolaLaplacian briscolaLaplacian.C)

# Dependencies
target_link_libraries(briscolaLaplacian PRIVATE briscolaFiniteVolume)

# OpenFOAM configuration
openfoam_configure_target(briscolaLaplacian)

# Install
install(TARGETS briscolaLaplacian RUNTIME DESTINATION bin)

# Add to solvers target
add_dependencies(solvers briscolaLaplacian)
