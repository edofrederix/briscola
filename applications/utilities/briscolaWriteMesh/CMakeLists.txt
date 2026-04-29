# Define executable
add_executable(briscolaWriteMesh briscolaWriteMesh.C)

# Dependencies
target_link_libraries(briscolaWriteMesh PRIVATE briscolaFiniteVolume)

# OpenFOAM configuration
openfoam_configure_target(briscolaWriteMesh)

# Install
install(TARGETS briscolaWriteMesh RUNTIME DESTINATION bin)

# Add to utilities target
add_dependencies(utilities briscolaWriteMesh)
