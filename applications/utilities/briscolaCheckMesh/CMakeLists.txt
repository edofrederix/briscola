# Define executable
add_executable(briscolaCheckMesh briscolaCheckMesh.C)

# Dependencies
target_link_libraries(briscolaCheckMesh PRIVATE briscolaFiniteVolume)

# OpenFOAM configuration
openfoam_configure_target(briscolaCheckMesh)

# Install
install(TARGETS briscolaCheckMesh RUNTIME DESTINATION bin)
