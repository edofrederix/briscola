# Source files

# Define executable
add_executable(briscolaCheckMesh ${LIB_SOURCES})

# Define sources
target_sources(briscolaCheckMesh PRIVATE briscolaCheckMesh.C)

# Dependencies
target_link_libraries(briscolaCheckMesh PRIVATE briscolaFiniteVolume)

# OpenFOAM configuration
openfoam_configure_target(briscolaCheckMesh)

# Install
install(TARGETS briscolaCheckMesh RUNTIME DESTINATION ${FOAM_USER_APPBIN})
