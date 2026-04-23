# Source files

# Define executable
add_executable(briscolaWriteMesh ${LIB_SOURCES})

# Define sources
target_sources(briscolaWriteMesh PRIVATE briscolaWriteMesh.C)

# Dependencies
target_link_libraries(briscolaWriteMesh PRIVATE briscolaFiniteVolume)

# OpenFOAM configuration
openfoam_configure_target(briscolaWriteMesh)

# Install
install(TARGETS briscolaWriteMesh RUNTIME DESTINATION ${FOAM_USER_APPBIN})
