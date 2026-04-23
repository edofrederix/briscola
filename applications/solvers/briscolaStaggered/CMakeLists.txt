# Define executable
add_executable(briscolaStaggered briscolaStaggered.C)

# Dependencies
target_link_libraries(briscolaStaggered PRIVATE briscolaFiniteVolume)

# OpenFOAM configuration
openfoam_configure_target(briscolaStaggered)

# Install
install(TARGETS briscolaStaggered RUNTIME DESTINATION $ENV{FOAM_USER_APPBIN})
