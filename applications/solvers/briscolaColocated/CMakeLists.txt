# Define executable
add_executable(briscolaColocated briscolaColocated.C)

# Dependencies
target_link_libraries(briscolaColocated PRIVATE briscolaFiniteVolume)

# OpenFOAM configuration
openfoam_configure_target(briscolaColocated)

# Install
install(TARGETS briscolaColocated RUNTIME DESTINATION $ENV{FOAM_USER_APPBIN})
