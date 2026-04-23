# Define executable
add_executable(writeProlongation writeProlongation.C)

# Dependencies
target_link_libraries(writeProlongation PRIVATE briscolaFiniteVolume)

# OpenFOAM configuration
openfoam_configure_target(writeProlongation)

# Install
install(TARGETS writeProlongation RUNTIME DESTINATION $ENV{FOAM_USER_APPBIN})
