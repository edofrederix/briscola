# Define executable
add_executable(briscolaStaggeredTwoPhase briscolaStaggeredTwoPhase.C)

# Dependencies
target_link_libraries(briscolaStaggeredTwoPhase PRIVATE briscolaTwoPhase)

# OpenFOAM configuration
openfoam_configure_target(briscolaStaggeredTwoPhase)

# Install
install(TARGETS briscolaStaggeredTwoPhase RUNTIME DESTINATION $ENV{FOAM_USER_APPBIN})
