# Define executable
add_executable(briscolaColocatedTwoPhase briscolaColocatedTwoPhase.C)

# Dependencies
target_link_libraries(briscolaColocatedTwoPhase PRIVATE briscolaTwoPhase)

# OpenFOAM configuration
openfoam_configure_target(briscolaColocatedTwoPhase)

# Install
install(TARGETS briscolaColocatedTwoPhase RUNTIME DESTINATION bin)
