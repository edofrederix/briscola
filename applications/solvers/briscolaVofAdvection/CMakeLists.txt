# Define executable
add_executable(briscolaVofAdvection briscolaVofAdvection.C)

# Dependencies
target_link_libraries(briscolaVofAdvection PRIVATE briscolaTwoPhase)

# OpenFOAM configuration
openfoam_configure_target(briscolaVofAdvection)

# Install
install(TARGETS briscolaVofAdvection RUNTIME DESTINATION bin)
