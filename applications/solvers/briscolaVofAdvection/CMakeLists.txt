# Source files

# Define executable
add_executable(briscolaVofAdvection ${LIB_SOURCES})

# Define sources
target_sources(briscolaVofAdvection PRIVATE briscolaVofAdvection.C)

# Dependencies
target_link_libraries(briscolaVofAdvection PRIVATE briscolaTwoPhase)

# OpenFOAM configuration
openfoam_configure_target(briscolaVofAdvection)

# Install
install(TARGETS briscolaVofAdvection RUNTIME DESTINATION ${FOAM_USER_APPBIN})
