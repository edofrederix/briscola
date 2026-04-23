# Source files

# Define executable
add_executable(briscolaColocatedTwoPhase ${LIB_SOURCES})

# Define sources
target_sources(briscolaColocatedTwoPhase PRIVATE briscolaColocatedTwoPhase.C)

# Dependencies
target_link_libraries(briscolaColocatedTwoPhase PRIVATE briscolaTwoPhase)

# OpenFOAM configuration
openfoam_configure_target(briscolaColocatedTwoPhase)

# Install
install(TARGETS briscolaColocatedTwoPhase RUNTIME DESTINATION ${FOAM_USER_APPBIN})
