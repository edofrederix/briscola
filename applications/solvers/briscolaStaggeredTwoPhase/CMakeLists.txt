# Source files

# Define executable
add_executable(briscolaStaggeredTwoPhase ${LIB_SOURCES})

# Define sources
target_sources(briscolaStaggeredTwoPhase PRIVATE briscolaStaggeredTwoPhase.C)

# Dependencies
target_link_libraries(briscolaStaggeredTwoPhase PRIVATE briscolaTwoPhase)

# OpenFOAM configuration
openfoam_configure_target(briscolaStaggeredTwoPhase)

# Install
install(TARGETS briscolaStaggeredTwoPhase RUNTIME DESTINATION ${FOAM_USER_APPBIN})
