# Source files

# Define executable
add_executable(writeProlongation ${LIB_SOURCES})

# Define sources
target_sources(writeProlongation PRIVATE writeProlongation.C)

# Dependencies
target_link_libraries(writeProlongation PRIVATE briscolaFiniteVolume)

# OpenFOAM configuration
openfoam_configure_target(writeProlongation)

# Install
install(TARGETS writeProlongation RUNTIME DESTINATION ${FOAM_USER_APPBIN})
