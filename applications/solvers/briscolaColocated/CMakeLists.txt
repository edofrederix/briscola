# Source files

# Define executable
add_executable(briscolaColocated ${LIB_SOURCES})

# Define sources
target_sources(briscolaColocated PRIVATE briscolaColocated.C)

# Dependencies
target_link_libraries(briscolaColocated PRIVATE briscolaFiniteVolume)

# OpenFOAM configuration
openfoam_configure_target(briscolaColocated)

# Install
install(TARGETS briscolaColocated RUNTIME DESTINATION ${FOAM_USER_APPBIN})
