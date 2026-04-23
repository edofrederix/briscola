# Source files

# Define executable
add_executable(briscolaStaggered ${LIB_SOURCES})

# Define sources
target_sources(briscolaStaggered PRIVATE briscolaStaggered.C)

# Dependencies
target_link_libraries(briscolaStaggered PRIVATE briscolaFiniteVolume)

# OpenFOAM configuration
openfoam_configure_target(briscolaStaggered)

# Install
install(TARGETS briscolaStaggered RUNTIME DESTINATION ${FOAM_USER_APPBIN})
