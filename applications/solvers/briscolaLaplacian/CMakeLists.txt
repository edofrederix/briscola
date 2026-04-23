# Source files

# Define executable
add_executable(briscolaLaplacian ${LIB_SOURCES})

# Define sources
target_sources(briscolaLaplacian PRIVATE briscolaLaplacian.C)

# Dependencies
target_link_libraries(briscolaLaplacian PRIVATE briscolaFiniteVolume)

# OpenFOAM configuration
openfoam_configure_target(briscolaLaplacian)

# Install
install(TARGETS briscolaLaplacian RUNTIME DESTINATION ${FOAM_USER_APPBIN})
