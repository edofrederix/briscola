# Define executable
add_executable(briscolaWriteRestriction briscolaWriteRestriction.C)

# Dependencies
target_link_libraries(briscolaWriteRestriction PRIVATE briscolaFiniteVolume)

# OpenFOAM configuration
openfoam_configure_target(briscolaWriteRestriction)

# Install
install(TARGETS briscolaWriteRestriction RUNTIME DESTINATION bin)
