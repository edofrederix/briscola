# Define executable
add_executable(writeRestriction writeRestriction.C)

# Dependencies
target_link_libraries(writeRestriction PRIVATE briscolaFiniteVolume)

# OpenFOAM configuration
openfoam_configure_target(writeRestriction)

# Install
install(TARGETS writeRestriction RUNTIME DESTINATION $ENV{FOAM_USER_APPBIN})
