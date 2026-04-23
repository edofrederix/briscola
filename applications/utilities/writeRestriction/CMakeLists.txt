# Source files

# Define executable
add_executable(writeRestriction ${LIB_SOURCES})

# Define sources
target_sources(writeRestriction PRIVATE writeRestriction.C)

# Dependencies
target_link_libraries(writeRestriction PRIVATE briscolaFiniteVolume)

# OpenFOAM configuration
openfoam_configure_target(writeRestriction)

# Install
install(TARGETS writeRestriction RUNTIME DESTINATION ${FOAM_USER_APPBIN})
