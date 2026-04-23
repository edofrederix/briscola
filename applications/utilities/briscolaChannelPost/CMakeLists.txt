# Source files

# Define executable
add_executable(briscolaChannelPost ${LIB_SOURCES})

# Define sources
target_sources(briscolaChannelPost PRIVATE briscolaChannelPost.C)

# Dependencies
target_link_libraries(briscolaChannelPost PRIVATE briscolaFiniteVolume)

# OpenFOAM configuration
openfoam_configure_target(briscolaChannelPost)

# Install
install(TARGETS briscolaChannelPost RUNTIME DESTINATION ${FOAM_USER_APPBIN})
