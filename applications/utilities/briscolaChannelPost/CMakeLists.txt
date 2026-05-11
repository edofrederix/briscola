# Define executable
add_executable(briscolaChannelPost briscolaChannelPost.C)

# Dependencies
target_link_libraries(briscolaChannelPost PRIVATE briscolaFiniteVolume)

# Install
install(TARGETS briscolaChannelPost RUNTIME DESTINATION bin)

# Add to utilities target
add_dependencies(utilities briscolaChannelPost)
