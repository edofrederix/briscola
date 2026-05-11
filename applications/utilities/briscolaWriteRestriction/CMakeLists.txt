# Define executable
add_executable(briscolaWriteRestriction briscolaWriteRestriction.C)

# Dependencies
target_link_libraries(briscolaWriteRestriction PRIVATE briscolaFiniteVolume)

# Install
install(TARGETS briscolaWriteRestriction RUNTIME DESTINATION bin)

# Add to utilities target
add_dependencies(utilities briscolaWriteRestriction)
