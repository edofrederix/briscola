# Define executable
add_executable(briscolaCheckMesh briscolaCheckMesh.C)

# Dependencies
target_link_libraries(briscolaCheckMesh PRIVATE briscolaFiniteVolume)

# Install
install(TARGETS briscolaCheckMesh RUNTIME DESTINATION bin)

# Add to utilities target
add_dependencies(utilities briscolaCheckMesh)
