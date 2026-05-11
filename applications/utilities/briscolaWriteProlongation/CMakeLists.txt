# Define executable
add_executable(briscolaWriteProlongation briscolaWriteProlongation.C)

# Dependencies
target_link_libraries(briscolaWriteProlongation PRIVATE briscolaFiniteVolume)

# Install
install(TARGETS briscolaWriteProlongation RUNTIME DESTINATION bin)

# Add to utilities target
add_dependencies(utilities briscolaWriteProlongation)
