# Define executable
add_executable(briscolaColocated briscolaColocated.C)

# Dependencies
target_link_libraries(briscolaColocated PRIVATE briscolaFiniteVolume)

# Install
install(TARGETS briscolaColocated RUNTIME DESTINATION bin)

# Add to solvers target
add_dependencies(solvers briscolaColocated)
