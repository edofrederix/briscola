# Define executable
add_executable(briscolaStaggered briscolaStaggered.C)

# Dependencies
target_link_libraries(briscolaStaggered PRIVATE briscolaFiniteVolume)

# Install
install(TARGETS briscolaStaggered RUNTIME DESTINATION bin)

# Add to solvers target
add_dependencies(solvers briscolaStaggered)
