# Define executable
add_executable(briscolaStaggeredTwoPhase briscolaStaggeredTwoPhase.C)

# Dependencies
target_link_libraries(briscolaStaggeredTwoPhase PRIVATE briscolaTwoPhase)

# Install
install(TARGETS briscolaStaggeredTwoPhase RUNTIME DESTINATION bin)

# Add to solvers target
add_dependencies(solvers briscolaStaggeredTwoPhase)
