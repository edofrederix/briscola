# Define executable
add_executable(briscolaColocatedTwoPhase briscolaColocatedTwoPhase.C)

# Dependencies
target_link_libraries(briscolaColocatedTwoPhase PRIVATE briscolaTwoPhase)

# Install
install(TARGETS briscolaColocatedTwoPhase RUNTIME DESTINATION bin)

# Add to solvers target
add_dependencies(solvers briscolaColocatedTwoPhase)
