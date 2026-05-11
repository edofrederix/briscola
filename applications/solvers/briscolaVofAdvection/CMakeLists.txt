# Define executable
add_executable(briscolaVofAdvection briscolaVofAdvection.C)

# Dependencies
target_link_libraries(briscolaVofAdvection PRIVATE briscolaTwoPhase)

# Install
install(TARGETS briscolaVofAdvection RUNTIME DESTINATION bin)

# Add to solvers target
add_dependencies(solvers briscolaVofAdvection)
