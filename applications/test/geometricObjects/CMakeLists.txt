# Define executable
add_executable(geometricObjects Test-geometricObjects.C)

# Dependencies
target_link_libraries(geometricObjects PRIVATE briscolaTwoPhase)

# OpenFOAM configuration
openfoam_configure_target(geometricObjects)

# Add to tests target
add_dependencies(tests geometricObjects)

# Add test
add_test(NAME geometricObjects COMMAND $<TARGET_FILE:geometricObjects>)
