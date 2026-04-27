# Define executable
add_executable(edgeScalar Test-edgeScalar.C)

# Dependencies
target_link_libraries(edgeScalar PRIVATE briscolaCore)

# OpenFOAM configuration
openfoam_configure_target(edgeScalar)

# Add to tests target
add_dependencies(tests edgeScalar)

# Add test
add_test(NAME edgeScalar COMMAND edgeScalar)
