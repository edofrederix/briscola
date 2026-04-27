# Define executable
add_executable(edgeVector Test-edgeVector.C)

# Dependencies
target_link_libraries(edgeVector PRIVATE briscolaCore)

# OpenFOAM configuration
openfoam_configure_target(edgeVector)

# Add to tests target
add_dependencies(tests edgeVector)

# Add test
add_test(NAME edgeVector COMMAND edgeVector)
