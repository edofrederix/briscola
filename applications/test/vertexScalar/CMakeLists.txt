# Define executable
add_executable(vertexScalar Test-vertexScalar.C)

# Dependencies
target_link_libraries(vertexScalar PRIVATE briscolaCore)

# OpenFOAM configuration
openfoam_configure_target(vertexScalar)

# Add to tests target
add_dependencies(tests vertexScalar)

# Add test
add_test(NAME vertexScalar COMMAND vertexScalar)
