# Define executable
add_executable(vertexScalar Test-vertexScalar.C)

# Dependencies
target_link_libraries(vertexScalar PRIVATE briscolaCore)

# Add to tests target
add_dependencies(tests vertexScalar)

# Add test
add_test(NAME vertexScalar COMMAND vertexScalar)
