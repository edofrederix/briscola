# Define executable
add_executable(vertexVector Test-vertexVector.C)

# Dependencies
target_link_libraries(vertexVector PRIVATE briscolaCore)

# Add to tests target
add_dependencies(tests vertexVector)

# Add test
add_test(NAME vertexVector COMMAND vertexVector)
