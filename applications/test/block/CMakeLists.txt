# Define executable
add_executable(block Test-block.C)

# Dependencies
target_link_libraries(block PRIVATE briscolaCore)

# Add to tests target
add_dependencies(tests block)

# Add test
add_test(NAME block COMMAND $<TARGET_FILE:block>)
