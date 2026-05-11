# Define executable
add_executable(offsetNumbers Test-offsetNumbers.C)

# Dependencies
target_link_libraries(offsetNumbers PRIVATE briscolaCore)

# Add to tests target
add_dependencies(tests offsetNumbers)

# Add test
add_test(NAME offsetNumbers COMMAND $<TARGET_FILE:offsetNumbers>)
