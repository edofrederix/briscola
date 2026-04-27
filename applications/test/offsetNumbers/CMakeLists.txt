# Define executable
add_executable(offsetNumbers Test-offsetNumbers.C)

# Dependencies
target_link_libraries(offsetNumbers PRIVATE briscolaCore)

# OpenFOAM configuration
openfoam_configure_target(offsetNumbers)

# Add to tests target
add_dependencies(tests offsetNumbers)

# Add test
add_test(NAME offsetNumbers COMMAND $<TARGET_FILE:offsetNumbers>)
