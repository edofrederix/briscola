# Define executable
add_executable(block Test-block.C)

# Dependencies
target_link_libraries(block PRIVATE briscolaCore)

# OpenFOAM configuration
openfoam_configure_target(block)

# Add to tests target
add_dependencies(tests block)

# Add test
add_test(NAME block COMMAND $<TARGET_FILE:block>)
