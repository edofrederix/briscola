# Define executable
add_executable(stencils Test-stencils.C)

# Dependencies
target_link_libraries(stencils PRIVATE briscolaCore)

# OpenFOAM configuration
openfoam_configure_target(stencils)

# Add to tests target
add_dependencies(tests stencils)

# Add test
add_test(NAME stencils COMMAND $<TARGET_FILE:stencils>)
