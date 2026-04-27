# Define executable
add_executable(faceScalar Test-faceScalar.C)

# Dependencies
target_link_libraries(faceScalar PRIVATE briscolaCore)

# OpenFOAM configuration
openfoam_configure_target(faceScalar)

# Add to tests target
add_dependencies(tests faceScalar)

# Add test
add_test(NAME faceScalar COMMAND faceScalar)
