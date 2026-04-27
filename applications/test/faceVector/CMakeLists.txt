# Define executable
add_executable(faceVector Test-faceVector.C)

# Dependencies
target_link_libraries(faceVector PRIVATE briscolaCore)

# OpenFOAM configuration
openfoam_configure_target(faceVector)

# Add to tests target
add_dependencies(tests faceVector)

# Add test
add_test(NAME faceVector COMMAND faceVector)
