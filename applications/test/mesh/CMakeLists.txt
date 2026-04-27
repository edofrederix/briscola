# Define executable
add_executable(mesh Test-mesh.C)

# Dependencies
target_link_libraries(mesh PRIVATE briscolaMesh)

# OpenFOAM configuration
openfoam_configure_target(mesh)

# Add to tests target
add_dependencies(tests mesh)

# Add tests

set(MESHES 2x 2y 2z 8)
set(NUMPROCS 4 4 4 16)

foreach(MESH NUMPROC IN ZIP_LISTS MESHES NUMPROCS)

    add_test(NAME mesh_${MESH}
        COMMAND
            ${MPIEXEC_EXECUTABLE}
            ${MPIEXEC_NUMPROC_FLAG} ${NUMPROC} --oversubscribe
            $<TARGET_FILE:mesh> -parallel briscolaMeshDict.${MESH}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    set_tests_properties(mesh_${MESH} PROPERTIES
        PROCESSORS ${NUMPROC}
    )

endforeach()
