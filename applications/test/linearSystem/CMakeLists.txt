# Define executable
add_executable(linearSystem Test-linearSystem.C)

# Dependencies
target_link_libraries(linearSystem PRIVATE briscolaFiniteVolume)

# OpenFOAM configuration
openfoam_configure_target(linearSystem)

# Add to tests target
add_dependencies(tests linearSystem)

# Add tests

set(MESHES box bigbox pipe pipeGrad 2x 2y 2z)
set(NUMPROCS 4 16 5 5 4 4 4)

foreach(MESH NUMPROC IN ZIP_LISTS MESHES NUMPROCS)

    add_test(NAME linearSystem_${MESH}
        COMMAND
            ${MPIEXEC_EXECUTABLE}
            ${MPIEXEC_NUMPROC_FLAG} ${NUMPROC} --oversubscribe
            $<TARGET_FILE:linearSystem> -parallel briscolaMeshDict.${MESH}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    set_tests_properties(linearSystem_${MESH} PROPERTIES
        PROCESSORS ${NUMPROC}
    )

endforeach()
