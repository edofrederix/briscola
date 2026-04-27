# Define executable
add_executable(meshDirection Test-meshDirection.C)

# Dependencies
target_link_libraries(meshDirection PRIVATE briscolaFiniteVolume)

# OpenFOAM configuration
openfoam_configure_target(meshDirection)

# Add to tests target
add_dependencies(tests meshDirection)

# Add tests

set(MESHES box bigbox pipe pipeGrad 2x 2y 2z)
set(NUMPROCS 4 16 5 5 4 4 4)

foreach(MESH NUMPROC IN ZIP_LISTS MESHES NUMPROCS)

    add_test(NAME meshDirection_${MESH}
        COMMAND
            ${MPIEXEC_EXECUTABLE}
            ${MPIEXEC_NUMPROC_FLAG} ${NUMPROC} --oversubscribe
            $<TARGET_FILE:meshDirection> -parallel briscolaMeshDict.${MESH}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    # Disable sigFpe
    set_tests_properties(meshDirection_${MESH} PROPERTIES
        PROCESSORS ${NUMPROC}
        ENVIRONMENT_MODIFICATION "FOAM_SIGFPE=unset:"
    )

endforeach()
