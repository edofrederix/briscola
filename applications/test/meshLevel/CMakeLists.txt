# Define executable
add_executable(meshLevel Test-meshLevel.C)

# Dependencies
target_link_libraries(meshLevel PRIVATE briscolaFiniteVolume)

# OpenFOAM configuration
openfoam_configure_target(meshLevel)

# Add to tests target
add_dependencies(tests meshLevel)

# Add tests

set(MESHES box bigbox pipe pipeGrad 2x 2y 2z)
set(NUMPROCS 4 16 5 5 4 4 4)

foreach(MESH NUMPROC IN ZIP_LISTS MESHES NUMPROCS)

    add_test(NAME meshLevel_${MESH}
        COMMAND
            ${MPIEXEC_EXECUTABLE}
            ${MPIEXEC_NUMPROC_FLAG} ${NUMPROC} --oversubscribe
            $<TARGET_FILE:meshLevel> -parallel briscolaMeshDict.${MESH}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    # Disable sigFpe
    set_tests_properties(meshLevel_${MESH} PROPERTIES
        PROCESSORS ${NUMPROC}
        ENVIRONMENT_MODIFICATION "FOAM_SIGFPE=unset:"
    )

endforeach()
