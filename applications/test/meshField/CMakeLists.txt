# Define executable
add_executable(meshField Test-meshField.C)

# Dependencies
target_link_libraries(meshField PRIVATE briscolaFiniteVolume)

# OpenFOAM configuration
openfoam_configure_target(meshField)

# Add to tests target
add_dependencies(tests meshField)

# Add tests

set(MESHES box bigbox pipe 2x 2y 2z)
set(NUMPROCS 4 16 5 4 4 4)

foreach(MESH NUMPROC IN ZIP_LISTS MESHES NUMPROCS)

    add_test(NAME meshField_${MESH}
        COMMAND
            ${MPIEXEC_EXECUTABLE}
            ${MPIEXEC_NUMPROC_FLAG} ${NUMPROC} --oversubscribe
            $<TARGET_FILE:meshField> -parallel briscolaMeshDict.${MESH}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    # Disable sigFpe
    set_tests_properties(meshField_${MESH} PROPERTIES
        PROCESSORS ${NUMPROC}
        ENVIRONMENT_MODIFICATION "FOAM_SIGFPE=unset:;FOAM_ABORT=unset:"
    )

endforeach()
