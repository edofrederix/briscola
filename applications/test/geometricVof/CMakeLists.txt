# Define executable
add_executable(geometricVof Test-geometricVof.C)

# Dependencies
target_link_libraries(geometricVof PRIVATE briscolaTwoPhase)

# OpenFOAM configuration
openfoam_configure_target(geometricVof)

# Add to tests target
add_dependencies(tests geometricVof)

# Add tests

set(MESHES 2x 2y 2z 8)
set(NUMPROCS 4 4 4 16)

foreach(MESH NUMPROC IN ZIP_LISTS MESHES NUMPROCS)

    add_test(NAME geometricVof_${MESH}
        COMMAND
            ${MPIEXEC_EXECUTABLE}
            ${MPIEXEC_NUMPROC_FLAG} ${NUMPROC} --oversubscribe
            $<TARGET_FILE:geometricVof> -parallel briscolaMeshDict.${MESH}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    set_tests_properties(geometricVof_${MESH} PROPERTIES
        PROCESSORS ${NUMPROC}
    )

endforeach()
