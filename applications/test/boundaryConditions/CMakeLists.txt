# Define executable
add_executable(boundaryConditions Test-boundaryConditions.C)

# Dependencies
target_link_libraries(boundaryConditions PRIVATE briscolaFiniteVolume)

# Add to tests target
add_dependencies(tests boundaryConditions)

# Add tests

set(MESHES box bigbox 2x 2y 2z)
set(NUMPROCS 4 16 4 4 4)

foreach(MESH NUMPROC IN ZIP_LISTS MESHES NUMPROCS)

    add_test(NAME boundaryConditions_${MESH}
        COMMAND
            ${MPIEXEC_EXECUTABLE}
            ${MPIEXEC_NUMPROC_FLAG} ${NUMPROC} --oversubscribe
            $<TARGET_FILE:boundaryConditions> -parallel briscolaMeshDict.${MESH}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    set_tests_properties(boundaryConditions_${MESH} PROPERTIES
        PROCESSORS ${NUMPROC}
    )

endforeach()

