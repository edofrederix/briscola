# Define executable
add_executable(pointDataExchange Test-pointDataExchange.C)

# Dependencies
target_link_libraries(pointDataExchange PRIVATE briscolaFiniteVolume)

# OpenFOAM configuration
openfoam_configure_target(pointDataExchange)

# Add to tests target
add_dependencies(tests pointDataExchange)

# Add tests

set(MESHES box bigbox 2x 2y 2z)
set(NUMPROCS 4 16 4 4 4)

foreach(MESH NUMPROC IN ZIP_LISTS MESHES NUMPROCS)

    add_test(NAME pointDataExchange_${MESH}
        COMMAND
            ${MPIEXEC_EXECUTABLE}
            ${MPIEXEC_NUMPROC_FLAG} ${NUMPROC} --oversubscribe
            $<TARGET_FILE:pointDataExchange> -parallel briscolaMeshDict.${MESH}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    set_tests_properties(pointDataExchange_${MESH} PROPERTIES
        PROCESSORS ${NUMPROC}
    )

endforeach()
