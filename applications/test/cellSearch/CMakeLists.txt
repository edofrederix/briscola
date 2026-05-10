# Define executable
add_executable(cellSearch Test-cellSearch.C)

# Dependencies
target_link_libraries(cellSearch PRIVATE briscolaFiniteVolume)

set(MESHES box bigbox 2x 2y 2z)
set(NUMPROCS 4 16 4 4 4)

# Add to tests target
add_dependencies(tests cellSearch)

# Add test

foreach(MESH NUMPROC IN ZIP_LISTS MESHES NUMPROCS)

    add_test(NAME cellSearch_${MESH}
        COMMAND
            ${MPIEXEC_EXECUTABLE}
            ${MPIEXEC_NUMPROC_FLAG} ${NUMPROC} --oversubscribe
            $<TARGET_FILE:cellSearch> -parallel briscolaMeshDict.${MESH}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    set_tests_properties(cellSearch_${MESH} PROPERTIES
        PROCESSORS ${NUMPROC}
    )

endforeach()
