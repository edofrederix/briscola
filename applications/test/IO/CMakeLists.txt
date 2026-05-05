# Define executable
add_executable(IO Test-IO.C)

# Dependencies
target_link_libraries(IO PRIVATE briscolaTwoPhase)

# Add to tests target
add_dependencies(tests IO)

# Add tests

# set(FORMATS ascii binary)
set(FORMATS binary)
set(MESHES box pipe pipeGrad bigbox 2x 2y 2z)
set(NUMPROCS 4 5 5 16 4 4 4)

# IO tests need to run sequentially because they write data to the same files.
# Let each test depend on the previous one.

set(prev_test "")

foreach(FORMAT IN LISTS FORMATS)
    foreach(MESH NUMPROC IN ZIP_LISTS MESHES NUMPROCS)

        add_test(NAME IO_${MESH}_${FORMAT}
            COMMAND
                ${MPIEXEC_EXECUTABLE}
                ${MPIEXEC_NUMPROC_FLAG} ${NUMPROC} --oversubscribe
                $<TARGET_FILE:IO>
                -parallel briscolaMeshDict.${MESH} ${FORMAT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )

        set_tests_properties(IO_${MESH}_${FORMAT} PROPERTIES
            PROCESSORS ${NUMPROC}
        )

        if(prev_test)
            set_tests_properties(IO_${MESH}_${FORMAT} PROPERTIES
                DEPENDS ${prev_test}
            )
        endif()

        set(prev_test IO_${MESH}_${FORMAT})

    endforeach()
endforeach()

# Add cleanup target
add_custom_target(clean-tests-IO
    COMMAND find ${CMAKE_CURRENT_SOURCE_DIR} -maxdepth 1
        -name "briscola_*.series" -delete
    COMMAND find ${CMAKE_CURRENT_SOURCE_DIR} -maxdepth 1
        -name "[1-9]" -delete
)

add_dependencies(clean-tests clean-tests-IO)
