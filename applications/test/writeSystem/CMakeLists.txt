# Define executable
add_executable(writeSystem Test-writeSystem.C)

# Dependencies
target_link_libraries(writeSystem PRIVATE briscolaFiniteVolume)

# Add to tests target
add_dependencies(tests writeSystem)

# Add test
add_test(NAME writeSystem
    COMMAND
        ${MPIEXEC_EXECUTABLE}
        ${MPIEXEC_NUMPROC_FLAG} 8 --oversubscribe
        $<TARGET_FILE:writeSystem> -parallel
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

set_tests_properties(writeSystem PROPERTIES PROCESSORS 8)

# Add cleanup target
add_custom_target(clean-tests-writeSystem
    COMMAND find ${CMAKE_CURRENT_SOURCE_DIR} -maxdepth 1
        -name "f-colocated-*" -delete
    COMMAND find ${CMAKE_CURRENT_SOURCE_DIR} -maxdepth 1
        -name "f-staggered-*" -delete
)

add_dependencies(clean-tests clean-tests-writeSystem)
