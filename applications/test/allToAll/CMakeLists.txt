# Define executable
add_executable(allToAll Test-allToAll.C)

# Dependencies
target_link_libraries(allToAll PRIVATE briscolaCore)

# Add to tests target
add_dependencies(tests allToAll)

# Add test
add_test(NAME allToAll
    COMMAND
        ${MPIEXEC_EXECUTABLE}
        ${MPIEXEC_NUMPROC_FLAG} 18 --oversubscribe
        $<TARGET_FILE:allToAll> -parallel
        "(32 48 32)"
        "(9 2 1)"
        "(3 3 2)"
        2
)

set_tests_properties(allToAll PROPERTIES PROCESSORS 18)
