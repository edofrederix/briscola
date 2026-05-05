# Define executable
add_executable(externalSolver Test-externalSolver.C)

# Dependencies
target_link_libraries(externalSolver PRIVATE briscolaFiniteVolume)

# Add to tests target
add_dependencies(tests externalSolver)

# Add tests

add_test(NAME externalSolver
    COMMAND
        ${MPIEXEC_EXECUTABLE}
        ${MPIEXEC_NUMPROC_FLAG} 4 --oversubscribe
        $<TARGET_FILE:externalSolver> -parallel
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

set_tests_properties(externalSolver PROPERTIES PROCESSORS 4)

if(Python3_FOUND AND NUMPY_RESULT EQUAL 0)
    add_test(NAME externalSolver_python
        COMMAND ${Python3_EXECUTABLE} check.py
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
else()
    message(STATUS "Skipping check_output test: "
        "Python3=${Python3_FOUND}, "
        "numpy=${NUMPY_RESULT}"
    )
endif()

set_tests_properties(externalSolver PROPERTIES PROCESSORS 8)

# Add cleanup target
add_custom_target(clean-tests-externalSolver
    COMMAND find ${CMAKE_CURRENT_SOURCE_DIR} -maxdepth 1
        -name "f-colocated-*" -delete
    COMMAND find ${CMAKE_CURRENT_SOURCE_DIR} -maxdepth 1
        -name "f-staggered-*" -delete
)

add_dependencies(clean-tests clean-tests-externalSolver)
