# Define executable
add_executable(linearSystemAggregation Test-linearSystemAggregation.C)

# Dependencies
target_link_libraries(linearSystemAggregation PRIVATE briscolaFiniteVolume)

if(EIGEN3_FOUND)
    target_link_libraries(linearSystemAggregation PRIVATE Eigen3)
endif()

# OpenFOAM configuration
openfoam_configure_target(linearSystemAggregation)

# Add to tests target
add_dependencies(tests linearSystemAggregation)

# Add test
add_test(NAME linearSystemAggregation
    COMMAND
        ${MPIEXEC_EXECUTABLE}
        ${MPIEXEC_NUMPROC_FLAG} 6 --oversubscribe
        $<TARGET_FILE:linearSystemAggregation> -parallel
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

set_tests_properties(linearSystemAggregation PROPERTIES PROCESSORS 6)

# Add cleanup target
add_custom_target(clean-tests-linearSystemAggregation
    COMMAND find ${CMAKE_CURRENT_SOURCE_DIR} -maxdepth 1
        -name "f-colocated-*" -delete
    COMMAND find ${CMAKE_CURRENT_SOURCE_DIR} -maxdepth 1
        -name "f-staggered-*" -delete
)

add_dependencies(clean-tests clean-tests-linearSystemAggregation)
