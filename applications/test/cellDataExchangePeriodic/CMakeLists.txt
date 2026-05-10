# Define executable
add_executable(cellDataExchangePeriodic Test-cellDataExchangePeriodic.C)

# Dependencies
target_link_libraries(cellDataExchangePeriodic PRIVATE briscolaFiniteVolume)

# Add to tests target
add_dependencies(tests cellDataExchangePeriodic)

# Add test

foreach(NUMPROC 1 2 3)

    math(EXPR NUMPROC_CUBED "${NUMPROC} * ${NUMPROC} * ${NUMPROC}")

    if(NUMPROC_CUBED EQUAL 1)
        add_test(NAME cellDataExchangePeriodic_${NUMPROC}
            COMMAND
                $<TARGET_FILE:cellDataExchangePeriodic>
                briscolaMeshDict.${NUMPROC}
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )
    else()
        add_test(NAME cellDataExchangePeriodic_${NUMPROC}
            COMMAND
                ${MPIEXEC_EXECUTABLE}
                ${MPIEXEC_NUMPROC_FLAG} ${NUMPROC_CUBED} --oversubscribe
                $<TARGET_FILE:cellDataExchangePeriodic> -parallel
                briscolaMeshDict.${NUMPROC}
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )
    endif()

    set_tests_properties(cellDataExchangePeriodic_${NUMPROC} PROPERTIES
        PROCESSORS ${NUMPROC_CUBED}
    )

endforeach()
