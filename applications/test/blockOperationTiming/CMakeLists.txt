# Define executable
add_executable(blockOperationTiming Test-blockOperationTiming.C)

# Dependencies
target_link_libraries(blockOperationTiming PRIVATE briscolaFiniteVolume)

# Add to tests target
add_dependencies(tests blockOperationTiming)

# Add test
add_test(NAME blockOperationTiming
    COMMAND blockOperationTiming
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
