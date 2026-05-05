add_custom_target(tests)
add_custom_target(clean-tests)

# Find Python3 and check for numpy availability -- needed for some tests

find_package(Python3 COMPONENTS Interpreter)

if(Python3_FOUND)
    execute_process(
        COMMAND ${Python3_EXECUTABLE} -c "import numpy"
        RESULT_VARIABLE NUMPY_RESULT
        OUTPUT_QUIET
        ERROR_QUIET
    )
endif()

# Add all sub-directories with a CMakeLists.txt file. Exclude tests from all.

file(GLOB dirs LIST_DIRECTORIES true "./*/")

foreach(dir ${dirs})
  if(EXISTS "${dir}/CMakeLists.txt")
    add_subdirectory(${dir} EXCLUDE_FROM_ALL)
  endif()
endforeach()
