add_custom_target(tests)
add_custom_target(clean-tests)

# Add all sub-directories with a CMakeLists.txt file. Exclude tests from all.

file(GLOB dirs LIST_DIRECTORIES true "./*/")

foreach(dir ${dirs})
  if(EXISTS "${dir}/CMakeLists.txt")
    add_subdirectory(${dir} EXCLUDE_FROM_ALL)
  endif()
endforeach()
