add_custom_target(libraries)

# Add all sub-directories with a CMakeLists.txt file

file(GLOB dirs LIST_DIRECTORIES true "./*/")

foreach(dir ${dirs})
  if(EXISTS "${dir}/CMakeLists.txt")
    add_subdirectory(${dir})
  endif()
endforeach()
