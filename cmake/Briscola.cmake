set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Build types

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type")
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
  "Debug" "Release" "Profile"
)

# Store built libraries and applications in lib and bin in the build directory

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Helper modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")

# Default prefix points to OpenFOAM's user install location

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    cmake_path(SET install_path NORMALIZE "$ENV{FOAM_USER_LIBBIN}/../")
    set(CMAKE_INSTALL_PREFIX ${install_path} CACHE PATH "Install prefix" FORCE)
    message(STATUS "Using OpenFOAM's user install location as prefix:
        ${install_path}")
endif()

# Function to add an additional lnInclude target to the given target

function(add_lnInclude_target TARGET_NAME)

    add_custom_target(${TARGET_NAME}_lnInclude ALL
        COMMAND $ENV{WM_PROJECT_DIR}/wmake/wmakeLnInclude -u -s
            ${CMAKE_CURRENT_SOURCE_DIR}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMENT "Running wmakeLnInclude for ${TARGET_NAME}"
    )

    add_dependencies(${TARGET_NAME} ${TARGET_NAME}_lnInclude)

endfunction()
