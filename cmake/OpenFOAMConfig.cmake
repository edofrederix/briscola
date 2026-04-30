# Apply OpenFOAM-compatible compiler settings to a target

function(openfoam_configure_target TARGET_NAME)

    # Default definitions
    target_compile_definitions(${TARGET_NAME} PRIVATE
        ${OPENFOAM_COMPILE_DEFINITIONS}
        OMPI_SKIP_MPICXX
    )

    # Default include directories
    target_include_directories(${TARGET_NAME} SYSTEM PRIVATE
        ${OPENFOAM_INCLUDE_DIRS}
    )

    # Default link libraries
    target_link_libraries(${TARGET_NAME} PRIVATE
        ${OPENFOAM_LIB}
    )

    # Ensure lnInclude is always updated for libraries

    get_target_property(TARGET_TYPE ${TARGET_NAME} TYPE)

    if(NOT TARGET_TYPE STREQUAL "EXECUTABLE")

        add_custom_target(${TARGET_NAME}_lnInclude ALL
            COMMAND $ENV{WM_PROJECT_DIR}/wmake/wmakeLnInclude -u -s
                ${CMAKE_CURRENT_SOURCE_DIR}
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            COMMENT "Running wmakeLnInclude for ${TARGET_NAME}"
        )

        add_dependencies(${TARGET_NAME} ${TARGET_NAME}_lnInclude)

    endif()

    # Flags
    target_compile_options(${TARGET_NAME} PRIVATE

        # Default OpenFOAM warning flags
        -Wall
        -Wextra
        -Wold-style-cast
        -Wnon-virtual-dtor
        -Wno-unused-parameter
        -Wno-invalid-offsetof
        -Wno-attributes
        -ftemplate-depth-100

        # Custom warning flags
        -Wshadow

        # Build type flags
        $<$<CONFIG:Release>:-O3>
        $<$<CONFIG:Debug>:-g -O0 -DFULLDEBUG>
        $<$<CONFIG:Profile>:-g -O2>

    )

endfunction()
