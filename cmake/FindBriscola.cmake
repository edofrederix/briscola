# Locate Briscola installation via environment variables
if(NOT DEFINED ENV{BRISCOLA} OR "$ENV{BRISCOLA}" STREQUAL "")
    message(FATAL_ERROR "BRISCOLA environment variable not set")
endif()

# Require OpenFOAM
find_package(OpenFOAM REQUIRED)

# Required include paths from Briscola
set(BRISCOLA_INCLUDE_DIRS
    $ENV{BRISCOLA}/src/briscolaCore/lnInclude
    $ENV{BRISCOLA}/src/briscolaMesh/lnInclude
    $ENV{BRISCOLA}/src/briscolaFiniteVolume/lnInclude
    $ENV{BRISCOLA}/src/briscolaTwoPhase/lnInclude
)

# Find Briscola
find_library(BRISCOLA_LIB_FV briscolaFiniteVolume HINTS $ENV{FOAM_USER_LIBBIN})
find_library(BRISCOLA_LIB_TP briscolaTwoPhase HINTS $ENV{FOAM_USER_LIBBIN})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Briscola DEFAULT_MSG
    BRISCOLA_LIB_FV BRISCOLA_LIB_TP)

# Set Briscola target. It depends on the fv, two-phase and OpenFOAM libraries.

if(BRISCOLA_FOUND AND NOT TARGET Briscola)
    add_library(Briscola INTERFACE IMPORTED)

    set_target_properties(Briscola PROPERTIES
        INTERFACE_LINK_LIBRARIES
            "${BRISCOLA_LIB_FV};${BRISCOLA_LIB_TP};OpenFOAM"
        INTERFACE_INCLUDE_DIRECTORIES "${BRISCOLA_INCLUDE_DIRS}"
    )
endif()
