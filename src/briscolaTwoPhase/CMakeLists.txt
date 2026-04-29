# Define library
add_library(briscolaTwoPhase SHARED)

# Define sources
target_sources(briscolaTwoPhase PRIVATE
    vof/vof/vof.C
    vof/geometricVof/geometricVof.C
    vof/splitAdvection/splitAdvection.C

    tagAlgorithms/tagAlgorithm/tagAlgorithm.C
    tagAlgorithms/twoPass/twoPass.C

    vofField/vofField.C

    LVE/LVE.C

    normalSchemes/normalScheme/normalScheme.C
    normalSchemes/Youngs/Youngs.C
    normalSchemes/LSGIR/LSGIR.C
    normalSchemes/LSFIR/LSFIR.C
    normalSchemes/MYC/MYC.C

    twoPhaseModels/twoPhaseModel/TwoPhaseModel.C
    twoPhaseModels/twoPhaseModel/twoPhaseModel.C
    twoPhaseModels/twoPhaseModel/twoPhaseModels.C

    surfaceTensionSchemes/surfaceTensionScheme/surfaceTensionScheme.C
    surfaceTensionSchemes/sigmaModels/constantSigma/constantSigma.C
    surfaceTensionSchemes/sigmaModels/zeroSigma/zeroSigma.C
    surfaceTensionSchemes/surfaceTensionScheme/surfaceTensionSchemes.C

    curvatureSchemes/curvatureScheme/curvatureScheme.C
    curvatureSchemes/CV/CV.C
    curvatureSchemes/SHF/SHF.C
    curvatureSchemes/parabolicFit/parabolicFit.C

    limiterSchemes/twoPhaseLimiterScheme/twoPhaseLimiterSchemes.C
)

# Dependencies
target_link_libraries(briscolaTwoPhase PUBLIC briscolaFiniteVolume)

# Add local includes
target_include_directories(briscolaTwoPhase PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/lnInclude
)

# OpenFOAM configuration
openfoam_configure_target(briscolaTwoPhase)

# Install
install(TARGETS briscolaTwoPhase LIBRARY DESTINATION lib)

# Add to libraries target
add_dependencies(libraries briscolaTwoPhase)
