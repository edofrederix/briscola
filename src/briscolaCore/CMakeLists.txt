# Define library
add_library(briscolaCore SHARED)

# Define sources
target_sources(briscolaCore PRIVATE

    arguments/arguments.C

    primitives/FaceSpace/faceScalar/faceScalar.C
    primitives/FaceSpace/faceLabel/faceLabel.C
    primitives/FaceSpace/faceVector/faceVector.C
    primitives/FaceSpace/faceLabelVector/faceLabelVector.C
    primitives/FaceSpace/faceLabelTensor/faceLabelTensor.C

    primitives/EdgeSpace/edgeScalar/edgeScalar.C
    primitives/EdgeSpace/edgeLabel/edgeLabel.C
    primitives/EdgeSpace/edgeVector/edgeVector.C
    primitives/EdgeSpace/edgeLabelVector/edgeLabelVector.C
    primitives/EdgeSpace/edgeLabelTensor/edgeLabelTensor.C

    primitives/VertexSpace/vertexScalar/vertexScalar.C
    primitives/VertexSpace/vertexLabel/vertexLabel.C
    primitives/VertexSpace/vertexVector/vertexVector.C
    primitives/VertexSpace/vertexLabelVector/vertexLabelVector.C
    primitives/VertexSpace/vertexLabelTensor/vertexLabelTensor.C

    primitives/StencilSpace/stencil/stencil.C
    primitives/StencilSpace/diagStencil/diagStencil.C

    primitives/geometricObjects/tetra.C
    primitives/geometricObjects/hexa.C

    blocks/labelBlock/labelBlock.C
    blocks/scalarBlock/scalarBlock.C

    blocks/vectorBlock/vectorBlock.C
    blocks/labelVectorBlock/labelVectorBlock.C
    blocks/tensorBlock/tensorBlock.C
    blocks/sphericalTensorBlock/sphericalTensorBlock.C
    blocks/symmTensorBlock/symmTensorBlock.C
    blocks/diagTensorBlock/diagTensorBlock.C

    blocks/faceScalarBlock/faceScalarBlock.C
    blocks/faceVectorBlock/faceVectorBlock.C

    blocks/edgeScalarBlock/edgeScalarBlock.C
    blocks/edgeVectorBlock/edgeVectorBlock.C

    blocks/vertexScalarBlock/vertexScalarBlock.C
    blocks/vertexVectorBlock/vertexVectorBlock.C

    blocks/stencilBlock/stencilBlock.C
    blocks/diagStencilBlock/diagStencilBlock.C

    functionObjects/briscolaFunctionObject/briscolaFunctionObject.C
)

# Generate the lnInclude directory to provide wmake-like behavior
add_lnInclude_target(briscolaCore)

# Dependencies
target_link_libraries(briscolaCore PUBLIC OpenFOAM Pstream)

# Add includes directories
target_include_directories(briscolaCore PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/lnInclude
)

# Install
install(TARGETS briscolaCore LIBRARY DESTINATION lib)

# Add to libraries target
add_dependencies(libraries briscolaCore)
