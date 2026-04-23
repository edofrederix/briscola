# Source files

# Define library
add_library(briscolaMesh SHARED ${LIB_SOURCES})

# Define sources
target_sources(briscolaMesh PRIVATE
    mesh/meshData.C
    mesh/meshObject.C

    edges/edge/edge.C
    edges/edgeKey/edgeKey.C
    edges/lineEdge/lineEdge.C
    edges/arcEdge/arcEdge.C

    gradings/grading/grading.C
    gradings/noGrading/noGrading.C
    gradings/geometricGrading/geometricGrading.C
    gradings/geometricEdgeGrading/geometricEdgeGrading.C

    vertex/vertex.C
    face/face.C

    brick/brick.C
    brick/brickLink.C
    brick/brickLinks.C
    brick/brickFaceLink.C
    brick/brickEdgeLink.C
    brick/brickVertexLink.C
    brick/brickTopology.C
    brick/brickTopologyMap.C

    patches/patch/patch.C
    patches/mappedPatch/mappedPatch.C
    patches/emptyPatch/emptyPatch.C
    patches/periodicPatch/periodicPatch.C

    patchPair/patchPair.C

    geometry/geometry.C

    decompositions/decomposition/decomposition.C
    decompositions/decomposition/decompositionMap.C
    decompositions/manualDecomp/manualDecomp.C

    decompositions/decompositionInterface/decompositionInterface.C

    boundaries/boundary/boundary.C
    boundaries/parallelBoundary/parallelBoundary.C
    boundaries/periodicBoundary/periodicBoundary.C
    boundaries/mappedBoundary/mappedBoundary.C
    boundaries/patchBoundary/patchBoundary.C
    boundaries/emptyBoundary/emptyBoundary.C

    level/level.C
    level/levelBoundaries.C
    level/levelComms.C
    level/levelPoints.C

    mesh/unstructuredMesh/unstructuredMesh.C
    mesh/structuredMesh/structuredMesh.C
    mesh/rectilinearMesh/rectilinearMesh.C
    mesh/uniformMesh/uniformMesh.C
    mesh/mesh.C
)

# Dependencies
target_link_libraries(briscolaMesh PUBLIC briscolaCore MPI::MPI_CXX)

# Add local includes
target_include_directories(briscolaMesh PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/lnInclude
)

# OpenFOAM configuration
openfoam_configure_target(briscolaMesh)

# Install
install(TARGETS briscolaMesh LIBRARY DESTINATION ${FOAM_USER_LIBBIN})
