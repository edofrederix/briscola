#ifndef fvMeshMetrics_H
#define fvMeshMetrics_H

#include "immersedBoundaryFwd.H"
#include "labelList.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class fvMesh;

/**

@brief Class that defines metrics of the mesh

This class generates important mesh metrics, such as cell/face/edge/vertex
centers, normals, volumes, weights and numbers. Also immersed boundaries are
generated. The metric fields can then be accessed by everything else, through
the fvMesh object of which this class is part.

*/

template<class MeshType>
class fvMeshMetrics
{
    // Reference to the fvMesh object
    const fvMesh& fvMsh_;

    // Face center field
    meshField<faceVector,MeshType> faceCenters_;

    // Edge center field
    meshField<edgeVector,MeshType> edgeCenters_;

    // Vertex center field
    meshField<vertexVector,MeshType> vertexCenters_;

    // Face outward normal field
    meshField<faceVector,MeshType> faceNormals_;

    // Face area field
    meshField<faceScalar,MeshType> faceAreas_;

    // Face outward normal with magnitude face area field
    meshField<faceVector,MeshType> faceAreaNormals_;

    // Cell center field
    meshField<vector,MeshType> cellCenters_;

    // Cell volume field
    meshField<scalar,MeshType> cellVolumes_;

    // Inverse cell volume field
    meshField<scalar,MeshType> inverseCellVolumes_;

    // Face delta coefficients (inverse of the cell-to-cell distance across a
    // face)
    meshField<faceScalar,MeshType> faceDeltas_;

    // Face interpolation weights for the cell center
    meshField<faceScalar,MeshType> faceWeightsCenter_;

    // Face interpolation weights for the neighbor
    meshField<faceScalar,MeshType> faceWeightsNeighbor_;

    // Global cell numbers
    meshField<label,MeshType> globalCellNumbers_;

    // Global cell counts per level
    List<labelList> globalCellCounts_;

    // List of immersed boundaries (if any)
    PtrList<immersedBoundary<MeshType>> immersedBoundaries_;

    // Global immersed boundary mask
    meshField<scalar,MeshType> mask_;

    // Calculate vertex centers
    void calculateVertexCenters();

    // Calculate face centers
    void calculateFaceCenters();

    // Calculate edge centers
    void calculateEdgeCenters();

    // Calculate face areas and normals
    void calculateFaceAreasAndNormals();

    // Calculate cell centers
    void calculateCellCenters();

    // Calculate cell volumes
    void calculateCellVolumes();

    // Calculate face deltas
    void calculateFaceDeltas();

    // Calculate face weights
    void calculateFaceWeights();

    // Set the global cell numbers
    void setGlobalCellNumbers();

    // Set the immersed boundaries
    void setImmersedBoundaries();

    // Set the immersed boundary mask
    void setImmersedBoundaryMask();

public:

    // Declare name of the class and its debug switch
    ClassName("fvMeshMetrics");

    // Construct from fvMesh
    fvMeshMetrics(const fvMesh&);

    // Destructor
    ~fvMeshMetrics();

    // Const access to immersed boundaries
    const PtrList<immersedBoundary<MeshType>>& immersedBoundaries() const
    {
        return immersedBoundaries_;
    }

    // Disallow default bitwise assignment
    void operator=(const fvMeshMetrics&);

    // Const access to fvMesh
    inline const fvMesh& fvMsh() const
    {
        return fvMsh_;
    }

    // Const access to face centers
    inline const meshField<faceVector,MeshType>& faceCenters() const
    {
        return faceCenters_;
    }

    // Const access to edge centers
    inline const meshField<edgeVector,MeshType>& edgeCenters() const
    {
        return edgeCenters_;
    }

    // Const access to vertex centers
    inline const meshField<vertexVector,MeshType>& vertexCenters() const
    {
        return vertexCenters_;
    }

    // Const access to face outward normals
    inline const meshField<faceVector,MeshType>& faceNormals() const
    {
        return faceNormals_;
    }

    // Const access to face areas
    inline const meshField<faceScalar,MeshType>& faceAreas() const
    {
        return faceAreas_;
    }

    // Const access to face outward normal with area magnitude
    inline const meshField<faceVector,MeshType>& faceAreaNormals() const
    {
        return faceAreaNormals_;
    }

    // Const access to cell centers
    inline const meshField<vector,MeshType>& cellCenters() const
    {
        return cellCenters_;
    }

    // Const access to cell volumes
    inline const meshField<scalar,MeshType>& cellVolumes() const
    {
        return cellVolumes_;
    }

    // Const access to inverse cell volumes
    inline const meshField<scalar,MeshType>& inverseCellVolumes() const
    {
        return inverseCellVolumes_;
    }

    // Const access to face deltas
    inline const meshField<faceScalar,MeshType>& faceDeltas() const
    {
        return faceDeltas_;
    }

    // Const access to face weights for the cell center
    inline const meshField<faceScalar,MeshType>& faceWeightsCenter() const
    {
        return faceWeightsCenter_;
    }

    // Const access to face weights for the neighbor
    inline const meshField<faceScalar,MeshType>& faceWeightsNeighbor() const
    {
        return faceWeightsNeighbor_;
    }

    // Const access to the global cell numbers
    inline const meshField<label,MeshType>& globalCellNumbers() const
    {
        return globalCellNumbers_;
    }

    // Global number of cells
    inline const List<labelList>& globalCellCounts() const
    {
        return globalCellCounts_;
    }

    // Immersed boundary mask
    inline const meshField<scalar,MeshType>& mask() const
    {
        return mask_;
    }

    // Fluid mask
    inline tmp<meshField<scalar,MeshType>> fluidMask() const
    {
        return 1.0 - mask_;
    }

    // Solid mask
    inline tmp<meshField<scalar,MeshType>> solidMask() const
    {
        return mask_;
    }
};

}

}

}

#endif
