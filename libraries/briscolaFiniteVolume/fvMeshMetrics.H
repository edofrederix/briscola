#ifndef fvMeshMetrics_H
#define fvMeshMetrics_H

#include "colocatedFieldsFwd.H"
#include "staggeredFieldsFwd.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class fvMesh;

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

    // Face delta coefficients (inverse of the cell-to-cell distance across a
    // face)
    meshField<faceScalar,MeshType> faceDeltas_;

    // Face interpolation weights for the cell center
    meshField<faceScalar,MeshType> faceWeightsCenter_;

    // Face interpolation weights for the neighbor
    meshField<faceScalar,MeshType> faceWeightsNeighbor_;

    // Global cell numbers
    meshField<label,MeshType> globalCellNumbers_;

    // Calculate face centers
    void calculateFaceCenters();

    // Calculate edge centers
    void calculateEdgeCenters();

    // Calculate vertex centers
    void calculateVertexCenters();

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

public:

    // Declare name of the class and its debug switch
    ClassName("fvMeshMetrics");

    // Construct from fvMesh
    fvMeshMetrics(const fvMesh&);

    // Destructor
    ~fvMeshMetrics();

    // Disallow default bitwise assignment
    void operator=(const fvMeshMetrics&);

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
};

}

}

}

#endif
