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

    // Face center field pointer
    meshField<hexVector,MeshType> faceCenters_;

    // Face normal field pointer
    meshField<hexVector,MeshType> faceNormals_;

    // Face area field pointer
    meshField<hexScalar,MeshType> faceAreas_;

    // Face normals with magnitude face area field pointer
    meshField<hexVector,MeshType> faceAreaNormals_;

    // Cell center field pointer
    meshField<vector,MeshType> cellCenters_;

    // Cell volume field pointer
    meshField<scalar,MeshType> cellVolumes_;

    // Face delta coefficients (inverse of the cell-to-cell distance across a
    // face)
    meshField<hexScalar,MeshType> faceDeltas_;

    // Calculate face centers
    void calculateFaceCenters();

    // Calculate face areas and normals
    void calculateFaceAreasAndNormals();

    // Calculate cell centers
    void calculateCellCenters();

    // Calculate cell volumes
    void calculateCellVolumes();

    // Calculate face deltas
    void calculateFaceDeltas();

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
    inline const meshField<hexVector,MeshType>& faceCenters() const
    {
        return faceCenters_;
    }

    // Const access to face normals
    inline const meshField<hexVector,MeshType>& faceNormals() const
    {
        return faceNormals_;
    }

    // Const access to face areas
    inline const meshField<hexScalar,MeshType>& faceAreas() const
    {
        return faceAreas_;
    }

    // Const access to face areas
    inline const meshField<hexVector,MeshType>& faceAreaNormals() const
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
    inline const meshField<hexScalar,MeshType>& faceDeltas() const
    {
        return faceDeltas_;
    }
};

}

}

}

#endif
