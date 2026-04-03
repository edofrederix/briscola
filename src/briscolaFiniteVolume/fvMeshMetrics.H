#ifndef fvMeshMetrics_H
#define fvMeshMetrics_H

#include "labelList.H"
#include "FastPtrList.H"

#include "colocatedFieldsFwd.H"
#include "staggeredFieldsFwd.H"
#include "faceFieldFwd.H"
#include "immersedBoundaryFwd.H"

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

    // Cell metric field pointers
    meshField<vector,MeshType>* cellCentersPtr_;
    meshField<scalar,MeshType>* cellVolumesPtr_;
    meshField<scalar,MeshType>* inverseCellVolumesPtr_;

    // Face metric field pointers
    faceField<vector,MeshType>* faceCentersPtr_;
    faceField<vector,MeshType>* faceNormalsPtr_;
    faceField<scalar,MeshType>* faceAreasPtr_;
    faceField<vector,MeshType>* faceAreaNormalsPtr_;
    faceField<scalar,MeshType>* faceDeltasPtr_;
    faceField<scalar,MeshType>* faceWeightsCenterPtr_;
    faceField<scalar,MeshType>* faceWeightsNeighborPtr_;

    // Vertex metric field pointers
    meshField<vertexVector,MeshType>* vertexCentersPtr_;

    // Global cell number field pointer
    meshField<label,MeshType>* globalCellNumbersPtr_;

    // Global cell counts per level
    List<labelList> globalCellCounts_;

    // List of immersed boundaries (if any)
    FastPtrList<immersedBoundary<MeshType>> immersedBoundaries_;

    // Global immersed boundary mask pointer
    meshField<scalar,MeshType>* maskPtr_;

public:

    // On-demand Array of Structure (AoS) metrics
    class AoS
    {
        // Reference to the mesh metrics object
        const fvMeshMetrics<MeshType>& metrics_;

    public:

        // Construct from fvMeshMetrics
        AoS(const fvMeshMetrics<MeshType>& metrics)
        :
            metrics_(metrics)
        {}

        // AoS metrics
        tmp<meshField<faceVector,MeshType>> faceCenters() const;
        tmp<meshField<faceVector,MeshType>> faceNormals() const;
        tmp<meshField<faceScalar,MeshType>> faceAreas() const;
        tmp<meshField<faceVector,MeshType>> faceAreaNormals() const;
        tmp<meshField<faceScalar,MeshType>> faceDeltas() const;
        tmp<meshField<faceScalar,MeshType>> faceWeightsCenter() const;
        tmp<meshField<faceScalar,MeshType>> faceWeightsNeighbor() const;
        tmp<meshField<edgeVector,MeshType>> edgeCenters() const;
    };

private:

    // On-demand Array of Structure (AoS) metrics
    autoPtr<AoS> aosPtr_;

    // Functions to set cell metrics
    void setCellCenters();
    void setCellVolumes();

    // Functions to set face metrics
    void setFaceCenters();
    void setFaceNormals();
    void setFaceAreas();
    void setFaceAreaNormals();
    void setFaceDeltas();
    void setFaceWeights();

    // Function to set vertex metrics
    void setVertexCenters();

    // Function to set global cell numbers
    void setGlobalCellNumbers();

    // Function to set immersed boundaries
    void setImmersedBoundaries();
    void setImmersedBoundaryMask();

public:

    // Construct from fvMesh
    fvMeshMetrics(const fvMesh&);

    // Destructor
    ~fvMeshMetrics();

    // Const access

    const FastPtrList<immersedBoundary<MeshType>>& immersedBoundaries() const
    {
        return immersedBoundaries_;
    }

    inline const fvMesh& fvMsh() const
    {
        return fvMsh_;
    }

    inline const meshField<vector,MeshType>& cellCenters() const
    {
        if (!cellCentersPtr_)
            const_cast<fvMeshMetrics<MeshType>*>(this)->setCellCenters();

        return *cellCentersPtr_;
    }

    inline const meshField<scalar,MeshType>& cellVolumes() const
    {
        if (!cellVolumesPtr_)
            const_cast<fvMeshMetrics<MeshType>*>(this)->setCellVolumes();

        return *cellVolumesPtr_;
    }

    inline const meshField<scalar,MeshType>& inverseCellVolumes() const
    {
        if (!inverseCellVolumesPtr_)
            const_cast<fvMeshMetrics<MeshType>*>(this)->setCellVolumes();

        return *inverseCellVolumesPtr_;
    }

    inline const faceField<vector,MeshType>& faceCenters() const
    {
        if (!faceCentersPtr_)
            const_cast<fvMeshMetrics<MeshType>*>(this)->setFaceCenters();

        return *faceCentersPtr_;
    }

    inline const faceField<vector,MeshType>& faceNormals() const
    {
        if (!faceNormalsPtr_)
            const_cast<fvMeshMetrics<MeshType>*>(this)->setFaceNormals();

        return *faceNormalsPtr_;
    }

    inline const faceField<scalar,MeshType>& faceAreas() const
    {
        if (!faceAreasPtr_)
            const_cast<fvMeshMetrics<MeshType>*>(this)->setFaceAreas();

        return *faceAreasPtr_;
    }

    inline const faceField<vector,MeshType>& faceAreaNormals() const
    {
        if (!faceAreaNormalsPtr_)
            const_cast<fvMeshMetrics<MeshType>*>(this)->setFaceAreaNormals();

        return *faceAreaNormalsPtr_;
    }

    inline const faceField<scalar,MeshType>& faceDeltas() const
    {
        if (!faceDeltasPtr_)
            const_cast<fvMeshMetrics<MeshType>*>(this)->setFaceDeltas();

        return *faceDeltasPtr_;
    }

    inline const faceField<scalar,MeshType>& faceWeightsCenter() const
    {
        if (!faceWeightsCenterPtr_)
            const_cast<fvMeshMetrics<MeshType>*>(this)->setFaceWeights();

        return *faceWeightsCenterPtr_;
    }

    inline const faceField<scalar,MeshType>& faceWeightsNeighbor() const
    {
        if (!faceWeightsNeighborPtr_)
            const_cast<fvMeshMetrics<MeshType>*>(this)->setFaceWeights();

        return *faceWeightsNeighborPtr_;
    }

    inline const meshField<vertexVector,MeshType>& vertexCenters() const
    {
        if (!vertexCentersPtr_)
            const_cast<fvMeshMetrics<MeshType>*>(this)->setVertexCenters();

        return *vertexCentersPtr_;
    }

    inline const meshField<label,MeshType>& globalCellNumbers() const
    {
        if (!globalCellNumbersPtr_)
            const_cast<fvMeshMetrics<MeshType>*>(this)->setGlobalCellNumbers();

        return *globalCellNumbersPtr_;
    }

    inline const List<labelList>& globalCellCounts() const
    {
        if (!globalCellNumbersPtr_)
            const_cast<fvMeshMetrics<MeshType>*>(this)->setGlobalCellNumbers();

        return globalCellCounts_;
    }

    inline const meshField<scalar,MeshType>& mask() const
    {
        if (!maskPtr_)
            const_cast<fvMeshMetrics<MeshType>*>(this)
          ->setImmersedBoundaryMask();

        return *maskPtr_;
    }

    inline tmp<meshField<scalar,MeshType>> fluidMask() const
    {
        return 1.0 - mask();
    }

    inline const AoS& aos() const
    {
        return aosPtr_;
    }

    // Befriend AoS
    friend class AoS;
};

}

}

}

#endif
