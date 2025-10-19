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

    // Metrics
    faceField<vector,MeshType> faceCenters_;
    meshField<vertexVector,MeshType> vertexCenters_;
    faceField<vector,MeshType> faceNormals_;
    faceField<scalar,MeshType> faceAreas_;
    faceField<vector,MeshType> faceAreaNormals_;
    meshField<vector,MeshType> cellCenters_;
    meshField<scalar,MeshType> cellVolumes_;
    meshField<scalar,MeshType> inverseCellVolumes_;
    faceField<scalar,MeshType> faceDeltas_;
    faceField<scalar,MeshType> faceWeightsCenter_;
    faceField<scalar,MeshType> faceWeightsNeighbor_;

    // Global cell numbers
    meshField<label,MeshType> globalCellNumbers_;

    // Global cell counts per level
    List<labelList> globalCellCounts_;

    // List of immersed boundaries (if any)
    FastPtrList<immersedBoundary<MeshType>> immersedBoundaries_;

    // Global immersed boundary mask
    meshField<scalar,MeshType> mask_;

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
        tmp<meshField<edgeVector,MeshType>> edgeCenters() const;
        tmp<meshField<faceVector,MeshType>> faceNormals() const;
        tmp<meshField<faceScalar,MeshType>> faceAreas() const;
        tmp<meshField<faceVector,MeshType>> faceAreaNormals() const;
        tmp<meshField<faceScalar,MeshType>> faceDeltas() const;
        tmp<meshField<faceScalar,MeshType>> faceWeightsCenter() const;
        tmp<meshField<faceScalar,MeshType>> faceWeightsNeighbor() const;
    };

private:

    // On-demand Array of Structure (AoS) metrics
    autoPtr<AoS> aosPtr_;

    // Calculate metrics
    void calculateVertexCenters();
    void calculateFaceCenters();
    void calculateFaceAreasAndNormals();
    void calculateCellCenters();
    void calculateCellVolumes();
    void calculateFaceDeltas();
    void calculateFaceWeights();
    void setGlobalCellNumbers();
    void setImmersedBoundaries();
    void setImmersedBoundaryMask();

public:

    // Construct from fvMesh
    fvMeshMetrics(const fvMesh&);

    // Const access

    const FastPtrList<immersedBoundary<MeshType>>& immersedBoundaries() const
    {
        return immersedBoundaries_;
    }

    inline const fvMesh& fvMsh() const
    {
        return fvMsh_;
    }

    inline const faceField<vector,MeshType>& faceCenters() const
    {
        return faceCenters_;
    }

    inline const meshField<vertexVector,MeshType>& vertexCenters() const
    {
        return vertexCenters_;
    }

    inline const faceField<vector,MeshType>& faceNormals() const
    {
        return faceNormals_;
    }

    inline const faceField<scalar,MeshType>& faceAreas() const
    {
        return faceAreas_;
    }

    inline const faceField<vector,MeshType>& faceAreaNormals() const
    {
        return faceAreaNormals_;
    }

    inline const meshField<vector,MeshType>& cellCenters() const
    {
        return cellCenters_;
    }

    inline const meshField<scalar,MeshType>& cellVolumes() const
    {
        return cellVolumes_;
    }

    inline const meshField<scalar,MeshType>& inverseCellVolumes() const
    {
        return inverseCellVolumes_;
    }

    inline const faceField<scalar,MeshType>& faceDeltas() const
    {
        return faceDeltas_;
    }

    inline const faceField<scalar,MeshType>& faceWeightsCenter() const
    {
        return faceWeightsCenter_;
    }

    inline const faceField<scalar,MeshType>& faceWeightsNeighbor() const
    {
        return faceWeightsNeighbor_;
    }

    inline const meshField<label,MeshType>& globalCellNumbers() const
    {
        return globalCellNumbers_;
    }

    inline const List<labelList>& globalCellCounts() const
    {
        return globalCellCounts_;
    }

    inline const meshField<scalar,MeshType>& mask() const
    {
        return mask_;
    }

    inline tmp<meshField<scalar,MeshType>> fluidMask() const
    {
        return 1.0 - mask_;
    }

    inline tmp<meshField<scalar,MeshType>> solidMask() const
    {
        return mask_;
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
