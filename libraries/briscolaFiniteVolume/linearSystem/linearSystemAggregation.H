#ifndef linearSystemAggregation_H
#define linearSystemAggregation_H

#include "fvMesh.H"
#include "meshDirection.H"
#include "List.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
class linearSystemAggregation
{
    // Full stencil type
    typedef typename SType::fullStencilType FullSType;

    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Level and direction number
    const label l_;

    // Number of parts
    label nParts_;

    // Singularity bool for each direction
    List<bool> singular_;

    // This processor's part number
    label myPartNum_;

    // This processor's part master number
    label myPartMasterNum_;

    // Number of processor's per part
    label nProcsPerPart_;

    // Global row/column start and end indices for each direction
    labelList globalStarts_;
    labelList globalEnds_;

    // Column numbers for each processor's row in this part, for each direction.
    // Construct for a full stencil type and not a symmetric one.
    List<List<List<FixedList<label,FullSType::nComponents>>>> colNums_;

public:

    // Construct from components
    linearSystemAggregation
    (
        const fvMesh& fvMsh,
        const label l,
        const label nParts,
        const List<bool>& singular
    );

    // Copy constructor
    linearSystemAggregation
    (
        const linearSystemAggregation<SType,Type,MeshType>&
    );

    // Destrictor
    ~linearSystemAggregation();

    // Return labels
    inline label levelNum() const { return l_; }
    inline label nParts() const { return nParts_; }
    inline label myPartNum() const { return myPartNum_; }
    inline label myPartMasterNum() const { return myPartMasterNum_; }
    inline label nProcsPerPart() const { return nProcsPerPart_; }
    inline const labelList& globalStart() const { return globalStarts_; }
    inline const labelList& globalEnd() const { return globalEnds_; }

    // Return constraint bool
    inline const List<bool>& singular() { return singular_; }

    // Return reference to the column number data
    inline const List<List<List<FixedList<label,FullSType::nComponents>>>>&
    colNums() const
    {
        return colNums_;
    }

    // Return the aggregated matrix dimension of the given direction
    inline label size(const label d) const
    {
        return globalEnds_[d] - globalStarts_[d];
    }

    // Return true if this processor is master
    inline bool master() const
    {
        return Pstream::myProcNo() == myPartMasterNum_;
    }

    // Set the row coefficients of the aggregated matrix per processor
    void rowCoeffs
    (
        List<List<FullSType>>& rows,
        const meshDirection<SType,MeshType>& A
    ) const;

    // Set the right-hand side source
    void rhsSource
    (
        List<Type>& rhs,
        const meshDirection<SType,MeshType>& A,
        const meshDirection<Type,MeshType>& x,
        const meshDirection<Type,MeshType>& b
    ) const;

    // Return the aggregated matrix in compressed row format, optionally sorted
    void compressedRowFormat
    (
        scalarList& values,
        labelList& inners,
        labelList& outers,
        const meshDirection<SType,MeshType>& A,
        const bool sort = false
    ) const;

    // Correct the aggregation boundary conditions of the given level
    void correctBoundaryConditions(meshLevel<Type,MeshType>& x) const;
};

}

}

}

#ifdef NoRepository
    #include "linearSystemAggregation.C"
#endif

#endif
