#ifndef EigenLinearSystem_H
#define EigenLinearSystem_H

#include "EigenSolver.H"
#include "meshDirection.H"
#include "linearSystemAggregation.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
class EigenLinearSystem
:
    public linearSystemAggregation<SType,Type,MeshType>
{
    // Reference to the linear system
    const linearSystem<SType,Type,MeshType>& sys_;

    // Row and column-major Eigen matrices
    FastPtrList<EigenSolverBase::RowMajorMatrixType> rowMajorMatrices_;
    FastPtrList<EigenSolverBase::ColMajorMatrixType> colMajorMatrices_;

    // Level number
    const label l_;

public:

    // Construct from components
    EigenLinearSystem
    (
        const linearSystem<SType,Type,MeshType>& sys,
        const label l,
        const label nParts
    );

    // Copy constructor
    EigenLinearSystem(const EigenLinearSystem<SType,Type,MeshType>&);

    // Destrictor
    ~EigenLinearSystem();

    // Prepare the matrix
    void prepare(const label order, const label d);

    // Return reference to the row-major matrices

    inline EigenSolverBase::RowMajorMatrixType&
    rowMajorMatrix(const label d)
    {
        return rowMajorMatrices_[d];
    }

    inline const EigenSolverBase::RowMajorMatrixType&
    rowMajorMatrix(const label d) const
    {
        return rowMajorMatrices_[d];
    }
    // Return reference to the column-major matrices

    inline EigenSolverBase::ColMajorMatrixType& colMajorMatrix(const label d)
    {
        return colMajorMatrices_[d];
    }

    inline const EigenSolverBase::ColMajorMatrixType&
    colMajorMatrix(const label d) const
    {
        return colMajorMatrices_[d];
    }
};

}

}

}

#ifdef NoRepository
    #include "EigenLinearSystem.C"
#endif

#endif
