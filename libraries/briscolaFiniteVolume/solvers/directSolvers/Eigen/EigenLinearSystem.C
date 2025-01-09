#include "EigenLinearSystem.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
EigenLinearSystem<SType,Type,MeshType>::EigenLinearSystem
(
    const linearSystem<SType,Type,MeshType>& sys,
    const label l,
    const label nParts
)
:
    linearSystemAggregation<SType,Type,MeshType>(sys, l, nParts),
    sys_(sys),
    rowMajorMatrices_(MeshType::numberOfDirections),
    colMajorMatrices_(MeshType::numberOfDirections),
    l_(l)
{}

template<class SType, class Type, class MeshType>
EigenLinearSystem<SType,Type,MeshType>::EigenLinearSystem
(
    const EigenLinearSystem<SType,Type,MeshType>& A
)
:
    linearSystemAggregation<SType,Type,MeshType>(A),
    sys_(A.sys_),
    rowMajorMatrices_(A.rowMajorMatrices_),
    colMajorMatrices_(A.colMajorMatrices_),
    l_(A.l_)
{}

template<class SType, class Type, class MeshType>
EigenLinearSystem<SType,Type,MeshType>::~EigenLinearSystem()
{}

template<class SType, class Type, class MeshType>
void EigenLinearSystem<SType,Type,MeshType>::prepare
(
    const label order,
    const label d
)
{
    const linearSystemAggregation<SType,Type,MeshType>& lsa = *this;

    scalarList values;
    labelList inners;
    labelList outers;

    lsa.compressedRowFormat
    (
        values,
        inners,
        outers,
        sys_,
        d,
        true
    );

    // Print inner indices for debugging

    // for (int j = 0; j < outers.size()-1; j++)
    //     for (int i = outers[j]; i < outers[j+1]; i++)
    //         Pout<< inners[i] << (i == outers[j+1]-1 ? nl : ' ');
    // Pout<< endl;

    if (lsa.master())
    {
        const label size = lsa.partSize(d);

        // Generate in row-major order

        EigenSolverBase::RowMajorMatrixType::Map M
        (
            size,
            size,
            values.size(),
            outers.begin(),
            inners.begin(),
            values.begin()
        );

        if (order == ::Eigen::RowMajor)
        {
            if (!rowMajorMatrices_.set(d))
                rowMajorMatrices_.set
                (
                    d,
                    new EigenSolverBase::RowMajorMatrixType()
                );

            rowMajorMatrices_[d] = M.eval();
            rowMajorMatrices_[d].makeCompressed();
        }
        else if (order == ::Eigen::ColMajor)
        {
            if (!colMajorMatrices_.set(d))
                colMajorMatrices_.set
                (
                    d,
                    new EigenSolverBase::ColMajorMatrixType()
                );

            colMajorMatrices_[d] = M.eval();
            colMajorMatrices_[d].makeCompressed();
        }
        else
        {
            FatalErrorInFunction
                << "Invalid order" << endl << abort(FatalError);
        }
    }
}

}

}

}
