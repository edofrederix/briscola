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
    const meshLevel<SType,MeshType>& A,
    const label nParts,
    const List<bool>& singular
)
:
    PtrList<EigenSolver::EigenMatrixType>(),
    linearSystemAggregation<SType,Type,MeshType>
    (
        A.fvMsh(),
        A.levelNum(),
        nParts,
        singular
    )
{
    this->update(A, singular);
}

template<class SType, class Type, class MeshType>
EigenLinearSystem<SType,Type,MeshType>::EigenLinearSystem
(
    const EigenLinearSystem<SType,Type,MeshType>& A
)
:
    PtrList<EigenSolver::EigenMatrixType>(A),
    linearSystemAggregation<SType,Type,MeshType>(A)
{}

template<class SType, class Type, class MeshType>
EigenLinearSystem<SType,Type,MeshType>::~EigenLinearSystem()
{}

template<class SType, class Type, class MeshType>
void EigenLinearSystem<SType,Type,MeshType>::update
(
    const meshLevel<SType,MeshType>& A,
    const List<bool>& singular
)
{
    linearSystemAggregation<SType,Type,MeshType>& lsa = *this;

    if (!this->size())
    {
        this->clear();
        this->resize(MeshType::numberOfDirections);

        forAll(*this, d)
            this->set(d, new EigenSolver::EigenMatrixType());
    }

    forAll(*this, d)
    {
        scalarList values;
        labelList inners;
        labelList outers;

        lsa.compressedRowFormat
        (
            values,
            inners,
            outers,
            A[d],
            true
        );

        if (lsa.master())
        {
            const label size = lsa.size(d);

            this->operator[](d) =
                ::Eigen::SparseMatrix<double, ::Eigen::RowMajor>::Map
                (
                    size,
                    size,
                    values.size(),
                    outers.begin(),
                    inners.begin(),
                    values.begin()
                ).eval();

            this->operator[](d).makeCompressed();
        }
    }
}

}

}

}
