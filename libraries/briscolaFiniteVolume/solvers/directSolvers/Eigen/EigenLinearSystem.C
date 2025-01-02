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
    PtrList<EigenSolver::EigenMatrixType>(),
    linearSystemAggregation<SType,Type,MeshType>(sys, l, nParts),
    l_(l)
{
    this->update(sys);
}

template<class SType, class Type, class MeshType>
EigenLinearSystem<SType,Type,MeshType>::EigenLinearSystem
(
    const EigenLinearSystem<SType,Type,MeshType>& A
)
:
    PtrList<EigenSolver::EigenMatrixType>(A),
    linearSystemAggregation<SType,Type,MeshType>(A),
    l_(A.l_),
    diagonal_(A.diagonal_)
{}

template<class SType, class Type, class MeshType>
EigenLinearSystem<SType,Type,MeshType>::~EigenLinearSystem()
{}

template<class SType, class Type, class MeshType>
void EigenLinearSystem<SType,Type,MeshType>::update
(
    const linearSystem<SType,Type,MeshType>& sys
)
{
    linearSystemAggregation<SType,Type,MeshType>& lsa = *this;

    diagonal_ =
        const_cast<linearSystem<SType,Type,MeshType>&>(sys).diagonal();

    if (!this->size())
    {
        this->clear();
        this->resize(MeshType::numberOfDirections);

        forAll(*this, d)
            if (!diagonal_[d])
                this->set(d, new EigenSolver::EigenMatrixType());
    }

    forAll(*this, d)
    {
        if (!diagonal_[d])
        {
            scalarList values;
            labelList inners;
            labelList outers;

            lsa.compressedRowFormat
            (
                values,
                inners,
                outers,
                sys.A()[l_][d],
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

}
