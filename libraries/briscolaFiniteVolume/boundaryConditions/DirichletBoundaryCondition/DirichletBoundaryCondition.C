#include "DirichletBoundaryCondition.H"

#include "colocated.H"
#include "staggered.H"

#include "meshLevel.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
DirichletBoundaryCondition<Type,MeshType>::DirichletBoundaryCondition
(
    const meshField<Type,MeshType>& mshField,
    const partPatch& patch
)
:
    boundaryCondition<Type,MeshType>(mshField, patch),
    boundaryValues_(this->dict().lookup("values"))
{}

template<class Type, class MeshType>
DirichletBoundaryCondition<Type,MeshType>::DirichletBoundaryCondition
(
    const DirichletBoundaryCondition<Type,MeshType>& bc
)
:
    boundaryCondition<Type,MeshType>(bc.mshField(), bc.patch()),
    boundaryValues_(bc.boundaryValues_)
{}

template<class Type, class MeshType>
DirichletBoundaryCondition<Type,MeshType>::DirichletBoundaryCondition
(
    const meshField<Type,MeshType>& field,
    const DirichletBoundaryCondition<Type,MeshType>& bc
)
:
    boundaryCondition<Type,MeshType>(field, bc.patch()),
    boundaryValues_(bc.boundaryValues_)
{}

template<class Type, class MeshType>
void DirichletBoundaryCondition<Type,MeshType>::initEvaluate(const label)
{}

template<class Type, class MeshType>
void DirichletBoundaryCondition<Type,MeshType>::evaluate(const label l)
{
    meshLevel<Type,MeshType>& field = this->mshField()[l];

    const scalar H = l == 0;

    const labelVector bo(this->boundaryOffset());

    forAll(field, d)
    {
        meshDirection<Type,MeshType>& fd = field[d];

        const labelVector S(this->S(l,d));
        const labelVector E(this->E(l,d));

        labelVector ijk;

        if (MeshType::shifted(d,bo))
        {
            // Ghost value not needed because the internal value is constrained.
            // Set to the internal value.

            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                fd(ijk+bo) = H*boundaryValues_[d];
            }
        }
        else
        {
            // Eliminated so infer value from boundary source

            const block<Type> B(this->boundarySources(l,d));

            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                fd(ijk+bo) = H*B(ijk-S) - fd(ijk);
            }
        }
    }
}

template<class Type, class MeshType>
tmp<block<Type>> DirichletBoundaryCondition<Type,MeshType>::internalValue
(
    const label l,
    const label d
)
{
    const scalar H = l == 0;

    return tmp<block<Type>>
    (
        new block<Type>(this->N(l,d), H*boundaryValues_[d])
    );
}

template<class Type, class MeshType>
tmp<block<Type>> DirichletBoundaryCondition<Type,MeshType>::boundarySources
(
    const label l,
    const label d
)
{
    const scalar H = l == 0;

    return tmp<block<Type>>
    (
        new block<Type>(this->N(l,d), H*2.0*boundaryValues_[d])
    );
}

}

}

}

