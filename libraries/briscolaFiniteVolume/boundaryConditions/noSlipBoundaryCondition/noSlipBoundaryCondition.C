#include "noSlipBoundaryCondition.H"

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
noSlipBoundaryCondition<Type,MeshType>::noSlipBoundaryCondition
(
    const meshField<Type,MeshType>& mshField,
    const partPatch& patch
)
:
    boundaryCondition<Type,MeshType>(mshField, patch)
{}

template<class Type, class MeshType>
noSlipBoundaryCondition<Type,MeshType>::noSlipBoundaryCondition
(
    const noSlipBoundaryCondition<Type,MeshType>& bc
)
:
    boundaryCondition<Type,MeshType>(bc.mshField(), bc.patch())
{}

template<class Type, class MeshType>
noSlipBoundaryCondition<Type,MeshType>::noSlipBoundaryCondition
(
    const meshField<Type,MeshType>& field,
    const noSlipBoundaryCondition<Type,MeshType>& bc
)
:
    boundaryCondition<Type,MeshType>(field, bc.patch())
{}

template<class Type, class MeshType>
void noSlipBoundaryCondition<Type,MeshType>::initEvaluate(const label)
{}

template<class Type, class MeshType>
void noSlipBoundaryCondition<Type,MeshType>::evaluate(const label l)
{
    meshLevel<Type,MeshType>& field = this->mshField()[l];

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
                fd(ijk+bo) = Zero;
            }
        }
        else
        {
            // Eliminated so infer value from boundary source

            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                fd(ijk+bo) = - fd(ijk);
            }
        }
    }
}

template<class Type, class MeshType>
tmp<block<Type>> noSlipBoundaryCondition<Type,MeshType>::internalValue
(
    const label l,
    const label d
)
{
    return tmp<block<Type>>
    (
        new block<Type>(this->N(l,d), Zero)
    );
}

template<class Type, class MeshType>
tmp<block<Type>> noSlipBoundaryCondition<Type,MeshType>::boundarySources
(
    const label l,
    const label d
)
{
    return tmp<block<Type>>
    (
        new block<Type>(this->N(l,d), Zero)
    );
}

}

}

}

