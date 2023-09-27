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
    boundaryValues_()
{
    const labelVector bo(this->boundaryOffset());
    List<Type> values(this->dict().lookup("values"));

    forAll(mshField, l)
    {
        forAll(mshField[l], d)
        {
            boundaryValues_.append
            (
                new block<Type>
                (
                    mshField.boundaryN(l,d,bo),
                    values[d]
                )
            );
        }
    }
}

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
void DirichletBoundaryCondition<Type,MeshType>::evaluate
(
    const label l,
    const bool homogeneous
)
{
    meshLevel<Type,MeshType>& field = this->mshField()[l];

    const scalar H = homogeneous ? 0.0 : 1.0;

    const labelVector bo(this->boundaryOffset());

    forAll(field, d)
    {
        meshDirection<Type,MeshType>& fd = field[d];

        const labelVector S(fd.boundaryStart(bo));
        const labelVector E(fd.boundaryEnd(bo));

        const block<Type>& val =
            boundaryValues_[l*MeshType::numberOfDirections + d];

        labelVector ijk;

        if (fd.shifted(bo))
        {
            // Ghost value not needed because the internal value is constrained.
            // Set to the internal value.

            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                fd(ijk+bo) = H*val(ijk-S);
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
    tmp<block<Type>> tv
    (
        new block<Type>(boundaryValues_[l*MeshType::numberOfDirections + d])
    );

    return tv;
}

template<class Type, class MeshType>
tmp<block<Type>> DirichletBoundaryCondition<Type,MeshType>::boundarySources
(
    const label l,
    const label d
)
{
    return 2.0*boundaryValues_[l*MeshType::numberOfDirections + d];
}

}

}

}

