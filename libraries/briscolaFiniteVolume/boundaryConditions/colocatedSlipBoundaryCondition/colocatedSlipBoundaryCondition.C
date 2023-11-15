#include "colocatedSlipBoundaryCondition.H"

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
colocatedSlipBoundaryCondition<Type,MeshType>::colocatedSlipBoundaryCondition
(
    const meshField<Type,MeshType>& mshField,
    const partPatch& patch
)
:
    boundaryCondition<Type,MeshType>(mshField, patch)
{}

template<class Type, class MeshType>
colocatedSlipBoundaryCondition<Type,MeshType>::colocatedSlipBoundaryCondition
(
    const colocatedSlipBoundaryCondition<Type,MeshType>& bc
)
:
    boundaryCondition<Type,MeshType>(bc.mshField(), bc.patch())
{}

template<class Type, class MeshType>
colocatedSlipBoundaryCondition<Type,MeshType>::colocatedSlipBoundaryCondition
(
    const meshField<Type,MeshType>& field,
    const colocatedSlipBoundaryCondition<Type,MeshType>& bc
)
:
    boundaryCondition<Type,MeshType>(field, bc.patch())
{}

template<class Type, class MeshType>
void colocatedSlipBoundaryCondition<Type,MeshType>::initEvaluate(const label)
{}

template<class Type, class MeshType>
void colocatedSlipBoundaryCondition<Type,MeshType>::evaluate(const label l)
{
    meshLevel<Type,MeshType>& field = this->mshField()[l];

    const labelVector bo(this->boundaryOffset());

    forAll(field, d)
    {
        meshDirection<Type,MeshType>& fd = field[d];

        const labelVector S(this->S(l,d));
        const labelVector E(this->E(l,d));
        const block<Type> B(this->boundarySources(l,d));

        labelVector ijk;

        for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
        for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
        for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
        {
            fd(ijk+bo) = fd(ijk) + B(ijk-S);
        }
    }
}

template<class Type, class MeshType>
tmp<block<Type>> colocatedSlipBoundaryCondition<Type,MeshType>::boundarySources
(
    const label l,
    const label d
)
{
    const meshField<Type,MeshType>& field = this->mshField();
    const meshField<faceVector,MeshType>& fn = this->faceNormals();

    const labelVector bo(this->boundaryOffset());
    const label fb(faceNumber(bo));

    const labelVector S(this->S(l,d));
    const labelVector E(this->E(l,d));

    tmp<block<Type>> tR(new block<Type>(E-S));
    block<Type>& R = tR.ref();

    if (l == 0)
    {
        labelVector ijk;

        for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
        for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
        for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
        {
            R(ijk-S) =
                - 2.0 * (fn(l,d,ijk)[fb] & field(l,d,ijk)) * fn(l,d,ijk)[fb];
        }
    }
    else
    {
        R = Zero;
    }

    return tR;
}

}

}

}

