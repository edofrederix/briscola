#include "slipBoundaryCondition.H"

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
slipBoundaryCondition<Type,MeshType>::slipBoundaryCondition
(
    const meshField<Type,MeshType>& mshField,
    const boundary& b
)
:
    boundaryCondition<Type,MeshType>(mshField, b)
{}

template<class Type, class MeshType>
slipBoundaryCondition<Type,MeshType>::slipBoundaryCondition
(
    const slipBoundaryCondition<Type,MeshType>& bc
)
:
    boundaryCondition<Type,MeshType>(bc)
{}

template<class Type, class MeshType>
slipBoundaryCondition<Type,MeshType>::slipBoundaryCondition
(
    const meshField<Type,MeshType>& field,
    const slipBoundaryCondition<Type,MeshType>& bc
)
:
    boundaryCondition<Type,MeshType>(field, bc)
{}

template<class Type, class MeshType>
void slipBoundaryCondition<Type,MeshType>::evaluate
(
    const label l,
    const label d
)
{
    const labelVector bo(this->offset());
    const label faceNum(faceNumber(bo));

    meshDirection<Type,MeshType>& fd = this->mshField_[l][d];

    const meshDirection<faceVector,colocated>& fn = this->faceNormals()[l][d];

    const labelVector S(this->S(l,d));
    const labelVector E(this->E(l,d));

    labelVector ijk;
    for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
    for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
    for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
    {
        fd(ijk+bo) =
            fd(ijk) - 2.0*(fn(ijk)[faceNum] & fd(ijk)) * fn(ijk)[faceNum];
    }
}

}

}

}

