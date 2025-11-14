#include "dummyBoundaryCondition.H"

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
dummyBoundaryCondition<Type,MeshType>::dummyBoundaryCondition
(
    const meshField<Type,MeshType>& mshField,
    const boundary& b
)
:
    boundaryCondition<Type,MeshType>(mshField, b)
{}

template<class Type, class MeshType>
dummyBoundaryCondition<Type,MeshType>::dummyBoundaryCondition
(
    const dummyBoundaryCondition<Type,MeshType>& bc
)
:
    boundaryCondition<Type,MeshType>(bc)
{}

template<class Type, class MeshType>
dummyBoundaryCondition<Type,MeshType>::dummyBoundaryCondition
(
    const meshField<Type,MeshType>& field,
    const dummyBoundaryCondition<Type,MeshType>& bc
)
:
    boundaryCondition<Type,MeshType>(field, bc)
{}

template<class Type, class MeshType>
void dummyBoundaryCondition<Type,MeshType>::evaluate
(
    const label l,
    const label d
)
{
    const labelVector bo(this->offset());

    meshDirection<Type,MeshType>& field = this->mshField_[l][d];

    const labelVector S(this->S(l,d));
    const labelVector E(this->E(l,d));

    labelVector ijk;

    // Set the ghost equal to the closest inner cell value

    for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
    for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
    for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
    {
        field(ijk+bo) = field(ijk);
    }
}

}

}

}

