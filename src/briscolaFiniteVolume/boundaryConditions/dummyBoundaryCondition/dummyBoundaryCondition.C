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
    const meshLevel<Type,MeshType>& level,
    const boundary& b
)
:
    boundaryCondition<Type,MeshType>(level, b)
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
    const meshLevel<Type,MeshType>& level,
    const dummyBoundaryCondition<Type,MeshType>& bc
)
:
    boundaryCondition<Type,MeshType>(level, bc)
{}

template<class Type, class MeshType>
void dummyBoundaryCondition<Type,MeshType>::evaluate(const label d)
{
    const labelVector bo(this->offset());

    meshDirection<Type,MeshType>& field = this->level_[d];

    const labelVector S(this->S(d));
    const labelVector E(this->E(d));

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

