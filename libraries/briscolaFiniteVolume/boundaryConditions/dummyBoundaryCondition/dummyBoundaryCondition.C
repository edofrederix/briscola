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

}

}

}

