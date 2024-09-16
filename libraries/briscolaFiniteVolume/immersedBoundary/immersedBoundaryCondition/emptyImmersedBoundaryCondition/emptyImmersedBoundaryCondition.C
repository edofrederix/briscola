#include "emptyImmersedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Constructor

template<class Type, class MeshType>
emptyImmersedBoundaryCondition<Type,MeshType>::emptyImmersedBoundaryCondition
(
    const meshField<Type,MeshType>& mshField,
    const immersedBoundary<MeshType>& ib
)
:
    immersedBoundaryCondition<Type,MeshType>(mshField,ib,ib.emptyField())
{}

// Destructor

template<class Type, class MeshType>
emptyImmersedBoundaryCondition<Type,MeshType>::~emptyImmersedBoundaryCondition()
{}

template<class Type, class MeshType>
void emptyImmersedBoundaryCondition<Type,MeshType>
::correctJacobiPoints
(
    meshLevel<Type,MeshType>& x
) const
{}

}

}

}
