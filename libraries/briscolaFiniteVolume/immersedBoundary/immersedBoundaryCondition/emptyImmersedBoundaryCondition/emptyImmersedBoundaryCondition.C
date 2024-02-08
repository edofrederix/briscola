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
    immersedBoundaryCondition<Type,MeshType>(mshField,ib)
{}

// Destructor

template<class Type, class MeshType>
emptyImmersedBoundaryCondition<Type,MeshType>::~emptyImmersedBoundaryCondition()
{}

}

}

}
