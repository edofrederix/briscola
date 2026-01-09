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
    const meshField<Type,MeshType>& field,
    const immersedBoundary<MeshType>& ib
)
:
    immersedBoundaryCondition<Type,MeshType>(field, ib, nullptr)
{}

// Destructor

template<class Type, class MeshType>
emptyImmersedBoundaryCondition<Type,MeshType>::~emptyImmersedBoundaryCondition()
{}

template<class Type, class MeshType>
void emptyImmersedBoundaryCondition<Type,MeshType>::evaluate
(
    const label,
    const label
)
{}

}

}

}
