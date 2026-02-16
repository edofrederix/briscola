#include "emptyImmersedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Constructors

template<class Type, class MeshType>
emptyImmersedBoundaryCondition<Type,MeshType>::emptyImmersedBoundaryCondition
(
    const meshField<Type,MeshType>& field,
    const immersedBoundary<MeshType>& ib
)
:
    immersedBoundaryCondition<Type,MeshType>(field, ib, nullptr)
{}

template<class Type, class MeshType>
emptyImmersedBoundaryCondition<Type,MeshType>::
emptyImmersedBoundaryCondition
(
    const emptyImmersedBoundaryCondition<Type,MeshType>& ibc
)
:
    immersedBoundaryCondition<Type,MeshType>(ibc)
{}

template<class Type, class MeshType>
emptyImmersedBoundaryCondition<Type,MeshType>::
emptyImmersedBoundaryCondition
(
    const emptyImmersedBoundaryCondition<Type,MeshType>& ibc,
    const meshField<Type,MeshType>& field
)
:
    immersedBoundaryCondition<Type,MeshType>(ibc, field)
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
