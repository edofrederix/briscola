#include "penalizationDirichletImmersedBoundaryCondition.H"
#include "immersedBoundary.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Constructor

template<class Type, class MeshType>
penalizationDirichletImmersedBoundaryCondition<Type,MeshType>::
penalizationDirichletImmersedBoundaryCondition
(
    const meshField<Type,MeshType>& mshField,
    const immersedBoundary<MeshType>& ib
)
:
    immersedBoundaryCondition<Type,MeshType>(mshField, ib, &ib.mask()),
    boundaryValues_(this->read("value"))
{}

// Destructor

template<class Type, class MeshType>
penalizationDirichletImmersedBoundaryCondition<Type,MeshType>::
~penalizationDirichletImmersedBoundaryCondition()
{}

template<class Type, class MeshType>
void penalizationDirichletImmersedBoundaryCondition<Type,MeshType>::evaluate
(
    const label l,
    const label d
)
{
    const scalar omega = this->omega_;

    meshDirection<Type,MeshType>& x = this->mshField_[l][d];
    const meshDirection<label,MeshType>& mask = this->forcingMask()[l][d];

    forAllCells(x,i,j,k)
        if (mask(i,j,k))
            x(i,j,k) = (1.0 - omega)*x(i,j,k) + omega*boundaryValues_[d];
}

}

}

}
