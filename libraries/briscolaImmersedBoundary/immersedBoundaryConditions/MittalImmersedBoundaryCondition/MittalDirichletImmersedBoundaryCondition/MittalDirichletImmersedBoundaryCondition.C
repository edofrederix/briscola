#include "MittalDirichletImmersedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Constructor

template<class Type, class MeshType>
MittalDirichletImmersedBoundaryCondition<Type,MeshType>::
MittalDirichletImmersedBoundaryCondition
(
    const meshField<Type,MeshType>& mshField,
    const immersedBoundary<MeshType>& ib
)
:
    MittalImmersedBoundaryCondition<Type,MeshType>(mshField, ib),
    boundaryValues_(this->dict().lookup("values"))
{}

// Destructor

template<class Type, class MeshType>
MittalDirichletImmersedBoundaryCondition<Type,MeshType>::
~MittalDirichletImmersedBoundaryCondition()
{}

template<class Type, class MeshType>
void MittalDirichletImmersedBoundaryCondition<Type,MeshType>::evaluate
(
    const label l,
    const label d
)
{
    const scalar omega = this->omega_;

    meshDirection<Type,MeshType>& x = this->mshField_[l][d];
    const meshDirection<label,MeshType>& mask = this->forcingMask()[l][d];

    List<Type> data
    (
        move(this->exchanges_[l][d].dataFunc(this->mshField_))
    );

    label c = 0;

    forAllCells(x,i,j,k)
    if (mask(i,j,k))
    {
        x(i,j,k) =
            (1.0 - omega)*x(i,j,k)
          + omega*(2.0*boundaryValues_[d] - data[c]);

        c++;
    }
}

}

}

}
