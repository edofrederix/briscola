#include "MittalDirichletImmersedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Constructors

template<class Type, class MeshType>
MittalDirichletImmersedBoundaryCondition<Type,MeshType>::
MittalDirichletImmersedBoundaryCondition
(
    const meshField<Type,MeshType>& field,
    const immersedBoundary<MeshType>& ib
)
:
    MittalImmersedBoundaryCondition<Type,MeshType>(field, ib),
    boundaryValues_(this->read("value"))
{}

template<class Type, class MeshType>
MittalDirichletImmersedBoundaryCondition<Type,MeshType>::
MittalDirichletImmersedBoundaryCondition
(
    const MittalDirichletImmersedBoundaryCondition<Type,MeshType>& ibc
)
:
    MittalImmersedBoundaryCondition<Type,MeshType>(ibc),
    boundaryValues_(ibc.boundaryValues_)
{}

template<class Type, class MeshType>
MittalDirichletImmersedBoundaryCondition<Type,MeshType>::
MittalDirichletImmersedBoundaryCondition
(
    const MittalDirichletImmersedBoundaryCondition<Type,MeshType>& ibc,
    const meshField<Type,MeshType>& field
)
:
    MittalImmersedBoundaryCondition<Type,MeshType>(ibc, field),
    boundaryValues_(ibc.boundaryValues_)
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

    meshDirection<Type,MeshType>& x = this->field_[l][d];
    const meshDirection<label,MeshType>& mask = this->forcingMask()[l][d];

    List<Type> data(this->exchanges_[l][d](this->field_));

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
