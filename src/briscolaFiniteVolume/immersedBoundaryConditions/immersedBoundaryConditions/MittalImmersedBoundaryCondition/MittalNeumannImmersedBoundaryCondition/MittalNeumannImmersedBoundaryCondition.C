#include "MittalNeumannImmersedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Constructors

template<class Type, class MeshType>
MittalNeumannImmersedBoundaryCondition<Type,MeshType>::
MittalNeumannImmersedBoundaryCondition
(
    const meshField<Type,MeshType>& field,
    const immersedBoundary<MeshType>& ib
)
:
    MittalImmersedBoundaryCondition<Type,MeshType>(field, ib),
    boundaryGradients_(this->read("gradient"))
{}

template<class Type, class MeshType>
MittalNeumannImmersedBoundaryCondition<Type,MeshType>::
MittalNeumannImmersedBoundaryCondition
(
    const MittalNeumannImmersedBoundaryCondition<Type,MeshType>& ibc
)
:
    MittalImmersedBoundaryCondition<Type,MeshType>(ibc),
    boundaryGradients_(ibc.boundaryGradients_)
{}

template<class Type, class MeshType>
MittalNeumannImmersedBoundaryCondition<Type,MeshType>::
MittalNeumannImmersedBoundaryCondition
(
    const MittalNeumannImmersedBoundaryCondition<Type,MeshType>& ibc,
    const meshField<Type,MeshType>& field
)
:
    MittalImmersedBoundaryCondition<Type,MeshType>(ibc, field),
    boundaryGradients_(ibc.boundaryGradients_)
{}

// Destructor

template<class Type, class MeshType>
MittalNeumannImmersedBoundaryCondition<Type,MeshType>::
~MittalNeumannImmersedBoundaryCondition()
{}

template<class Type, class MeshType>
void MittalNeumannImmersedBoundaryCondition<Type,MeshType>::evaluate
(
    const label l,
    const label d
)
{
    const scalar omega = this->omega_;

    meshDirection<Type,MeshType>& x = this->field_[l][d];

    const meshDirection<label,MeshType>& mask = this->forcingMask()[l][d];
    const meshDirection<vector,MeshType>& cc =
        this->fvMsh_.template metrics<MeshType>().cellCenters()[l][d];

    List<Type> data(move(this->exchanges_[l][d](this->field_)));

    const vectorList& mps = this->exchanges_[l][d].interp().points();

    label c = 0;

    forAllCells(x,i,j,k)
    if (mask(i,j,k))
    {
        const scalar dist = mag(mps[c] - cc(i,j,k));

        x(i,j,k) =
            (1.0 - omega)*x(i,j,k)
          + omega*(dist*boundaryGradients_[d] + data[c]);

        c++;
    }
}

}

}

}
