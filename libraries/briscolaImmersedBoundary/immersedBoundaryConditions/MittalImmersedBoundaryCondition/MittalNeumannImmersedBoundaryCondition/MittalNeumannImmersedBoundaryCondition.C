#include "MittalNeumannImmersedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Constructor

template<class Type, class MeshType>
MittalNeumannImmersedBoundaryCondition<Type,MeshType>::
MittalNeumannImmersedBoundaryCondition
(
    const meshField<Type,MeshType>& mshField,
    const immersedBoundary<MeshType>& ib
)
:
    MittalImmersedBoundaryCondition<Type,MeshType>(mshField, ib),
    boundaryGradients_(this->dict().lookup("gradients"))
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

    meshDirection<Type,MeshType>& x = this->mshField_[l][d];

    const meshDirection<label,MeshType>& mask = this->forcingMask()[l][d];
    const meshDirection<vector,MeshType>& cc =
        this->fvMsh_.template metrics<MeshType>().cellCenters()[l][d];

    List<Type> data
    (
        move(this->exchanges_[l][d].dataFunc(this->mshField_))
    );

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
