#include "outflowBoundaryCondition.H"

#include "colocated.H"
#include "staggered.H"
#include "meshField.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Staggered

template<class Type>
void outflowBoundaryCondition<Type,staggered>::eliminateGhosts
(
    linearSystem<stencil,Type,staggered>& sys,
    const label l,
    const label d
)
{
    // Eliminate only for non-shifted boundaries

    if (faceNumber(this->offset())/2 != d)
        NeumannBoundaryCondition<Type,staggered>::eliminateGhosts(sys,l,d);
}

template<class Type>
void outflowBoundaryCondition<Type,staggered>::evaluate
(
    const label l,
    const label d
)
{
    const labelVector bo(this->offset());

    meshDirection<Type,staggered>& field = this->mshField_[l][d];

    const labelVector S(this->S(l,d));
    const labelVector E(this->E(l,d));

    labelVector ijk;
    for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
    for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
    for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
    {
        // Set velocity to zero if reverse flow
        if(Foam::mag(bo[d]) == 1 && field(ijk)*bo[d] < 0)
            field(ijk) = Zero;

        // Zero gradient
        field(ijk+bo) = field(ijk);
    }
}

}

}

}

