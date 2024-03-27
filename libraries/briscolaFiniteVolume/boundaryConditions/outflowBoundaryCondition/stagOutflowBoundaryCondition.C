#include "stagOutflowBoundaryCondition.H"

#include "colocated.H"
#include "staggered.H"
#include "meshLevel.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type>
stagOutflowBoundaryCondition<Type>::stagOutflowBoundaryCondition
(
    const meshField<Type,staggered>& mshField,
    const boundary& b
)
:
    outflowBoundaryCondition<Type,staggered>(mshField, b)
{}

template<class Type>
stagOutflowBoundaryCondition<Type>::stagOutflowBoundaryCondition
(
    const stagOutflowBoundaryCondition<Type>& bc
)
:
    outflowBoundaryCondition<Type,staggered>(bc.mshField(), bc.mshBoundary())
{}

template<class Type>
stagOutflowBoundaryCondition<Type>::stagOutflowBoundaryCondition
(
    const meshField<Type,staggered>& field,
    const stagOutflowBoundaryCondition<Type>& bc
)
:
    outflowBoundaryCondition<Type,staggered>(field, bc.mshBoundary())
{}

template<class Type>
void stagOutflowBoundaryCondition<Type>::eliminateGhosts
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
void stagOutflowBoundaryCondition<Type>::evaluate
(
    const label l,
    const label d
)
{
    const labelVector bo(this->offset());

    meshDirection<Type,staggered>& fd = this->mshField_[l][d];

    const labelVector S(this->S(l,d));
    const labelVector E(this->E(l,d));

    labelVector ijk;
    for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
    for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
    for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
    {
        // Set velocity to zero if reverse flow
        if
        (
               (Foam::mag(bo[d]) == 1)
            && (fd(ijk)*bo[d] < 0)
        )
        {
            fd(ijk) = Zero;
        }

        // Zero gradient
        fd(ijk+bo) = fd(ijk);
    }
}

}

}

}

