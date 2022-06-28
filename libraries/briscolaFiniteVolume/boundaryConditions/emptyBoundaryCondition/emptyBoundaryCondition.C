#include "emptyBoundaryCondition.H"

#include "colocated.H"
#include "staggered.H"

#include "meshLevel.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
emptyBoundaryCondition<Type,MeshType>::emptyBoundaryCondition
(
    const meshField<Type,MeshType>& mshField,
    const partPatch& patch
)
:
    boundaryCondition<Type,MeshType>(mshField, patch)
{
    // Empty boundaries can only sit normal to unit mesh size. Check.

    const partLevel& l = this->fvMsh_[0];
    const labelVector bo(this->boundaryOffset());

    for (label d = 0; d < 3; d++)
        if (Foam::mag(bo[d]) != 0 && l.N()[d] != 1)
            FatalErrorInFunction
                << "The empty boundary condition can only be applied "
                << "normal to mesh dimensions with only 1 cell thickness."
                << endl << exit(FatalError);

}

template<class Type, class MeshType>
emptyBoundaryCondition<Type,MeshType>::emptyBoundaryCondition
(
    const emptyBoundaryCondition<Type,MeshType>& bc
)
:
    boundaryCondition<Type,MeshType>(bc.mshField(), bc.patch())
{}

template<class Type, class MeshType>
emptyBoundaryCondition<Type,MeshType>::emptyBoundaryCondition
(
    const meshField<Type,MeshType>& field,
    const emptyBoundaryCondition<Type,MeshType>& bc
)
:
    boundaryCondition<Type,MeshType>(field, bc.patch())
{}

template<class Type, class MeshType>
void emptyBoundaryCondition<Type,MeshType>::initEvaluate(const label)
{}

template<class Type, class MeshType>
void emptyBoundaryCondition<Type,MeshType>::evaluate
(
    const label l,
    const bool
)
{
    meshLevel<Type,MeshType>& field = this->mshField()[l];

    const labelVector bo(this->boundaryOffset());

    forAll(field, d)
    {
        meshDirection<Type,MeshType>& fd = field[d];

        const labelVector S(fd.boundaryStart(bo));
        const labelVector E(fd.boundaryEnd(bo));

        // For shifted boundaries, the boundary values are non-eliminated and
        // set to zero. For non-shifted boundaries, apply homogeneous Neumann

        labelVector ijk;

        if (fd.shifted(bo))
        {
            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                fd(ijk) = Zero;
            }
        }
        else
        {
            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                fd(ijk+bo) = fd(ijk);
            }
        }
    }
}

makeBoundaryConditionType(empty,label,colocated)
makeBoundaryConditionType(empty,label,staggered)

makeBoundaryConditionType(empty,scalar,colocated)
makeBoundaryConditionType(empty,scalar,staggered)

makeBoundaryConditionType(empty,faceScalar,colocated)
makeBoundaryConditionType(empty,faceScalar,staggered)

makeBoundaryConditionType(empty,vector,colocated)
makeBoundaryConditionType(empty,vector,staggered)

makeBoundaryConditionType(empty,faceVector,colocated)
makeBoundaryConditionType(empty,faceVector,staggered)

makeBoundaryConditionType(empty,tensor,colocated)
makeBoundaryConditionType(empty,tensor,staggered)

makeBoundaryConditionType(empty,sphericalTensor,colocated)
makeBoundaryConditionType(empty,sphericalTensor,staggered)

makeBoundaryConditionType(empty,symmTensor,colocated)
makeBoundaryConditionType(empty,symmTensor,staggered)

makeBoundaryConditionType(empty,diagTensor,colocated)
makeBoundaryConditionType(empty,diagTensor,staggered)

makeBoundaryConditionType(empty,stencil,colocated)
makeBoundaryConditionType(empty,stencil,staggered)

makeBoundaryConditionType(empty,diagStencil,colocated)
makeBoundaryConditionType(empty,diagStencil,staggered)

}

}

}

