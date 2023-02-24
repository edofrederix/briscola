#include "DirichletBoundaryCondition.H"

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
DirichletBoundaryCondition<Type,MeshType>::DirichletBoundaryCondition
(
    const meshField<Type,MeshType>& mshField,
    const partPatch& patch
)
:
    boundaryCondition<Type,MeshType>(mshField, patch),
    boundaryValues_()
{
    const labelVector bo(this->boundaryOffset());
    List<Type> values(this->dict().lookup("values"));

    forAll(mshField, l)
    {
        forAll(mshField[l], d)
        {
            boundaryValues_.append
            (
                new block<Type>
                (
                    mshField[l][d].boundaryN(bo),
                    values[d]
                )
            );
        }
    }
}

template<class Type, class MeshType>
DirichletBoundaryCondition<Type,MeshType>::DirichletBoundaryCondition
(
    const DirichletBoundaryCondition<Type,MeshType>& bc
)
:
    boundaryCondition<Type,MeshType>(bc.mshField(), bc.patch()),
    boundaryValues_(bc.boundaryValues_)
{}

template<class Type, class MeshType>
DirichletBoundaryCondition<Type,MeshType>::DirichletBoundaryCondition
(
    const meshField<Type,MeshType>& field,
    const DirichletBoundaryCondition<Type,MeshType>& bc
)
:
    boundaryCondition<Type,MeshType>(field, bc.patch()),
    boundaryValues_(bc.boundaryValues_)
{}

template<class Type, class MeshType>
void DirichletBoundaryCondition<Type,MeshType>::initEvaluate(const label)
{}

template<class Type, class MeshType>
void DirichletBoundaryCondition<Type,MeshType>::evaluate
(
    const label l,
    const bool homogeneousBC
)
{
    meshLevel<Type,MeshType>& field = this->mshField()[l];

    const scalar H = homogeneousBC ? 0.0 : 1.0;

    const labelVector bo(this->boundaryOffset());

    forAll(field, d)
    {
        meshDirection<Type,MeshType>& fd = field[d];

        const block<Type>& val =
            boundaryValues_[l*field.size()+d];

        const labelVector S(fd.boundaryStart(bo));
        const labelVector E(fd.boundaryEnd(bo));

        // For shifted boundaries, the boundary values are directly set. For
        // non-shifted boundaries, set the ghost cell values appropriately.

        labelVector ijk;

        if (fd.shifted(bo))
        {
            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                fd(ijk+bo) = H*val(ijk-S);
            }
        }
        else
        {
            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                fd(ijk+bo) = 2.0*H*val(ijk-S) - fd(ijk);
            }
        }
    }
}

makeBoundaryConditionType(Dirichlet,label,colocated)
makeBoundaryConditionType(Dirichlet,label,staggered)

makeBoundaryConditionType(Dirichlet,scalar,colocated)
makeBoundaryConditionType(Dirichlet,scalar,staggered)

makeBoundaryConditionType(Dirichlet,faceScalar,colocated)
makeBoundaryConditionType(Dirichlet,faceScalar,staggered)

makeBoundaryConditionType(Dirichlet,vector,colocated)
makeBoundaryConditionType(Dirichlet,vector,staggered)

makeBoundaryConditionType(Dirichlet,faceVector,colocated)
makeBoundaryConditionType(Dirichlet,faceVector,staggered)

makeBoundaryConditionType(Dirichlet,tensor,colocated)
makeBoundaryConditionType(Dirichlet,tensor,staggered)

makeBoundaryConditionType(Dirichlet,sphericalTensor,colocated)
makeBoundaryConditionType(Dirichlet,sphericalTensor,staggered)

makeBoundaryConditionType(Dirichlet,symmTensor,colocated)
makeBoundaryConditionType(Dirichlet,symmTensor,staggered)

makeBoundaryConditionType(Dirichlet,diagTensor,colocated)
makeBoundaryConditionType(Dirichlet,diagTensor,staggered)

makeBoundaryConditionType(Dirichlet,stencil,colocated)
makeBoundaryConditionType(Dirichlet,stencil,staggered)

makeBoundaryConditionType(Dirichlet,diagStencil,colocated)
makeBoundaryConditionType(Dirichlet,diagStencil,staggered)

}

}

}

