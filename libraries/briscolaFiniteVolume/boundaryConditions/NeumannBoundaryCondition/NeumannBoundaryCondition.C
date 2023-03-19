#include "NeumannBoundaryCondition.H"

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
NeumannBoundaryCondition<Type,MeshType>::NeumannBoundaryCondition
(
    const meshField<Type,MeshType>& mshField,
    const partPatch& patch
)
:
    boundaryCondition<Type,MeshType>(mshField, patch),
    boundaryGradients_()
{
    const labelVector bo(this->boundaryOffset());
    List<Type> gradients(this->dict().lookup("gradients"));

    forAll(mshField, l)
    {
        forAll(mshField[l], d)
        {
            boundaryGradients_.append
            (
                new block<Type>
                (
                    mshField[l][d].boundaryN(bo),
                    gradients[d]
                )
            );
        }
    }
}

template<class Type, class MeshType>
NeumannBoundaryCondition<Type,MeshType>::NeumannBoundaryCondition
(
    const NeumannBoundaryCondition<Type,MeshType>& bc
)
:
    boundaryCondition<Type,MeshType>(bc.mshField(), bc.patch()),
    boundaryGradients_(bc.boundaryGradients_)
{}

template<class Type, class MeshType>
NeumannBoundaryCondition<Type,MeshType>::NeumannBoundaryCondition
(
    const meshField<Type,MeshType>& field,
    const NeumannBoundaryCondition<Type,MeshType>& bc
)
:
    boundaryCondition<Type,MeshType>(field, bc.patch()),
    boundaryGradients_(bc.boundaryGradients_)
{}

template<class Type, class MeshType>
void NeumannBoundaryCondition<Type,MeshType>::initEvaluate(const label)
{}

template<class Type, class MeshType>
void NeumannBoundaryCondition<Type,MeshType>::evaluate
(
    const label l,
    const bool homogeneousBC
)
{
    meshLevel<Type,MeshType>& field = this->mshField()[l];
    const meshLevel<faceScalar,MeshType>& fdl = this->faceDeltas()[l];

    const scalar H = homogeneousBC ? 0.0 : 1.0;

    const labelVector bo(this->boundaryOffset());
    const label fb(faceNumber(bo));
    const label fi(faceNumber(-bo));

    forAll(field, d)
    {
        meshDirection<Type,MeshType>& fd = field[d];
        const meshDirection<faceScalar,MeshType>& fdld = fdl[d];

        const block<Type>& grad =
            boundaryGradients_[l*field.size()+d];

        const labelVector S(fd.activeBoundaryStart(bo));
        const labelVector E(fd.activeBoundaryEnd(bo));

        labelVector ijk;

        if (fd.shifted(bo))
        {
            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                const scalar delta =
                    1.0/fdld(ijk)[fb] + 1.0/fdld(ijk)[fi];

                fd(ijk+bo) =
                    fd(ijk-bo) + H*grad(ijk-S)*delta;
            }
        }
        else
        {
            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                fd(ijk+bo) =
                    fd(ijk) + H*grad(ijk-S)/fdld(ijk)[fb];
            }
        }
    }
}

makeBoundaryConditionTypes(Neumann)


}

}

}

