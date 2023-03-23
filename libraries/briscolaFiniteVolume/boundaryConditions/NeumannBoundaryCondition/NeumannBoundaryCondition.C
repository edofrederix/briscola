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
    const bool homogeneous
)
{
    meshLevel<Type,MeshType>& field = this->mshField()[l];

    const scalar H = homogeneous ? 0.0 : 1.0;

    const labelVector bo(this->boundaryOffset());

    forAll(field, d)
    {
        meshDirection<Type,MeshType>& fd = field[d];

        const labelVector S(fd.boundaryStart(bo));
        const labelVector E(fd.boundaryEnd(bo));
        const block<Type> B(this->boundarySources(l,d));

        labelVector ijk;

        // Eliminated so infer value from boundary source

        if (fd.shifted(bo))
        {
            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                fd(ijk+bo) = fd(ijk-bo) + H*B(ijk-S);
            }
        }
        else
        {
            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                fd(ijk+bo) = fd(ijk) + H*B(ijk-S);
            }
        }
    }
}

template<class Type, class MeshType>
stencil NeumannBoundaryCondition<Type,MeshType>::boundaryCoeff
(
    const label l,
    const label d
)
{
    // For shifted boundaries, add the boundary coefficient to the opposite
    // coefficient. For non-shifted boundaries, add the boundary coefficient to
    // the center coefficient.

    const labelVector bo(this->boundaryOffset());

    if (this->mshField_[l][d].shifted(bo))
    {
        stencil S(Zero);
        S[faceNumber(-bo)+1] = 1.0;

        return S;
    }
    else
    {
        return pTraits<diagStencil>::one;
    }
}

template<class Type, class MeshType>
tmp<block<Type>> NeumannBoundaryCondition<Type,MeshType>::boundarySources
(
    const label l,
    const label d
)
{
    const meshLevel<Type,MeshType>& fl = this->mshField_[l];
    const meshDirection<Type,MeshType>& fld = fl[d];
    const meshDirection<faceScalar,MeshType>& fdld = this->faceDeltas()[l][d];

    const block<Type>& grad =
        boundaryGradients_[l*fl.size()+d];

    const labelVector bo(this->boundaryOffset());
    const label fb(faceNumber(bo));
    const label fi(faceNumber(-bo));

    const labelVector S(fld.boundaryStart(bo));
    const labelVector E(fld.boundaryEnd(bo));

    tmp<block<Type>> tR(new block<Type>(E-S));
    block<Type>& R = tR.ref();

    labelVector ijk;

    if (fld.shifted(bo))
    {
        for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
        for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
        for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
        {
            const scalar delta =
                1.0/fdld(ijk)[fb] + 1.0/fdld(ijk)[fi];

            R(ijk-S) = grad(ijk-S)*delta;
        }
    }
    else
    {
        for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
        for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
        for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
        {
            R(ijk-S) = grad(ijk-S)/fdld(ijk)[fb];
        }
    }

    return tR;
}

}

}

}

