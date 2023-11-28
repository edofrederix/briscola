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
    const boundary& b
)
:
    boundaryCondition<Type,MeshType>(mshField, b),
    boundaryGradients_(this->dict().lookup("gradients"))
{}

template<class Type, class MeshType>
NeumannBoundaryCondition<Type,MeshType>::NeumannBoundaryCondition
(
    const NeumannBoundaryCondition<Type,MeshType>& bc
)
:
    boundaryCondition<Type,MeshType>(bc.mshField(), bc.mshBoundary()),
    boundaryGradients_(bc.boundaryGradients_)
{}

template<class Type, class MeshType>
NeumannBoundaryCondition<Type,MeshType>::NeumannBoundaryCondition
(
    const meshField<Type,MeshType>& field,
    const NeumannBoundaryCondition<Type,MeshType>& bc
)
:
    boundaryCondition<Type,MeshType>(field, bc.mshBoundary()),
    boundaryGradients_(bc.boundaryGradients_)
{}

template<class Type, class MeshType>
void NeumannBoundaryCondition<Type,MeshType>::initEvaluate(const label)
{}

template<class Type, class MeshType>
void NeumannBoundaryCondition<Type,MeshType>::evaluate(const label l)
{
    meshLevel<Type,MeshType>& field = this->mshField()[l];

    const labelVector bo(this->offset());

    forAll(field, d)
    {
        meshDirection<Type,MeshType>& fd = field[d];

        const labelVector S(this->S(l,d));
        const labelVector E(this->E(l,d));
        const block<Type> B(this->boundarySources(l,d));

        labelVector ijk;

        // Eliminated so infer value from boundary source

        if (MeshType::shifted(d,bo))
        {
            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                fd(ijk+bo) = fd(ijk-bo) + B(ijk-S);
            }
        }
        else
        {
            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                fd(ijk+bo) = fd(ijk) + B(ijk-S);
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

    const labelVector bo(this->offset());

    if (MeshType::shifted(d,bo))
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
    const meshField<faceScalar,MeshType>& fd = this->faceDeltas();

    const labelVector bo(this->offset());
    const label fb(faceNumber(bo));
    const label fi(faceNumber(-bo));

    const labelVector S(this->S(l,d));
    const labelVector E(this->E(l,d));

    tmp<block<Type>> tR(new block<Type>(E-S));
    block<Type>& R = tR.ref();

    if (l == 0)
    {
        labelVector ijk;

        if (MeshType::shifted(d,bo))
        {
            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                const scalar twoDelta =
                    1.0/fd(l,d,ijk)[fb] + 1.0/fd(l,d,ijk)[fi];

                R(ijk-S) = boundaryGradients_[d]*twoDelta;
            }
        }
        else
        {
            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                R(ijk-S) = boundaryGradients_[d]/fd(l,d,ijk)[fb];
            }
        }
    }
    else
    {
        R = Zero;
    }

    return tR;
}

}

}

}

