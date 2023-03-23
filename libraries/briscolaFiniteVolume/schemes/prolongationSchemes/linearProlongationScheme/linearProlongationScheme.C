#include "linearProlongationScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
void linearProlongationScheme<Type,MeshType>::setWeights()
{
    for (label l = 0; l < this->fvMsh().size()-1; l++)
    {
        const labelVector R(this->fvMsh()[l+1].R());

        for (label d = 0; d < MeshType::numberOfDirections; d++)
        {
            const meshDirection<vector,MeshType>& ccf =
                this->fvMsh().template metrics<MeshType>()
               .cellCenters()[l][d];

            const meshDirection<vector,MeshType>& ccc =
                this->fvMsh().template metrics<MeshType>()
               .cellCenters()[l+1][d];

            meshDirection<tensor,MeshType>& weights = weights_[l][d];

            weights = Zero;

            forAllCells(ccf, i, j, k)
            {
                scalar sum = 0.0;

                const label ox = (i % 2)*2 - 1;
                const label oy = (j % 2)*2 - 1;
                const label oz = (k % 2)*2 - 1;

                for (label a = 0; a < R.x(); a++)
                for (label b = 0; b < R.y(); b++)
                for (label c = 0; c < R.z(); c++)
                {
                    const label i0 = i/R.x() + a*ox;
                    const label j0 = j/R.y() + b*oy;
                    const label k0 = k/R.z() + c*oz;

                    const label i1 = i/R.x() + (1-a)*ox;
                    const label j1 = j/R.y() + (1-b)*oy;
                    const label k1 = k/R.z() + (1-c)*oz;

                    const vector dx = ccc(i1,j0,k0) - ccc(i0,j0,k0);
                    const vector dy = ccc(i0,j1,k0) - ccc(i0,j0,k0);
                    const vector dz = ccc(i0,j0,k1) - ccc(i0,j0,k0);

                    const scalar dxf = Foam::mag(dx);
                    const scalar dyf = Foam::mag(dy);
                    const scalar dzf = Foam::mag(dz);

                    const scalar wx =
                        R.x() > 1
                      ? Foam::mag((ccf(i,j,k)-ccc(i0,j0,k0)) & dx/dxf)/dxf
                      : 0;

                    const scalar wy =
                        R.y() > 1
                      ? Foam::mag((ccf(i,j,k)-ccc(i0,j0,k0)) & dy/dyf)/dyf
                      : 0;

                    const scalar wz =
                        R.z() > 1
                      ? Foam::mag((ccf(i,j,k)-ccc(i0,j0,k0)) & dz/dzf)/dzf
                      : 0;

                    const scalar w = (1.0-wx)*(1.0-wy)*(1.0-wz);

                    sum += w;

                    weights(i,j,k)[a*4+b*2+c] = w;
                }

                weights(i,j,k) /= sum;
            }
        }
    }
}

template<class Type, class MeshType>
template<template<class> class OpType>
void linearProlongationScheme<Type,MeshType>::prolong
(
    meshDirection<Type,MeshType>& fine,
    const meshDirection<Type,MeshType>& coarse,
    const OpType<Type>& bop
)
{
    const labelVector R(coarse.level().R());

    meshDirection<tensor,MeshType>& weights =
        weights_[fine.levelNum()][fine.directionNum()];

    forAllCells(fine, i, j, k)
    {
        const label il = i/R.x();
        const label jl = j/R.y();
        const label kl = k/R.z();

        const label ox = (i % 2)*2 - 1;
        const label oy = (j % 2)*2 - 1;
        const label oz = (k % 2)*2 - 1;

        const label iu = il + (R.x() == 1 ? 0 : ox);
        const label ju = jl + (R.y() == 1 ? 0 : oy);
        const label ku = kl + (R.z() == 1 ? 0 : oz);

        bop
        (
            fine(i,j,k),
            weights(i,j,k)[0] * coarse(il,jl,kl)
          + weights(i,j,k)[1] * coarse(il,jl,ku)
          + weights(i,j,k)[2] * coarse(il,ju,kl)
          + weights(i,j,k)[3] * coarse(il,ju,ku)
          + weights(i,j,k)[4] * coarse(iu,jl,kl)
          + weights(i,j,k)[5] * coarse(iu,jl,ku)
          + weights(i,j,k)[6] * coarse(iu,ju,kl)
          + weights(i,j,k)[7] * coarse(iu,ju,ku)
        );
    }
}

template<class Type, class MeshType>
linearProlongationScheme<Type,MeshType>::linearProlongationScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    prolongationScheme<Type,MeshType>(dict,fvMsh),
    weights_
    (
        "weights",
        this->fvMsh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    )
{
    setWeights();
}

template<class Type, class MeshType>
linearProlongationScheme<Type,MeshType>::linearProlongationScheme
(
    const fvMesh& fvMsh
)
:
    prolongationScheme<Type,MeshType>(dictionary(),fvMsh),
    weights_
    (
        "weights",
        this->fvMsh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    )
{
    setWeights();
}

template<class Type, class MeshType>
void linearProlongationScheme<Type,MeshType>::prolong
(
    meshDirection<Type,MeshType>& fine,
    const meshDirection<Type,MeshType>& coarse,
    const eqOp<Type>& bop
)
{
    this->prolong<eqOp>(fine, coarse, bop);
}

template<class Type, class MeshType>
void linearProlongationScheme<Type,MeshType>::prolong
(
    meshDirection<Type,MeshType>& fine,
    const meshDirection<Type,MeshType>& coarse,
    const plusEqOp<Type>& bop
)
{
    this->prolong<plusEqOp>(fine, coarse, bop);
}

template<class Type, class MeshType>
void linearProlongationScheme<Type,MeshType>::prolong
(
    meshDirection<Type,MeshType>& fine,
    const meshDirection<Type,MeshType>& coarse,
    const minusEqOp<Type>& bop
)
{
    this->prolong<minusEqOp>(fine, coarse, bop);
}

}

}

}
