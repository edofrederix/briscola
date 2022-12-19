#include "linearRestrictionScheme.H"

#include "colocated.H"
#include "staggered.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
void linearRestrictionScheme<Type,MeshType>::setWeights()
{
    for (label l = 0; l < this->fvMsh().size()-1; l++)
    {
        const labelVector R(this->fvMsh()[l+1].R());

        for (label d = 0; d < MeshType::numberOfDirections; d++)
        {
            const vector shift(MeshType::shift[d]);

            const labelVector R2
            (
                shift.x() != 0.0 ? 1 : R.x(),
                shift.y() != 0.0 ? 1 : R.y(),
                shift.z() != 0.0 ? 1 : R.z()
            );

            const meshDirection<vector,MeshType>& ccf =
                this->fvMsh().template metrics<MeshType>()
               .cellCenters()[l][d];

            const meshDirection<vector,MeshType>& ccc =
                this->fvMsh().template metrics<MeshType>()
               .cellCenters()[l+1][d];

            meshDirection<tensor,MeshType>& weights =
                weights_[l+1][d];

            weights = Zero;

            forAllBlock(ccc, i, j, k)
            {
                scalar sum = 0.0;

                for (label a = 0; a < R2.x(); a++)
                for (label b = 0; b < R2.y(); b++)
                for (label c = 0; c < R2.z(); c++)
                {
                    const label i0 = i*R.x() + a;
                    const label j0 = j*R.y() + b;
                    const label k0 = k*R.z() + c;

                    const label i1 = i*R.x() + (1-a);
                    const label j1 = j*R.y() + (1-b);
                    const label k1 = k*R.z() + (1-c);

                    const vector dx = ccf(i1,j0,k0) - ccf(i0,j0,k0);
                    const vector dy = ccf(i0,j1,k0) - ccf(i0,j0,k0);
                    const vector dz = ccf(i0,j0,k1) - ccf(i0,j0,k0);

                    const scalar dxf = Foam::mag(dx);
                    const scalar dyf = Foam::mag(dy);
                    const scalar dzf = Foam::mag(dz);

                    const scalar wx =
                        R.x() > 1 ? Foam::mag((ccc(i,j,k)-ccf(i0,j0,k0)) & dx/dxf)/dxf : 0;

                    const scalar wy =
                        R.y() > 1 ? Foam::mag((ccc(i,j,k)-ccf(i0,j0,k0)) & dy/dyf)/dyf : 0;

                    const scalar wz =
                        R.z() > 1 ? Foam::mag((ccc(i,j,k)-ccf(i0,j0,k0)) & dz/dzf)/dzf : 0;

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
linearRestrictionScheme<Type,MeshType>::linearRestrictionScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    restrictionScheme<Type,MeshType>(dict,fvMsh),
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
linearRestrictionScheme<Type,MeshType>::linearRestrictionScheme(const fvMesh& fvMsh)
:
    restrictionScheme<Type,MeshType>(dictionary(),fvMsh),
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
void linearRestrictionScheme<Type,MeshType>::restrict
(
    meshDirection<Type,MeshType>& coarse,
    const meshDirection<Type,MeshType>& fine,
    const bool scale
)
{
    const labelVector R(coarse.level().R());
    const vector shift(MeshType::shift[coarse.directionNum()]);

    const labelVector R2
    (
        shift.x() != 0.0 ? 1 : R.x(),
        shift.y() != 0.0 ? 1 : R.y(),
        shift.z() != 0.0 ? 1 : R.z()
    );

    meshDirection<tensor,MeshType>& weights =
        weights_[coarse.levelNum()][coarse.directionNum()];

    if (scale)
    {
        const meshDirection<scalar,MeshType>& cvc =
            this->fvMsh().template
            metrics<MeshType>().cellVolumes()
            [coarse.levelNum()][coarse.directionNum()];

        const meshDirection<scalar,MeshType>& cvf =
            this->fvMsh().template
            metrics<MeshType>().cellVolumes()
            [fine.levelNum()][fine.directionNum()];

        forAllBlock(coarse, i, j, k)
        {
            const label il = i*R.x();
            const label jl = j*R.y();
            const label kl = k*R.z();

            const label iu = i*R.x() + (R2.x() == 2);
            const label ju = j*R.y() + (R2.y() == 2);
            const label ku = k*R.z() + (R2.z() == 2);

            coarse(i,j,k) =
                fine(il,jl,kl)/cvf(il,jl,kl)*weights(i,j,k)[0]
              + fine(il,jl,ku)/cvf(il,jl,ku)*weights(i,j,k)[1]
              + fine(il,ju,kl)/cvf(il,ju,kl)*weights(i,j,k)[2]
              + fine(il,ju,ku)/cvf(il,ju,ku)*weights(i,j,k)[3]
              + fine(iu,jl,kl)/cvf(iu,jl,kl)*weights(i,j,k)[4]
              + fine(iu,jl,ku)/cvf(iu,jl,ku)*weights(i,j,k)[5]
              + fine(iu,ju,kl)/cvf(iu,ju,kl)*weights(i,j,k)[6]
              + fine(iu,ju,ku)/cvf(iu,ju,ku)*weights(i,j,k)[7];

            coarse(i,j,k) *= cvc(i,j,k);
        }
    }
    else
    {
        forAllBlock(coarse, i, j, k)
        {
            const label il = i*R.x();
            const label jl = j*R.y();
            const label kl = k*R.z();

            const label iu = i*R.x() + (R2.x() == 2);
            const label ju = j*R.y() + (R2.y() == 2);
            const label ku = k*R.z() + (R2.z() == 2);

            coarse(i,j,k) =
                fine(il,jl,kl)*weights(i,j,k)[0]
              + fine(il,jl,ku)*weights(i,j,k)[1]
              + fine(il,ju,kl)*weights(i,j,k)[2]
              + fine(il,ju,ku)*weights(i,j,k)[3]
              + fine(iu,jl,kl)*weights(i,j,k)[4]
              + fine(iu,jl,ku)*weights(i,j,k)[5]
              + fine(iu,ju,kl)*weights(i,j,k)[6]
              + fine(iu,ju,ku)*weights(i,j,k)[7];
        }
    }
}

}

}

}
