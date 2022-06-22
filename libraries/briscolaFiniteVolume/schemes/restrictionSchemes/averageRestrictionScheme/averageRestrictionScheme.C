#include "averageRestrictionScheme.H"

#include "colocated.H"
#include "staggered.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
averageRestrictionScheme<Type,MeshType>::averageRestrictionScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    restrictionScheme<Type,MeshType>(dict,fvMsh)
{}

template<class Type, class MeshType>
averageRestrictionScheme<Type,MeshType>::averageRestrictionScheme(const fvMesh& fvMsh)
:
    restrictionScheme<Type,MeshType>(dictionary(),fvMsh)
{}

template<class Type, class MeshType>
void averageRestrictionScheme<Type,MeshType>::restrict
(
    meshDirection<Type,MeshType>& coarse,
    const meshDirection<Type,MeshType>& fine,
    const bool scale
)
{
    const labelVector R(coarse.level().R());
    const vector shift(MeshType::shift[coarse.directionNum()]);

    // In a shifted direction, coarse grid and fine grid cell centers have the
    // same coordinate in that direction. So in shifted directions, no
    // interpolation is required. In the limit of vertex fields (shifted in
    // three directions), this implies that coarse grid vertices coincide with
    // fine grid vertices, which they indeed do.

    const labelVector R2
    (
        shift.x() != 0.0 ? 1 : R.x(),
        shift.y() != 0.0 ? 1 : R.y(),
        shift.z() != 0.0 ? 1 : R.z()
    );

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

        forAllCells(coarse, i, j, k)
        {
            const label il = i*R.x();
            const label jl = j*R.y();
            const label kl = k*R.z();

            const label iu = i*R.x() + (R2.x() == 2);
            const label ju = j*R.y() + (R2.y() == 2);
            const label ku = k*R.z() + (R2.z() == 2);

            coarse(i,j,k) =
                fine(il,jl,kl)/cvf(il,jl,kl)
              + fine(il,jl,ku)/cvf(il,jl,ku)
              + fine(il,ju,kl)/cvf(il,ju,kl)
              + fine(il,ju,ku)/cvf(il,ju,ku)
              + fine(iu,jl,kl)/cvf(iu,jl,kl)
              + fine(iu,jl,ku)/cvf(iu,jl,ku)
              + fine(iu,ju,kl)/cvf(iu,ju,kl)
              + fine(iu,ju,ku)/cvf(iu,ju,ku);

            coarse(i,j,k) *= cvc(i,j,k)/8.0;
        }
    }
    else
    {
        forAllCells(coarse, i, j, k)
        {
            const label il = i*R.x();
            const label jl = j*R.y();
            const label kl = k*R.z();

            const label iu = i*R.x() + (R2.x() == 2);
            const label ju = j*R.y() + (R2.y() == 2);
            const label ku = k*R.z() + (R2.z() == 2);

            coarse(i,j,k) =
                fine(il,jl,kl)
              + fine(il,jl,ku)
              + fine(il,ju,kl)
              + fine(il,ju,ku)
              + fine(iu,jl,kl)
              + fine(iu,jl,ku)
              + fine(iu,ju,kl)
              + fine(iu,ju,ku);

            coarse(i,j,k) /= scalar(8.0);
        }
    }
}

}

}

}
