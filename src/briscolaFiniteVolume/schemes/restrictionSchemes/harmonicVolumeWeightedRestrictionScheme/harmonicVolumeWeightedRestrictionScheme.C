#include "harmonicVolumeWeightedRestrictionScheme.H"

#include "colocated.H"
#include "staggered.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
void harmonicVolumeWeightedRestrictionScheme<Type,MeshType>::restrict
(
    meshDirection<Type,MeshType>& coarse,
    const meshDirection<Type,MeshType>& fine,
    const bool scale
)
{
    this->errorNoScaling(scale);

    const labelVector R(coarse.mshPart().R());
    const vector shift(MeshType::shift[coarse.directionNum()]);

    const labelVector R2
    (
        shift.x() != 0.0 ? 1 : R.x(),
        shift.y() != 0.0 ? 1 : R.y(),
        shift.z() != 0.0 ? 1 : R.z()
    );

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
            cvf(il,jl,kl)/fine(il,jl,kl)
          + cvf(il,jl,ku)/fine(il,jl,ku)
          + cvf(il,ju,kl)/fine(il,ju,kl)
          + cvf(il,ju,ku)/fine(il,ju,ku)
          + cvf(iu,jl,kl)/fine(iu,jl,kl)
          + cvf(iu,jl,ku)/fine(iu,jl,ku)
          + cvf(iu,ju,kl)/fine(iu,ju,kl)
          + cvf(iu,ju,ku)/fine(iu,ju,ku);

        coarse(i,j,k) = cvc(i,j,k)/coarse(i,j,k);
    }
}

}

}

}
