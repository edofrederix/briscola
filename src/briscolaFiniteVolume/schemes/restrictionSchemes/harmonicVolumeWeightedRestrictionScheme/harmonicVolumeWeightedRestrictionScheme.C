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
    meshLevel<Type,MeshType>& coarse,
    const meshLevel<Type,MeshType>& fine,
    const bool scale
)
{
    this->errorNoScaling(scale);

    const labelVector R(coarse.lvl().R());

    const_cast<meshLevel<Type,MeshType>&>(fine).correctAggData();

    const meshLevel<scalar,MeshType>& cvc =
        this->fvMsh().template
        metrics<MeshType>().cellVolumes()[coarse.levelNum()];

    const meshLevel<scalar,MeshType>& cvf =
        this->fvMsh().template
        metrics<MeshType>().cellVolumes()[fine.levelNum()];

    forAll(coarse, d)
    {
        const vector shift(MeshType::shift[d]);

        const labelVector R2
        (
            shift.x() != 0.0 ? 1 : R.x(),
            shift.y() != 0.0 ? 1 : R.y(),
            shift.z() != 0.0 ? 1 : R.z()
        );

        forAllCells(coarse[d], i, j, k)
        {
            const label il = i*R.x();
            const label jl = j*R.y();
            const label kl = k*R.z();

            const label iu = i*R.x() + (R2.x() == 2);
            const label ju = j*R.y() + (R2.y() == 2);
            const label ku = k*R.z() + (R2.z() == 2);

            coarse(d,i,j,k) =
                cvf(d,il,jl,kl)/fine(d,il,jl,kl)
              + cvf(d,il,jl,ku)/fine(d,il,jl,ku)
              + cvf(d,il,ju,kl)/fine(d,il,ju,kl)
              + cvf(d,il,ju,ku)/fine(d,il,ju,ku)
              + cvf(d,iu,jl,kl)/fine(d,iu,jl,kl)
              + cvf(d,iu,jl,ku)/fine(d,iu,jl,ku)
              + cvf(d,iu,ju,kl)/fine(d,iu,ju,kl)
              + cvf(d,iu,ju,ku)/fine(d,iu,ju,ku);

            coarse(d,i,j,k) = cvc(d,i,j,k)/coarse(d,i,j,k);
        }
    }
}

}

}

}
