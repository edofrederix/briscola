#include "averageRestrictionScheme.H"

#include "colocated.H"
#include "staggered.H"
#include "meshFields.H"
#include "faceFields.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
void averageRestrictionScheme<Type,MeshType>::restrict
(
    meshLevel<Type,MeshType>& coarse,
    const meshLevel<Type,MeshType>& fine,
    const bool scale
)
{
    const labelVector R(coarse.lvl().R());

    const_cast<meshLevel<Type,MeshType>&>(fine).correctAggData();

    // In a shifted direction, coarse grid and fine grid cell centers have the
    // same coordinate in that direction. So in shifted directions, no
    // interpolation is required. In the limit of vertex fields (shifted in
    // three directions), this implies that coarse grid vertices coincide with
    // fine grid vertices, which they indeed do.

    if (scale)
    {
        const meshLevel<scalar,MeshType>& cvc =
            this->fvMsh().template
            metrics<MeshType>().cellVolumes()[coarse.levelNum()];

        const meshLevel<scalar,MeshType>& icvf =
            this->fvMsh().template
            metrics<MeshType>().inverseCellVolumes()[fine.levelNum()];

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
                    fine(d,il,jl,kl)*icvf(d,il,jl,kl)
                  + fine(d,il,jl,ku)*icvf(d,il,jl,ku)
                  + fine(d,il,ju,kl)*icvf(d,il,ju,kl)
                  + fine(d,il,ju,ku)*icvf(d,il,ju,ku)
                  + fine(d,iu,jl,kl)*icvf(d,iu,jl,kl)
                  + fine(d,iu,jl,ku)*icvf(d,iu,jl,ku)
                  + fine(d,iu,ju,kl)*icvf(d,iu,ju,kl)
                  + fine(d,iu,ju,ku)*icvf(d,iu,ju,ku);

                coarse(d,i,j,k) *= 0.125*cvc(d,i,j,k);
            }
        }
    }
    else
    {
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
                    fine(d,il,jl,kl)
                  + fine(d,il,jl,ku)
                  + fine(d,il,ju,kl)
                  + fine(d,il,ju,ku)
                  + fine(d,iu,jl,kl)
                  + fine(d,iu,jl,ku)
                  + fine(d,iu,ju,kl)
                  + fine(d,iu,ju,ku);

                coarse(d,i,j,k) *= 0.125;
            }
        }
    }
}

}

}

}
