#include "volumeWeightedRestrictionScheme.H"

#include "colocated.H"
#include "staggered.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
void volumeWeightedRestrictionScheme<Type,MeshType>::restrict
(
    meshLevel<Type,MeshType>& coarse,
    const meshLevel<Type,MeshType>& fine,
    const bool scale
)
{
    this->errorNoScaling(scale);

    const labelVector R(coarse.lvl().R());

    const_cast<meshLevel<Type,MeshType>&>(fine).correctAggData();

    const meshLevel<scalar,MeshType>& icvc =
        this->fvMsh().template
        metrics<MeshType>().inverseCellVolumes()[coarse.levelNum()];

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
            coarse(d,i,j,k) = Zero;

            for (int ii = i*R.x(); ii < i*R.x()+R2.x(); ii++)
            for (int jj = j*R.y(); jj < j*R.y()+R2.y(); jj++)
            for (int kk = k*R.z(); kk < k*R.z()+R2.z(); kk++)
            {
                coarse(d,i,j,k) += fine(d,ii,jj,kk)*cvf(d,ii,jj,kk);
            }

            coarse(d,i,j,k) *= icvc(d,i,j,k);
        }
    }
}

}

}

}
