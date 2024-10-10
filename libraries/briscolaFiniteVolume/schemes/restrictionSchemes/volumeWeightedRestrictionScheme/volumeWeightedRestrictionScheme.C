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
volumeWeightedRestrictionScheme<Type,MeshType>::volumeWeightedRestrictionScheme
(
    const fvMesh& fvMsh,
    Istream& is
)
:
    restrictionScheme<Type,MeshType>(fvMsh, is)
{}

template<class Type, class MeshType>
void volumeWeightedRestrictionScheme<Type,MeshType>::restrict
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
        coarse(i,j,k) = Zero;

        for (int ii = i*R.x(); ii < i*R.x()+R2.x(); ii++)
        for (int jj = j*R.y(); jj < j*R.y()+R2.y(); jj++)
        for (int kk = k*R.z(); kk < k*R.z()+R2.z(); kk++)
        {
            coarse(i,j,k) += fine(ii,jj,kk)*cvf(ii,jj,kk);
        }

        coarse(i,j,k) /= cvc(i,j,k);
    }
}

}

}

}
