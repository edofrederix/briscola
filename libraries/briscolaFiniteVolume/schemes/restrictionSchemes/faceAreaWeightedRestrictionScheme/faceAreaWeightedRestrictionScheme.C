#include "faceAreaWeightedRestrictionScheme.H"

#include "colocated.H"
#include "staggered.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
faceAreaWeightedRestrictionScheme<Type,MeshType>::
faceAreaWeightedRestrictionScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    restrictionScheme<Type,MeshType>(dict,fvMsh)
{}

template<class Type, class MeshType>
faceAreaWeightedRestrictionScheme<Type,MeshType>::
faceAreaWeightedRestrictionScheme
(
    const fvMesh& fvMsh
)
:
    restrictionScheme<Type,MeshType>(dictionary(),fvMsh)
{}

template<class Type, class MeshType>
void faceAreaWeightedRestrictionScheme<Type,MeshType>::restrict
(
    meshDirection<Type,MeshType>& coarse,
    const meshDirection<Type,MeshType>& fine,
    const bool scale
)
{
    this->errorNoScaling(scale);

    const labelVector R(coarse.level().R());

    coarse = Zero;

    const meshDirection<faceScalar,MeshType>& fac =
        this->fvMsh().template
        metrics<MeshType>().faceAreas()
        [coarse.levelNum()][coarse.directionNum()];

    const meshDirection<faceScalar,MeshType>& faf =
        this->fvMsh().template
        metrics<MeshType>().faceAreas()
        [fine.levelNum()][fine.directionNum()];

    // Face area weighted average of corresponding fine grid faces

    forAllCells(coarse, i, j, k)
    {
        const label fi = i*R.x();
        const label fj = j*R.y();
        const label fk = k*R.z();

        int ii = (R.x() == 2);
        for (int jj = 0; jj < R.y(); jj++)
        for (int kk = 0; kk < R.z(); kk++)
        {
            coarse(i,j,k).left() +=
                fine(fi, fj+jj, fk+kk).left()
              * faf(fi, fj+jj, fk+kk).left();

            coarse(i,j,k).right() +=
                fine(fi+ii, fj+jj, fk+kk).right()
              * faf(fi+ii, fj+jj, fk+kk).right();
        }

        coarse(i,j,k).left() /= fac(i,j,k).left();
        coarse(i,j,k).right() /= fac(i,j,k).right();

        int jj = (R.y() == 2);
        for (int ii = 0; ii < R.x(); ii++)
        for (int kk = 0; kk < R.z(); kk++)
        {
            coarse(i,j,k).bottom() +=
                fine(fi+ii, fj, fk+kk).bottom()
              * faf(fi+ii, fj, fk+kk).bottom();

            coarse(i,j,k).top() +=
                fine(fi+ii, fj+jj, fk+kk).top()
              * faf(fi+ii, fj+jj, fk+kk).top();
        }

        coarse(i,j,k).bottom() /= fac(i,j,k).bottom();
        coarse(i,j,k).top() /= fac(i,j,k).top();

        int kk = (R.z() == 2);
        for (int ii = 0; ii < R.x(); ii++)
        for (int jj = 0; jj < R.y(); jj++)
        {
            coarse(i,j,k).aft() +=
                fine(fi+ii, fj+jj, fk).aft()
              * faf(fi+ii, fj+jj, fk).aft();

            coarse(i,j,k).fore() +=
                fine(fi+ii, fj+jj, fk+kk).fore()
              * faf(fi+ii, fj+jj, fk+kk).fore();
        }

        coarse(i,j,k).aft() /= fac(i,j,k).aft();
        coarse(i,j,k).fore() /= fac(i,j,k).fore();
    }
}

}

}

}
