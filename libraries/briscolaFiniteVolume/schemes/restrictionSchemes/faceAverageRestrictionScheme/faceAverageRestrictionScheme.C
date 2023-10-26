#include "faceAverageRestrictionScheme.H"

#include "colocated.H"
#include "staggered.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
faceAverageRestrictionScheme<Type,MeshType>::faceAverageRestrictionScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    restrictionScheme<Type,MeshType>(dict,fvMsh)
{}

template<class Type, class MeshType>
faceAverageRestrictionScheme<Type,MeshType>::faceAverageRestrictionScheme
(
    const fvMesh& fvMsh
)
:
    restrictionScheme<Type,MeshType>(dictionary(),fvMsh)
{}

template<class Type, class MeshType>
void faceAverageRestrictionScheme<Type,MeshType>::restrict
(
    meshDirection<Type,MeshType>& coarse,
    const meshDirection<Type,MeshType>& fine,
    const bool scale
)
{
    this->errorNoScaling(scale);

    const labelVector R(coarse.level().R());

    coarse = Zero;

    // Average face values of corresponding fine grid faces

    forAllCells(coarse, i, j, k)
    {
        const label fi = i*R.x();
        const label fj = j*R.y();
        const label fk = k*R.z();

        int ii = (R.x() == 2);
        scalar nFaces = R.y()*R.z();
        for (int jj = 0; jj < R.y(); jj++)
        for (int kk = 0; kk < R.z(); kk++)
        {
            coarse(i,j,k).left() +=
                fine(fi, fj+jj, fk+kk).left()/nFaces;

            coarse(i,j,k).right() +=
                fine(fi+ii, fj+jj, fk+kk).right()/nFaces;
        }

        int jj = (R.y() == 2);
        nFaces = R.x()*R.z();
        for (int ii = 0; ii < R.x(); ii++)
        for (int kk = 0; kk < R.z(); kk++)
        {
            coarse(i,j,k).bottom() +=
                fine(fi+ii, fj, fk+kk).bottom()/nFaces;

            coarse(i,j,k).top() +=
                fine(fi+ii, fj+jj, fk+kk).top()/nFaces;
        }

        int kk = (R.z() == 2);
        nFaces = R.x()*R.y();
        for (int ii = 0; ii < R.x(); ii++)
        for (int jj = 0; jj < R.y(); jj++)
        {
            coarse(i,j,k).aft() +=
                fine(fi+ii, fj+jj, fk).aft()/nFaces;

            coarse(i,j,k).fore() +=
                fine(fi+ii, fj+jj, fk+kk).fore()/nFaces;
        }
    }
}

}

}

}
