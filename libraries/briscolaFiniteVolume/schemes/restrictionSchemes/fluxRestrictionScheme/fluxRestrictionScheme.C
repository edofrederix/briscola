#include "fluxRestrictionScheme.H"

#include "colocated.H"
#include "staggered.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class MeshType>
fluxRestrictionScheme<MeshType>::fluxRestrictionScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    restrictionScheme<faceScalar,MeshType>(dict,fvMsh)
{}

template<class MeshType>
fluxRestrictionScheme<MeshType>::fluxRestrictionScheme(const fvMesh& fvMsh)
:
    restrictionScheme<faceScalar,MeshType>(dictionary(),fvMsh)
{}

template<class MeshType>
void fluxRestrictionScheme<MeshType>::restrict
(
    meshDirection<faceScalar,MeshType>& coarse,
    const meshDirection<faceScalar,MeshType>& fine,
    const bool scale
)
{
    this->errorNoScaling(scale);

    const labelVector R(coarse.mshPart().R());

    coarse = Zero;

    // Sum fluxes from corresponding fine grid faces

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
                fine(fi, fj+jj, fk+kk).left();

            coarse(i,j,k).right() +=
                fine(fi+ii, fj+jj, fk+kk).right();
        }

        int jj = (R.y() == 2);
        for (int ii = 0; ii < R.x(); ii++)
        for (int kk = 0; kk < R.z(); kk++)
        {
            coarse(i,j,k).bottom() +=
                fine(fi+ii, fj, fk+kk).bottom();

            coarse(i,j,k).top() +=
                fine(fi+ii, fj+jj, fk+kk).top();
        }

        int kk = (R.z() == 2);
        for (int ii = 0; ii < R.x(); ii++)
        for (int jj = 0; jj < R.y(); jj++)
        {
            coarse(i,j,k).aft() +=
                fine(fi+ii, fj+jj, fk).aft();

            coarse(i,j,k).fore() +=
                fine(fi+ii, fj+jj, fk+kk).fore();
        }
    }
}

}

}

}
