#include "stagFaceAverageRestrictionScheme.H"

#include "staggered.H"
#include "staggered.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type>
stagFaceAverageRestrictionScheme<Type>::stagFaceAverageRestrictionScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    restrictionScheme<Type,staggered>(dict,fvMsh)
{}

template<class Type>
stagFaceAverageRestrictionScheme<Type>::stagFaceAverageRestrictionScheme
(
    const fvMesh& fvMsh
)
:
    restrictionScheme<Type,staggered>(dictionary(),fvMsh)
{}

template<class Type>
void stagFaceAverageRestrictionScheme<Type>::restrict
(
    meshDirection<Type,staggered>& coarse,
    const meshDirection<Type,staggered>& fine,
    const bool scale
)
{
    this->errorNoScaling(scale);

    const labelVector R(coarse.mshPart().R());
    const label d = coarse.directionNum();

    coarse = Zero;

    // Average face values of corresponding fine grid faces

    forAllFaces(coarse, fd, i, j, k)
    {
        labelVector ijk(i,j,k);
        labelVector nei(ijk-units[fd]);

        labelVector fijk(i*R.x(), j*R.y(), k*R.z());

        labelVector o(Zero);

        if (fd == d)
            o[d] = -1;

        label nFaces = 0;
        for (o.x(); o.x() < (fd == 0 || d == 0 ? 1 : R.x()); o.x()++)
        for (o.y(); o.y() < (fd == 1 || d == 1 ? 1 : R.y()); o.y()++)
        for (o.z(); o.z() < (fd == 2 || d == 2 ? 1 : R.z()); o.z()++)
        {
            coarse(ijk)[fd] += fine(fijk+o)[fd];
            nFaces++;
        }

        coarse(ijk)[fd] /= scalar(nFaces);
    }
}

}

}

}
