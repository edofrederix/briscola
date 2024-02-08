#include "coloFaceAverageRestrictionScheme.H"

#include "colocated.H"
#include "staggered.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type>
coloFaceAverageRestrictionScheme<Type>::coloFaceAverageRestrictionScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    restrictionScheme<Type,colocated>(dict,fvMsh)
{}

template<class Type>
coloFaceAverageRestrictionScheme<Type>::coloFaceAverageRestrictionScheme
(
    const fvMesh& fvMsh
)
:
    restrictionScheme<Type,colocated>(dictionary(),fvMsh)
{}

template<class Type>
void coloFaceAverageRestrictionScheme<Type>::restrict
(
    meshDirection<Type,colocated>& coarse,
    const meshDirection<Type,colocated>& fine,
    const bool scale
)
{
    this->errorNoScaling(scale);

    const labelVector R(coarse.mshPart().R());

    coarse = Zero;

    // Average face values of corresponding fine grid faces

    forAllFaces(coarse, fd, i, j, k)
    {
        labelVector ijk(i,j,k);
        labelVector nei(ijk-units[fd]);

        labelVector fijk(i*R.x(), j*R.y(), k*R.z());
        labelVector fnei(fijk-units[fd]);

        label nFaces = R[(fd+1)%3]*R[(fd+2)%3];

        labelVector o;
        for (o.x() = 0; o.x() < (fd == 0 ? 1 : R.x()); o.x()++)
        for (o.y() = 0; o.y() < (fd == 1 ? 1 : R.y()); o.y()++)
        for (o.z() = 0; o.z() < (fd == 2 ? 1 : R.z()); o.z()++)
        {
            coarse(ijk)[fd] += fine(fijk+o)[fd]/nFaces;
        }
    }
}

}

}

}
