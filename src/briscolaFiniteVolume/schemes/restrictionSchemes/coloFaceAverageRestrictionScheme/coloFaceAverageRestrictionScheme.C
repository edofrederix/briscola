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
    const fvMesh& fvMsh,
    Istream& is
)
:
    restrictionScheme<FaceSpace<Type>,colocated>(fvMsh, is)
{}

template<class Type>
void coloFaceAverageRestrictionScheme<Type>::restrict
(
    meshDirection<FaceSpace<Type>,colocated>& coarse,
    const meshDirection<FaceSpace<Type>,colocated>& fine,
    const bool scale
)
{
    this->errorNoScaling(scale);

    const labelVector R(coarse.mshPart().R());

    // Average face values of corresponding fine grid faces

    forAllFaces(coarse, fd, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        labelVector fijk(briscola::cmptMultiply(ijk,R));
        labelVector fnei(briscola::cmptMultiply(nei,R));

        label nFaces = R[(fd+1)%3]*R[(fd+2)%3];

        coarse(ijk)[fd*2  ] = Zero;
        coarse(nei)[fd*2+1] = Zero;

        labelVector o;
        for (o.x() = 0; o.x() < (fd == 0 ? 1 : R.x()); o.x()++)
        for (o.y() = 0; o.y() < (fd == 1 ? 1 : R.y()); o.y()++)
        for (o.z() = 0; o.z() < (fd == 2 ? 1 : R.z()); o.z()++)
        {
            coarse(ijk)[fd*2  ] += fine(fijk+o)[fd*2  ]/nFaces;
            coarse(nei)[fd*2+1] += fine(fnei+o)[fd*2+1]/nFaces;
        }
    }
}

}

}

}
