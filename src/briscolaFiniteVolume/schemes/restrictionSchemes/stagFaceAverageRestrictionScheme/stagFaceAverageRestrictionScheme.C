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
    const fvMesh& fvMsh,
    Istream& is
)
:
    restrictionScheme<FaceSpace<Type>,staggered>(fvMsh, is)
{}

template<class Type>
void stagFaceAverageRestrictionScheme<Type>::restrict
(
    meshDirection<FaceSpace<Type>,staggered>& coarse,
    const meshDirection<FaceSpace<Type>,staggered>& fine,
    const bool scale
)
{
    this->errorNoScaling(scale);

    const labelVector R(coarse.mshPart().R());
    const label d = coarse.directionNum();

    // Average face values of corresponding fine grid faces

    forAllFaces(coarse, fd, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        const labelVector fijk(briscola::cmptMultiply(ijk,R));
        const labelVector fnei(lowerNeighbor(fijk,fd));

        labelVector o(Zero);

        if (fd == d)
            o[d] = -1;

        coarse(ijk)[fd*2  ] = Zero;
        coarse(nei)[fd*2+1] = Zero;

        label nFaces = 0;
        for (o.x(); o.x() < (fd == 0 || d == 0 ? 1 : R.x()); o.x()++)
        for (o.y(); o.y() < (fd == 1 || d == 1 ? 1 : R.y()); o.y()++)
        for (o.z(); o.z() < (fd == 2 || d == 2 ? 1 : R.z()); o.z()++)
        {
            // Don't use cells with negative indices
            const labelVector fijko(briscola::cmptMax(fijk+o,zeroXYZ));
            const labelVector fneio(briscola::cmptMax(fnei+o,zeroXYZ));

            coarse(ijk)[fd*2  ] += fine(fijko)[fd*2 ];
            coarse(nei)[fd*2+1] += fine(fneio)[fd*2+1];

            nFaces++;
        }

        coarse(ijk)[fd*2  ] /= scalar(nFaces);
        coarse(nei)[fd*2+1] /= scalar(nFaces);
    }
}

}

}

}
