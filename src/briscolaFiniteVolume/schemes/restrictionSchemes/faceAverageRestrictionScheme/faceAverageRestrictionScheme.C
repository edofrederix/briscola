#include "faceAverageRestrictionScheme.H"

#include "colocated.H"
#include "staggered.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Colocated

template<class Type>
void faceAverageRestrictionScheme<Type,colocated>::restrict
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

// Staggered

template<class Type>
void faceAverageRestrictionScheme<Type,staggered>::restrict
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
