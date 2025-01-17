#include "coloHarmonicFaceAreaWeightedRestrictionScheme.H"

#include "colocated.H"
#include "staggered.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type>
coloHarmonicFaceAreaWeightedRestrictionScheme<Type>::
coloHarmonicFaceAreaWeightedRestrictionScheme
(
    const fvMesh& fvMsh,
    Istream& is
)
:
    restrictionScheme<Type,colocated>(fvMsh, is)
{}

template<class Type>
void coloHarmonicFaceAreaWeightedRestrictionScheme<Type>::restrict
(
    meshDirection<Type,colocated>& coarse,
    const meshDirection<Type,colocated>& fine,
    const bool scale
)
{
    this->errorNoScaling(scale);

    const labelVector R(coarse.mshPart().R());

    const meshDirection<faceScalar,colocated>& faf =
        this->fvMsh().template
        metrics<colocated>().faceAreas()
        [fine.levelNum()][fine.directionNum()];

    // Harmonic face area weighted average of corresponding fine grid faces

    forAllFaces(coarse, fd, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        const labelVector fijk(briscola::cmptMultiply(ijk,R));
        const labelVector fnei(lowerNeighbor(fijk,fd));

        scalar area = 0.0;

        coarse(ijk)[fd*2  ] = Zero;
        coarse(nei)[fd*2+1] = Zero;

        labelVector o;
        for (o.x() = 0; o.x() < (fd == 0 ? 1 : R.x()); o.x()++)
        for (o.y() = 0; o.y() < (fd == 1 ? 1 : R.y()); o.y()++)
        for (o.z() = 0; o.z() < (fd == 2 ? 1 : R.z()); o.z()++)
        {
            coarse(ijk)[fd*2  ] += faf(fijk+o)[fd*2]/fine(fijk+o)[fd*2  ];
            coarse(nei)[fd*2+1] += faf(fijk+o)[fd*2]/fine(fnei+o)[fd*2+1];

            area += faf(fijk+o)[fd*2];
        }

        coarse(ijk)[fd*2  ] = area/coarse(ijk)[fd*2  ];
        coarse(nei)[fd*2+1] = area/coarse(nei)[fd*2+1];
    }
}

}

}

}
