#include "stagHarmonicFaceAreaWeightedRestrictionScheme.H"

#include "staggered.H"
#include "staggered.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type>
stagHarmonicFaceAreaWeightedRestrictionScheme<Type>::
stagHarmonicFaceAreaWeightedRestrictionScheme
(
    const fvMesh& fvMsh,
    Istream& is
)
:
    restrictionScheme<Type,staggered>(fvMsh, is)
{}

template<class Type>
void stagHarmonicFaceAreaWeightedRestrictionScheme<Type>::restrict
(
    meshDirection<Type,staggered>& coarse,
    const meshDirection<Type,staggered>& fine,
    const bool scale
)
{
    this->errorNoScaling(scale);

    const label l = fine.levelNum();

    // First interpolate to colocated faces

    meshDirection<Type,colocated> coloFine(fine.fvMsh(), l, 0);
    const meshLevel<Type,staggered>& finel = fine.mshLevel();

    forAllFacesInDirection(coloFine, fd, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(ijk,fd));

        coloFine(ijk)[fd*2] =
            0.5*(finel(fd,ijk)[fd*2] + finel(fd,ijk)[fd*2+1]);

        coloFine(nei)[fd*2+1] = coloFine(ijk)[fd*2];
    }

    // Homogeneous Neumann for required colocated ghost cell faces

    for (int f = 0; f < 6; f++)
    {
        const labelVector bo(faceOffsets[f]);

        const labelVector S(fine.fvMsh().template S<colocated>(l,0,bo));
        const labelVector E(fine.fvMsh().template E<colocated>(l,0,bo));

        labelVector ijk;

        for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
        for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
        for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
        {
            coloFine(ijk+bo)[f] = coloFine(ijk)[(f%2) ? f-1 : f+1];
        }
    }

    // Harmonic face area weighted average of corresponding colocated fine grid
    // faces. We must use a cell iterator in order to avoid accessing
    // non-existent ghost cells of colocated ghost cells.

    const label d = fine.directionNum();

    const labelVector R(coarse.mshPart().R());

    const meshDirection<faceScalar,colocated>& faf =
        this->fvMsh().template metrics<colocated>().faceAreas()[l][0];

    forAllCellsInDirection(coarse, i, j, k)
    {
        const labelVector ijk(i,j,k);

        // The lower colocated fine grid cell contained by the staggered coarse
        // grid cell has an index scaled by R and padding subtracted

        const labelVector fijk
        (
            briscola::cmptMultiply(ijk,R)
          - staggered::padding[d]
        );

        for (int f = 0; f < 6; f++)
        {
            const int fd = f/2;
            const labelVector s((f%2) ? units[fd] : zeroXYZ);

            scalar area = 0.0;
            coarse(ijk)[f] = Zero;

            labelVector o;
            for (o.x() = 0; o.x() < (fd == 0 ? 1 : R.x()); o.x()++)
            for (o.y() = 0; o.y() < (fd == 1 ? 1 : R.y()); o.y()++)
            for (o.z() = 0; o.z() < (fd == 2 ? 1 : R.z()); o.z()++)
            {
                coarse(ijk)[f] +=
                    faf(fijk+o+s)[f]/stabilise(coloFine(fijk+o+s)[f], 1e-12);

                area += faf(fijk+o+s)[f];
            }

            coarse(ijk)[f] = area/coarse(ijk)[f];
        }
    }
}

}

}

}
