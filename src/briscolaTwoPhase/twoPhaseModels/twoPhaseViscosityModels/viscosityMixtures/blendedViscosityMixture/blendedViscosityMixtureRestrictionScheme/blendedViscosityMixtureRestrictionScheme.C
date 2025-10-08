#include "blendedViscosityMixtureRestrictionScheme.H"

#include "colocated.H"
#include "staggered.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Colocated

void blendedViscosityMixtureRestrictionScheme<colocated>::restrict
(
    meshDirection<faceScalar,colocated>& coarse,
    const meshDirection<faceScalar,colocated>& fine,
    const bool scale
)
{
    this->errorNoScaling(scale);

    const labelVector R(coarse.mshPart().R());

    const meshDirection<faceScalar,colocated>& faf =
        this->fvMsh().template
        metrics<colocated>().faceAreas()
        [fine.levelNum()][fine.directionNum()];

    forAllFaces(coarse, fd, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        const labelVector fijk(briscola::cmptMultiply(ijk,R));

        // Reconstruct the face alpha

        scalar alpha = 0.0;
        scalar area = 0.0;

        labelVector o;
        for (o.x() = 0; o.x() < (fd == 0 ? 1 : R.x()); o.x()++)
        for (o.y() = 0; o.y() < (fd == 1 ? 1 : R.y()); o.y()++)
        for (o.z() = 0; o.z() < (fd == 2 ? 1 : R.z()); o.z()++)
        {
            alpha += faf(fijk+o)[fd*2]*this->inv(fine(fijk+o)[fd*2]);
            area +=  faf(fijk+o)[fd*2];
        }

        alpha /= area;

        // Apply the blending function. Viscosity has no sign so we can copy
        // directly to the upper face of the lower neighbor.

        coarse(ijk)[fd*2  ] = this->blend(alpha);
        coarse(nei)[fd*2+1] = coarse(ijk)[fd*2];
    }
}

// Staggered

void blendedViscosityMixtureRestrictionScheme<staggered>::restrict
(
    meshDirection<faceScalar,staggered>& coarse,
    const meshDirection<faceScalar,staggered>& fine,
    const bool scale
)
{
    this->errorNoScaling(scale);

    const label l = fine.levelNum();

    // First interpolate to colocated faces

    meshDirection<faceScalar,colocated> coloFine(fine.fvMsh(), l, 0);
    const meshLevel<faceScalar,staggered>& finel = fine.mshLevel();

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

    // Blended viscosity calculation from the reconstructed volume fractions at
    // the colocated faces. We must use a cell iterator in order to avoid
    // accessing non-existent ghost cells of colocated ghost cells.

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
            scalar alpha = 0.0;

            labelVector o;
            for (o.x() = 0; o.x() < (fd == 0 ? 1 : R.x()); o.x()++)
            for (o.y() = 0; o.y() < (fd == 1 ? 1 : R.y()); o.y()++)
            for (o.z() = 0; o.z() < (fd == 2 ? 1 : R.z()); o.z()++)
            {
                alpha += faf(fijk+o+s)[f]*this->inv(coloFine(fijk+o+s)[f]);
                area += faf(fijk+o+s)[f];
            }

            alpha /= area;

            // Apply the blending function from the reconstructed alpha

            coarse(ijk)[f] = this->blend(alpha);
        }
    }
}

}

}

}
