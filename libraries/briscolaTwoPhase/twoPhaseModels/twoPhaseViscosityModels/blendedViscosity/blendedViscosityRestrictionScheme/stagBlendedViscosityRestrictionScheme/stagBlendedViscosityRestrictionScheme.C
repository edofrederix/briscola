#include "stagBlendedViscosityRestrictionScheme.H"

#include "staggered.H"
#include "restrictionSchemes.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeRestrictionSchemeNoTemplate(stagBlendedViscosity,lowerFaceScalar,staggered);

defineTemplateTypeNameAndDebugWithName
(
    blendedViscosityRestrictionScheme<staggered>,
    "blendedViscosity",
    0
);

stagBlendedViscosityRestrictionScheme::stagBlendedViscosityRestrictionScheme
(
    const fvMesh& fvMsh,
    Istream& is
)
:
    blendedViscosityRestrictionScheme<staggered>(fvMsh, is)
{}

void stagBlendedViscosityRestrictionScheme::restrict
(
    meshDirection<lowerFaceScalar,staggered>& coarse,
    const meshDirection<lowerFaceScalar,staggered>& fine,
    const bool scale
)
{
    this->errorNoScaling(scale);

    const labelVector R(coarse.mshPart().R());
    const label d = coarse.directionNum();

    coarse = Zero;

    forAllFaces(coarse, fd, i, j, k)
    {
        labelVector ijk(i,j,k);
        labelVector fijk(i*R.x(), j*R.y(), k*R.z());

        labelVector o(Zero);

        if (fd == d)
            o[d] = -1;

        // Reconstruct the face alpha

        scalar alphaf = 0;
        label nFaces = 0;

        for (o.x(); o.x() < (fd == 0 || d == 0 ? 1 : R.x()); o.x()++)
        for (o.y(); o.y() < (fd == 1 || d == 1 ? 1 : R.y()); o.y()++)
        for (o.z(); o.z() < (fd == 2 || d == 2 ? 1 : R.z()); o.z()++)
        {
            // Don't use cells with negative indices
            labelVector fijko(briscola::cmptMax(fijk+o,zeroXYZ));

            alphaf += this->inv(fine(fijko)[fd]);
            nFaces++;
        }

        alphaf /= nFaces;

        // Apply the blending function

        coarse(ijk)[fd] = this->blend(alphaf);
    }
}

}

}

}
